use clap::{Parser, ValueEnum};
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use rayon::prelude::*;

// Simple LCG-based RNG for deterministic seeding
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn from_entropy() -> Self {
        use std::time::{SystemTime, UNIX_EPOCH};
        let seed = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos() as u64;
        Self::new(seed)
    }

    fn next_u64(&mut self) -> u64 {
        // LCG constants from Numerical Recipes
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1);
        self.state
    }

    fn gen_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn gen_range(&mut self, min: usize, max: usize) -> usize {
        min + (self.next_u64() as usize % (max - min + 1))
    }

    fn gen_bool(&mut self) -> bool {
        self.next_u64() & 1 == 1
    }

    // Box-Muller transform for normal distribution
    fn box_muller(&mut self, mean: f64, std_dev: f64) -> f64 {
        let u1 = self.gen_f64();
        let u2 = self.gen_f64();
        let z0 = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
        mean + z0 * std_dev
    }

    // Geometric distribution
    fn geometric(&mut self, p: f64) -> usize {
        let u = self.gen_f64();
        ((u.ln() / (1.0 - p).ln()).floor() as usize).max(1)
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum SeqTech {
    Illumina,
    Nanopore,
    /// PacBio HiFi (Circular Consensus Sequencing)
    #[value(name = "pacbio")]
    PacBio,
    /// Illumina NovaSeq 6000
    #[value(name = "novaseq")]
    NovaSeq,
    /// Oxford Nanopore MinION R10.4
    #[value(name = "minion")]
    MinION,
}

#[derive(Debug, Clone)]
struct TechParams {
    name: String,
    read_length: usize,
    paired_end: bool,
    error_rate: f64,
    insertion_rate: f64,
    deletion_rate: f64,
    substitution_rate: f64,
    gc_bias_strength: f64,
}

impl TechParams {
    fn illumina(paired_end: bool) -> Self {
        Self {
            name: "Illumina".to_string(),
            read_length: 150,
            paired_end,
            error_rate: 0.001,
            insertion_rate: 0.0001,
            deletion_rate: 0.0001,
            substitution_rate: 0.0008,
            gc_bias_strength: 0.3,
        }
    }

    fn nanopore() -> Self {
        Self {
            name: "Nanopore".to_string(),
            read_length: 10000,
            paired_end: false,
            error_rate: 0.05,
            insertion_rate: 0.02,
            deletion_rate: 0.02,
            substitution_rate: 0.01,
            gc_bias_strength: 0.15,
        }
    }

    // Quick Win: Preset profiles for common platforms
    fn pacbio() -> Self {
        Self {
            name: "PacBio".to_string(),
            read_length: 15000,
            paired_end: false,
            error_rate: 0.002, // HiFi has very low error rate
            insertion_rate: 0.0008,
            deletion_rate: 0.0008,
            substitution_rate: 0.0004,
            gc_bias_strength: 0.1,
        }
    }

    fn novaseq() -> Self {
        Self {
            name: "NovaSeq".to_string(),
            read_length: 150,
            paired_end: true,
            error_rate: 0.0008, // Slightly better than standard Illumina
            insertion_rate: 0.00008,
            deletion_rate: 0.00008,
            substitution_rate: 0.00064,
            gc_bias_strength: 0.25,
        }
    }

    fn minion() -> Self {
        Self {
            name: "MinION".to_string(),
            read_length: 12000,
            paired_end: false,
            error_rate: 0.04, // R10.4 chemistry
            insertion_rate: 0.015,
            deletion_rate: 0.015,
            substitution_rate: 0.01,
            gc_bias_strength: 0.12,
        }
    }
}

// Phase 1 Optimization: Pre-computed lookup tables and buffers
struct SimulatorBuffers {
    complement_table: [u8; 256],
    bases: [u8; 4],
    alternatives: [[u8; 3]; 4],
    base_to_idx: [u8; 256],  // Quick Win: Base to index lookup
    // Phase 3 Optimization: Quality score pools for fast generation
    illumina_qual_pools: Vec<Vec<u8>>,
    nanopore_qual_pool: Vec<u8>,
}

impl SimulatorBuffers {
    fn new() -> Self {
        let mut complement_table = [0u8; 256];
        complement_table[b'A' as usize] = b'T';
        complement_table[b'T' as usize] = b'A';
        complement_table[b'G' as usize] = b'C';
        complement_table[b'C' as usize] = b'G';
        for i in 0..256 {
            if complement_table[i] == 0 {
                complement_table[i] = b'N';
            }
        }

        let bases = [b'A', b'C', b'G', b'T'];

        let alternatives = [
            [b'C', b'G', b'T'],
            [b'A', b'G', b'T'],
            [b'A', b'C', b'T'],
            [b'A', b'C', b'G'],
        ];

        // Quick Win: Base to index lookup table
        let mut base_to_idx = [0u8; 256];
        base_to_idx[b'A' as usize] = 0;
        base_to_idx[b'C' as usize] = 1;
        base_to_idx[b'G' as usize] = 2;
        base_to_idx[b'T' as usize] = 3;

        // Pre-compute quality score pools (power-of-2 size for fast bit masking)
        let mut illumina_qual_pools = Vec::new();
        let mut rng = SimpleRng::new(42);
        const POOL_SIZE: usize = 16384; // 2^14 for fast bit masking

        for position_ratio in [0.0f64, 0.3, 0.6, 0.9] {
            let mut pool = Vec::with_capacity(POOL_SIZE);
            for _ in 0..POOL_SIZE {
                let base_qual = 38i32;
                let degradation = (3.0 * position_ratio) as i32;
                let noise = (rng.gen_f64() * 4.0 - 2.0) as i32;
                let qual = (base_qual - degradation + noise).clamp(20, 41);
                pool.push((qual + 33) as u8);
            }
            illumina_qual_pools.push(pool);
        }

        let mut nanopore_qual_pool = Vec::with_capacity(POOL_SIZE);
        for _ in 0..POOL_SIZE {
            let qual = rng.box_muller(15.0, 5.0).round() as i32;
            let qual = qual.clamp(5, 30);
            nanopore_qual_pool.push((qual + 33) as u8);
        }

        Self {
            complement_table,
            bases,
            alternatives,
            base_to_idx,
            illumina_qual_pools,
            nanopore_qual_pool,
        }
    }
}

struct ReadSimulator {
    tech: TechParams,
    coverage: f64,
    rng: SimpleRng,
    buffers: SimulatorBuffers,
}

impl ReadSimulator {
    fn new(tech: TechParams, coverage: f64, seed: Option<u64>) -> Self {
        let rng = seed
            .map(SimpleRng::new)
            .unwrap_or_else(SimpleRng::from_entropy);

        Self {
            tech,
            coverage,
            rng,
            buffers: SimulatorBuffers::new(),
        }
    }

    #[inline]
    fn reverse_complement_in_place(&self, seq: &[u8], out: &mut Vec<u8>) {
        out.clear();
        out.reserve(seq.len());
        for &b in seq.iter().rev() {
            out.push(self.buffers.complement_table[b as usize]);
        }
    }

    fn calculate_gc_content(&self, seq: &[u8]) -> f64 {
        if seq.is_empty() {
            return 0.5;
        }
        let gc_count = seq.iter().filter(|&&b| b == b'G' || b == b'C').count();
        gc_count as f64 / seq.len() as f64
    }

    fn gc_bias_probability(&self, gc_content: f64) -> f64 {
        let bias = self.tech.gc_bias_strength;
        let deviation = (gc_content - 0.5).abs();
        (-bias * 10.0 * deviation.powi(2)).exp()
    }

    #[inline]
    fn introduce_errors_in_place(&mut self, seq: &[u8], out: &mut Vec<u8>) {
        out.clear();
        out.reserve(seq.len());

        let total_error =
            self.tech.insertion_rate + self.tech.deletion_rate + self.tech.substitution_rate;
        let ins_threshold = self.tech.insertion_rate / total_error;
        let del_threshold = ins_threshold + self.tech.deletion_rate / total_error;

        for &base in seq {
            if self.rng.gen_f64() < self.tech.error_rate {
                let error_type = self.rng.gen_f64();

                if error_type < ins_threshold {
                    out.push(base);
                    out.push(self.buffers.bases[self.rng.gen_range(0, 3)]);
                } else if error_type < del_threshold {
                    // Deletion - skip base
                } else {
                    // Substitution with pre-computed alternatives
                    let base_idx = self.buffers.base_to_idx[base as usize];
                    let alternatives = self.buffers.alternatives[base_idx as usize];
                    out.push(alternatives[self.rng.gen_range(0, 2)]);
                }
            } else {
                out.push(base);
            }
        }
    }

    #[inline]
    fn generate_quality_scores_fast(&mut self, length: usize, out: &mut Vec<u8>) {
        out.clear();
        out.reserve(length);

        if self.tech.name == "Illumina" {
            // Pool-based quality generation for Illumina with bit masking
            const POOL_MASK: u64 = 0x3FFF; // 16384 - 1 = 0x3FFF for fast modulo
            for i in 0..length {
                let position_ratio = (i as f64) / (length as f64);
                let pool_idx = if position_ratio < 0.25 {
                    0
                } else if position_ratio < 0.5 {
                    1
                } else if position_ratio < 0.75 {
                    2
                } else {
                    3
                };
                let pool = &self.buffers.illumina_qual_pools[pool_idx];
                // Fast bit masking instead of modulo
                let pool_idx_rand = (self.rng.next_u64() & POOL_MASK) as usize;
                out.push(pool[pool_idx_rand]);
            }
        } else {
            // Nanopore with pre-generated pool and bit masking
            const POOL_MASK: u64 = 0x3FFF;
            for _ in 0..length {
                let pool_idx = (self.rng.next_u64() & POOL_MASK) as usize;
                out.push(self.buffers.nanopore_qual_pool[pool_idx]);
            }
        }
    }

    fn sample_read_position(&mut self, seq: &[u8], read_length: usize) -> (usize, bool) {
        let max_start = if seq.len() > read_length {
            seq.len() - read_length
        } else {
            return (0, self.rng.gen_bool());
        };

        if self.tech.gc_bias_strength > 0.1 {
            for _ in 0..5 {
                let pos = self.rng.gen_range(0, max_start);
                let end = (pos + read_length).min(seq.len());
                let gc = self.calculate_gc_content(&seq[pos..end]);
                let prob = self.gc_bias_probability(gc);

                if self.rng.gen_f64() < prob {
                    return (pos, self.rng.gen_bool());
                }
            }
        }

        (self.rng.gen_range(0, max_start), self.rng.gen_bool())
    }

    fn generate_paired_end_reads(
        &mut self,
        seq: &[u8],
        seq_name: &str,
        read_id: usize,
    ) -> Option<(String, Vec<u8>, Vec<u8>, String, Vec<u8>, Vec<u8>)> {
        let read_len = self.tech.read_length;
        let insert_size = self.rng.box_muller(500.0, 50.0).round() as usize;
        let insert_size = insert_size.max(read_len * 2 + 50);

        if seq.len() < insert_size {
            return None;
        }

        let max_start = seq.len() - insert_size;
        let (start_pos, is_reverse) = self.sample_read_position(seq, insert_size);
        let start_pos = start_pos.min(max_start);

        let fragment_end = start_pos + insert_size;
        let frag_slice = &seq[start_pos..fragment_end];

        let mut fragment = if is_reverse {
            let mut temp = Vec::with_capacity(frag_slice.len());
            self.reverse_complement_in_place(frag_slice, &mut temp);
            temp
        } else {
            frag_slice.to_vec()
        };

        // R1: forward from fragment start
        let read1_end = read_len.min(fragment.len());
        let mut read1_seq = Vec::with_capacity(read1_end);
        self.introduce_errors_in_place(&fragment[..read1_end], &mut read1_seq);
        let mut qual1 = Vec::with_capacity(read1_seq.len());
        self.generate_quality_scores_fast(read1_seq.len(), &mut qual1);

        // R2: reverse complement from fragment end
        let r2_start = fragment.len().saturating_sub(read_len);
        let mut r2_rev = Vec::with_capacity(fragment.len() - r2_start);
        self.reverse_complement_in_place(&fragment[r2_start..], &mut r2_rev);
        let mut read2_seq = Vec::with_capacity(r2_rev.len());
        self.introduce_errors_in_place(&r2_rev, &mut read2_seq);
        let mut qual2 = Vec::with_capacity(read2_seq.len());
        self.generate_quality_scores_fast(read2_seq.len(), &mut qual2);

        let name1 = format!("@{}.{}/1", seq_name, read_id);
        let name2 = format!("@{}.{}/2", seq_name, read_id);

        Some((name1, read1_seq, qual1, name2, read2_seq, qual2))
    }

    fn generate_single_end_read(
        &mut self,
        seq: &[u8],
        seq_name: &str,
        read_id: usize,
    ) -> Option<(String, Vec<u8>, Vec<u8>)> {
        let mut read_len = self.tech.read_length;

        if self.tech.name == "Nanopore" {
            let length_variation = self.rng.geometric(1.0 / read_len as f64);
            read_len = ((read_len as f64 * 0.5) as usize + length_variation)
                .max(1000)
                .min(seq.len());
        }

        if seq.len() < read_len {
            return None;
        }

        let (start_pos, is_reverse) = self.sample_read_position(seq, read_len);
        let end_pos = (start_pos + read_len).min(seq.len());

        let mut read_data = if is_reverse {
            let mut temp = Vec::with_capacity(end_pos - start_pos);
            self.reverse_complement_in_place(&seq[start_pos..end_pos], &mut temp);
            temp
        } else {
            seq[start_pos..end_pos].to_vec()
        };

        let mut read_seq = Vec::with_capacity(read_data.len());
        self.introduce_errors_in_place(&read_data, &mut read_seq);

        let mut qual = Vec::with_capacity(read_seq.len());
        self.generate_quality_scores_fast(read_seq.len(), &mut qual);

        let name = format!("@{}.{}", seq_name, read_id);

        Some((name, read_seq, qual))
    }

    fn simulate_from_sequence(&mut self, seq: &[u8], seq_name: &str) -> Vec<ReadData> {
        let seq_length = seq.len();
        let read_length = self.tech.read_length;
        let num_reads = ((seq_length as f64 * self.coverage) / read_length as f64) as usize;

        let mut reads = Vec::new();

        for i in 0..num_reads {
            if self.tech.paired_end {
                if let Some((n1, s1, q1, n2, s2, q2)) =
                    self.generate_paired_end_reads(seq, seq_name, i)
                {
                    reads.push(ReadData::PairedEnd {
                        name1: n1,
                        seq1: s1,
                        qual1: q1,
                        name2: n2,
                        seq2: s2,
                        qual2: q2,
                    });
                }
            } else if let Some((name, seq_data, qual)) =
                self.generate_single_end_read(seq, seq_name, i)
            {
                reads.push(ReadData::SingleEnd {
                    name,
                    seq: seq_data,
                    qual,
                });
            }
        }

        reads
    }
}

#[derive(Debug)]
enum ReadData {
    PairedEnd {
        name1: String,
        seq1: Vec<u8>,
        qual1: Vec<u8>,
        name2: String,
        seq2: Vec<u8>,
        qual2: Vec<u8>,
    },
    SingleEnd {
        name: String,
        seq: Vec<u8>,
        qual: Vec<u8>,
    },
}

fn parse_fasta(path: &PathBuf) -> io::Result<Vec<(String, Vec<u8>)>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut sequences = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        if line.starts_with('>') {
            if let Some(name) = current_name.take() {
                sequences.push((name, current_seq.clone()));
                current_seq.clear();
            }
            current_name = Some(
                line[1..]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string(),
            );
        } else {
            current_seq.extend_from_slice(line.to_uppercase().as_bytes());
        }
    }

    if let Some(name) = current_name {
        sequences.push((name, current_seq));
    }

    Ok(sequences)
}

// Phase 1 Optimization: Direct byte writing instead of writeln!
// Quick Win: Gzip compression support
fn write_fastq(
    reads: &[ReadData],
    output_prefix: &str,
    paired_end: bool,
    gzip_output: bool,
) -> io::Result<()> {
    if paired_end {
        let (r1_file, r2_file) = if gzip_output {
            (
                format!("{}_R1.fastq.gz", output_prefix),
                format!("{}_R2.fastq.gz", output_prefix),
            )
        } else {
            (
                format!("{}_R1.fastq", output_prefix),
                format!("{}_R2.fastq", output_prefix),
            )
        };

        if gzip_output {
            let r1_file_handle = File::create(&r1_file)?;
            let r2_file_handle = File::create(&r2_file)?;
            let mut r1_writer = BufWriter::with_capacity(
                1024 * 1024,
                GzEncoder::new(r1_file_handle, Compression::fast()),
            );
            let mut r2_writer = BufWriter::with_capacity(
                1024 * 1024,
                GzEncoder::new(r2_file_handle, Compression::fast()),
            );

            for read in reads {
                if let ReadData::PairedEnd {
                    name1,
                    seq1,
                    qual1,
                    name2,
                    seq2,
                    qual2,
                } = read
                {
                    r1_writer.write_all(name1.as_bytes())?;
                    r1_writer.write_all(b"\n")?;
                    r1_writer.write_all(seq1)?;
                    r1_writer.write_all(b"\n+\n")?;
                    r1_writer.write_all(qual1)?;
                    r1_writer.write_all(b"\n")?;

                    r2_writer.write_all(name2.as_bytes())?;
                    r2_writer.write_all(b"\n")?;
                    r2_writer.write_all(seq2)?;
                    r2_writer.write_all(b"\n+\n")?;
                    r2_writer.write_all(qual2)?;
                    r2_writer.write_all(b"\n")?;
                }
            }

            r1_writer.flush()?;
            r2_writer.flush()?;
        } else {
            let mut r1_writer = BufWriter::with_capacity(1024 * 1024, File::create(&r1_file)?);
            let mut r2_writer = BufWriter::with_capacity(1024 * 1024, File::create(&r2_file)?);

            for read in reads {
                if let ReadData::PairedEnd {
                    name1,
                    seq1,
                    qual1,
                    name2,
                    seq2,
                    qual2,
                } = read
                {
                    r1_writer.write_all(name1.as_bytes())?;
                    r1_writer.write_all(b"\n")?;
                    r1_writer.write_all(seq1)?;
                    r1_writer.write_all(b"\n+\n")?;
                    r1_writer.write_all(qual1)?;
                    r1_writer.write_all(b"\n")?;

                    r2_writer.write_all(name2.as_bytes())?;
                    r2_writer.write_all(b"\n")?;
                    r2_writer.write_all(seq2)?;
                    r2_writer.write_all(b"\n+\n")?;
                    r2_writer.write_all(qual2)?;
                    r2_writer.write_all(b"\n")?;
                }
            }

            r1_writer.flush()?;
            r2_writer.flush()?;
        }

        println!("Written paired-end reads to {} and {}", r1_file, r2_file);
    } else {
        let out_file = if gzip_output {
            format!("{}.fastq.gz", output_prefix)
        } else {
            format!("{}.fastq", output_prefix)
        };

        if gzip_output {
            let file_handle = File::create(&out_file)?;
            let mut writer = BufWriter::with_capacity(
                1024 * 1024,
                GzEncoder::new(file_handle, Compression::fast()),
            );

            for read in reads {
                if let ReadData::SingleEnd { name, seq, qual } = read {
                    writer.write_all(name.as_bytes())?;
                    writer.write_all(b"\n")?;
                    writer.write_all(seq)?;
                    writer.write_all(b"\n+\n")?;
                    writer.write_all(qual)?;
                    writer.write_all(b"\n")?;
                }
            }

            writer.flush()?;
        } else {
            let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(&out_file)?);

            for read in reads {
                if let ReadData::SingleEnd { name, seq, qual } = read {
                    writer.write_all(name.as_bytes())?;
                    writer.write_all(b"\n")?;
                    writer.write_all(seq)?;
                    writer.write_all(b"\n+\n")?;
                    writer.write_all(qual)?;
                    writer.write_all(b"\n")?;
                }
            }

            writer.flush()?;
        }

        println!("Written single-end reads to {}", out_file);
    }

    Ok(())
}

#[derive(Parser, Debug, Clone)]
#[command(name = "fastq_simulator")]
#[command(version = "0.2.0")]
#[command(about = "Simulate FASTQ reads from genomic FASTA files", long_about = None)]
struct Cli {
    /// Input FASTA file
    #[arg(short, long)]
    input: PathBuf,

    /// Output prefix for FASTQ files
    #[arg(short, long)]
    output: String,

    /// Sequencing technology
    #[arg(short, long, value_enum, default_value = "illumina")]
    tech: SeqTech,

    /// Target coverage depth
    #[arg(short, long, default_value = "10.0")]
    coverage: f64,

    /// Generate single-end reads (default: paired-end for Illumina)
    #[arg(long)]
    single_end: bool,

    /// Gzip compress output files (slower but saves disk space)
    #[arg(long)]
    gzip: bool,

    /// Random seed for reproducibility
    #[arg(long)]
    seed: Option<u64>,
}

fn main() -> io::Result<()> {
    use std::time::Instant;

    let cli = Cli::parse();

    let tech = match cli.tech {
        SeqTech::Illumina => TechParams::illumina(!cli.single_end),
        SeqTech::Nanopore => TechParams::nanopore(),
        SeqTech::PacBio => TechParams::pacbio(),
        SeqTech::NovaSeq => TechParams::novaseq(),
        SeqTech::MinION => TechParams::minion(),
    };

    println!("Simulating {} reads...", tech.name);
    println!("  Read length: {} bp", tech.read_length);
    println!("  Paired-end: {}", tech.paired_end);
    println!("  Error rate: {:.4}", tech.error_rate);
    println!("  Target coverage: {:.1}X", cli.coverage);

    println!("\nReading sequences from {:?}...", cli.input);
    let start = Instant::now();
    let sequences = parse_fasta(&cli.input)?;
    eprintln!("[TIMING] FASTA parsing: {:?}", start.elapsed());
    println!("Found {} sequence(s)", sequences.len());

    // Quick Win: Conditional parallelization - only parallelize for multi-sequence files
    let start = Instant::now();
    let all_reads: Vec<ReadData> = if sequences.len() == 1 {
        // Single sequence: sequential processing (avoid Rayon overhead)
        sequences
            .into_iter()
            .flat_map(|(seq_name, seq)| {
                eprintln!("  Generating reads for {} ({} bp)...", seq_name, seq.len());
                let mut sim = ReadSimulator::new(tech.clone(), cli.coverage, cli.seed);
                sim.simulate_from_sequence(&seq, &seq_name)
            })
            .collect()
    } else {
        // Multiple sequences: parallel processing with Rayon
        sequences
            .into_par_iter()
            .flat_map(|(seq_name, seq)| {
                eprintln!("  Generating reads for {} ({} bp)...", seq_name, seq.len());
                let mut sim = ReadSimulator::new(tech.clone(), cli.coverage, cli.seed);
                sim.simulate_from_sequence(&seq, &seq_name)
            })
            .collect()
    };
    eprintln!("[TIMING] Read generation: {:?}", start.elapsed());

    println!(
        "\nGenerated {} read{}(s)",
        all_reads.len(),
        if tech.paired_end { " pair" } else { "" }
    );

    let start = Instant::now();
    write_fastq(&all_reads, &cli.output, tech.paired_end, cli.gzip)?;
    eprintln!("[TIMING] FASTQ writing: {:?}", start.elapsed());
    println!("Done!");

    Ok(())
}
