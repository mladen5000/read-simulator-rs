use clap::{Parser, ValueEnum};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

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
}

struct ReadSimulator {
    tech: TechParams,
    coverage: f64,
    rng: SimpleRng,
    bases: Vec<u8>,
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
            bases: vec![b'A', b'C', b'G', b'T'],
        }
    }

    fn complement(&self, base: u8) -> u8 {
        match base {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            _ => b'N',
        }
    }

    fn reverse_complement(&self, seq: &[u8]) -> Vec<u8> {
        seq.iter().rev().map(|&b| self.complement(b)).collect()
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

    fn introduce_errors(&mut self, seq: &[u8]) -> Vec<u8> {
        let mut result = Vec::with_capacity(seq.len());
        let total_error =
            self.tech.insertion_rate + self.tech.deletion_rate + self.tech.substitution_rate;
        let ins_threshold = self.tech.insertion_rate / total_error;
        let del_threshold = ins_threshold + self.tech.deletion_rate / total_error;

        for &base in seq {
            if self.rng.gen_f64() < self.tech.error_rate {
                let error_type = self.rng.gen_f64();

                if error_type < ins_threshold {
                    // Insertion
                    result.push(base);
                    result.push(self.bases[self.rng.gen_range(0, 3)]);
                } else if error_type < del_threshold {
                    // Deletion - skip base
                } else {
                    // Substitution
                    let alternatives: Vec<u8> =
                        self.bases.iter().filter(|&&b| b != base).copied().collect();
                    if !alternatives.is_empty() {
                        result.push(alternatives[self.rng.gen_range(0, alternatives.len() - 1)]);
                    }
                }
            } else {
                result.push(base);
            }
        }

        result
    }

    fn generate_quality_scores(&mut self, length: usize) -> Vec<u8> {
        let mut scores = Vec::with_capacity(length);

        if self.tech.name == "Illumina" {
            let base_qual = 38i32;
            for i in 0..length {
                let degradation = (3.0 * (i as f64 / length as f64)) as i32;
                let noise = (self.rng.gen_f64() * 4.0 - 2.0) as i32;
                let qual = (base_qual - degradation + noise).clamp(20, 41);
                scores.push((qual + 33) as u8);
            }
        } else {
            // Nanopore
            for _ in 0..length {
                let qual = self.rng.box_muller(15.0, 5.0).round() as i32;
                let qual = qual.clamp(5, 30);
                scores.push((qual + 33) as u8);
            }
        }

        scores
    }

    fn sample_read_position(&mut self, seq: &[u8], read_length: usize) -> (usize, bool) {
        let max_start = if seq.len() > read_length {
            seq.len() - read_length
        } else {
            return (0, self.rng.gen_bool());
        };

        // GC bias rejection sampling
        if self.tech.gc_bias_strength > 0.0 {
            for _ in 0..10 {
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

        let mut fragment = seq[start_pos..start_pos + insert_size].to_vec();

        if is_reverse {
            fragment = self.reverse_complement(&fragment);
        }

        // R1: forward from fragment start
        let read1_seq = self.introduce_errors(&fragment[..read_len.min(fragment.len())]);
        let qual1 = self.generate_quality_scores(read1_seq.len());

        // R2: reverse complement from fragment end
        let r2_start = fragment.len().saturating_sub(read_len);
        let read2_seq = self.reverse_complement(&fragment[r2_start..]);
        let read2_seq = self.introduce_errors(&read2_seq);
        let qual2 = self.generate_quality_scores(read2_seq.len());

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

        // Variable length for Nanopore
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
        let mut read_seq = seq[start_pos..end_pos].to_vec();

        if is_reverse {
            read_seq = self.reverse_complement(&read_seq);
        }

        let read_seq = self.introduce_errors(&read_seq);
        let qual = self.generate_quality_scores(read_seq.len());
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

fn write_fastq(
    reads: &[ReadData],
    output_prefix: &str,
    paired_end: bool,
    _gzip_output: bool,
) -> io::Result<()> {
    if paired_end {
        let r1_file = format!("{}_R1.fastq", output_prefix);
        let r2_file = format!("{}_R2.fastq", output_prefix);

        let mut r1_writer = BufWriter::new(File::create(&r1_file)?);
        let mut r2_writer = BufWriter::new(File::create(&r2_file)?);

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
                writeln!(
                    r1_writer,
                    "{}\n{}\n+\n{}",
                    name1,
                    String::from_utf8_lossy(seq1),
                    String::from_utf8_lossy(qual1)
                )?;
                writeln!(
                    r2_writer,
                    "{}\n{}\n+\n{}",
                    name2,
                    String::from_utf8_lossy(seq2),
                    String::from_utf8_lossy(qual2)
                )?;
            }
        }

        println!("Written paired-end reads to {} and {}", r1_file, r2_file);
    } else {
        let out_file = format!("{}.fastq", output_prefix);
        let mut writer = BufWriter::new(File::create(&out_file)?);

        for read in reads {
            if let ReadData::SingleEnd { name, seq, qual } = read {
                writeln!(
                    writer,
                    "{}\n{}\n+\n{}",
                    name,
                    String::from_utf8_lossy(seq),
                    String::from_utf8_lossy(qual)
                )?;
            }
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

    /// Do not gzip output files
    #[arg(long)]
    no_gzip: bool,

    /// Random seed for reproducibility
    #[arg(long)]
    seed: Option<u64>,
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();

    let tech = match cli.tech {
        SeqTech::Illumina => TechParams::illumina(!cli.single_end),
        SeqTech::Nanopore => TechParams::nanopore(),
    };

    println!("Simulating {} reads...", tech.name);
    println!("  Read length: {} bp", tech.read_length);
    println!("  Paired-end: {}", tech.paired_end);
    println!("  Error rate: {:.4}", tech.error_rate);
    println!("  Target coverage: {:.1}X", cli.coverage);

    let mut simulator = ReadSimulator::new(tech.clone(), cli.coverage, cli.seed);

    println!("\nReading sequences from {:?}...", cli.input);
    let sequences = parse_fasta(&cli.input)?;
    println!("Found {} sequence(s)", sequences.len());

    let mut all_reads = Vec::new();

    for (seq_name, seq) in sequences {
        println!("  Generating reads for {} ({} bp)...", seq_name, seq.len());
        let reads = simulator.simulate_from_sequence(&seq, &seq_name);
        all_reads.extend(reads);
    }

    println!(
        "\nGenerated {} read{}(s)",
        all_reads.len(),
        if tech.paired_end { " pair" } else { "" }
    );

    write_fastq(&all_reads, &cli.output, tech.paired_end, !cli.no_gzip)?;
    println!("Done!");

    Ok(())
}
