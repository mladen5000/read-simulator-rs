# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a genomic FASTQ read simulator written in Rust. It generates synthetic sequencing reads from FASTA reference files, supporting multiple sequencing technologies (Illumina and Nanopore) with realistic error profiles, GC bias, and quality scores.

## Build and Run Commands

```bash
# Build the project
cargo build

# Build optimized release version
cargo build --release

# Run the simulator (use --help for options)
cargo run -- --help

# Example: Generate Illumina paired-end reads with 30X coverage
cargo run -- -i input.fasta -o output_prefix -t illumina -c 30.0

# Example: Generate Nanopore single-end reads
cargo run -- -i input.fasta -o output_prefix -t nanopore -c 20.0 --single-end

# Run tests (currently no tests defined)
cargo test

# Run clippy for linting
cargo clippy

# Format code
cargo fmt
```

## Architecture

### Core Components

The codebase is organized around a single-file architecture in `src/main.rs` with these key structures:

1. **SimpleRng** - Custom LCG-based RNG implementation
   - Provides deterministic seeding for reproducibility
   - Implements various distributions (uniform, normal via Box-Muller, geometric)
   - Avoids external RNG dependencies for portability

2. **SeqTech & TechParams** - Technology profiles
   - `SeqTech`: Enum defining supported sequencing technologies
   - `TechParams`: Technology-specific parameters (read length, error rates, GC bias)
   - Each technology has preset parameters matching real-world characteristics

3. **ReadSimulator** - Main simulation engine
   - Generates reads with technology-specific error profiles
   - Implements GC bias through rejection sampling
   - Handles both paired-end and single-end read generation
   - Key methods:
     - `simulate_from_sequence()`: Main entry point for generating reads from a sequence
     - `introduce_errors()`: Injects insertions, deletions, and substitutions
     - `generate_quality_scores()`: Creates realistic Phred quality scores

4. **ReadData** - Output representation
   - Enum distinguishing paired-end vs single-end reads
   - Stores read name, sequence, and quality scores

### Data Flow

```
FASTA Input → parse_fasta() → ReadSimulator → ReadData collection → write_fastq() → FASTQ Output
```

1. Parse input FASTA file into (name, sequence) tuples
2. For each sequence, calculate number of reads based on coverage depth
3. Sample read positions (with optional GC bias)
4. Generate reads with error injection and quality scores
5. Write reads to FASTQ files (separate R1/R2 files for paired-end)

### Important Implementation Details

- **Coverage calculation**: `num_reads = (seq_length × coverage) / read_length`
- **Error injection**: Uses technology-specific error rates with weighted random selection among insertion/deletion/substitution
- **GC bias**: Rejection sampling with exponential probability based on deviation from 50% GC
- **Paired-end reads**: Simulates fragment insert size with normal distribution (mean=500bp, std=50bp)
- **Nanopore reads**: Variable length reads using geometric distribution
- **Quality scores**:
  - Illumina: High quality (Q38) with 3' degradation
  - Nanopore: Lower quality (Q15 mean) with higher variance

## Dependencies

- **clap 4.5.50** with `derive` feature: Command-line argument parsing using derive macros

## Project Configuration

- Rust edition: 2021
- No external RNG dependencies (uses custom SimpleRng for portability and determinism)
