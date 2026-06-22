# any2fasta (Rust Version)
Convert various biological sequence formats to standard FASTA format.
This project is a native Rust rewrite inspired by the original Perl `any2fasta` tool from https://github.com/tseemann/any2fasta.
Eliminates the dependency on the Perl runtime, supports native execution on Windows, Linux and macOS, with better runtime performance and cross-platform compatibility.

## Quick start
```bash
# Build from source
cargo build --release
# Basic usage
./target/release/any2fasta genome.gbk > genome.fasta
# Auto decompress gzip input
./target/release/any2fasta seq.gbk.gz > seq.fasta
# Auto decompress bzip2 input
./target/release/any2fasta protein.pdb.bz2 > protein.fasta
```

## Motivation
The original Perl `any2fasta` is widely used in bioinformatics pipelines for converting diverse sequence formats to standard FASTA. Traditional tools like EMBOSS `seqret` / `readseq` often corrupt sequence IDs containing special characters such as `|` and `.`, leading to inconsistent identifiers between GenBank and FASTA files in automated workflows.

While BioPerl and BioPython can implement similar parsing functions, they are bulky third-party libraries and have slow parsing performance for large GenBank files. The original Perl version only relies on core modules without extra dependencies, featuring lightweight and fast parsing characteristics.

This Rust port inherits all features and parameter behaviors of the original Perl version, and brings the following improvements:
1. No Perl runtime required, natively executable on Windows, Linux and macOS
2. Higher parsing performance, multi-threaded batch processing of multiple input files
3. Strictly consistent format parsing & parameter logic with the original tool
4. Native support for mainstream compressed formats without extra Perl compression modules

## Supported Input Formats
- Genbank flat file (`.gb`, `.gbk`, `.gbff`), starts with `LOCUS`
- EMBL flat file (`.embl`), starts with `ID`
- GFF with trailing FASTA sequence (`.gff`, `.gff3`), starts with `##gff`
- FASTA (`.fasta`, `.fa`, `.fna`, `.ffn`), starts with `>`
- FASTQ (`.fastq`, `.fq`), starts with `@`
- CLUSTAL multiple sequence alignment (`.clw`, `.clu`), starts with `CLUSTAL` / `MUSCLE`
- Stockholm alignment format (`.sth`), starts with `# STOCKHOLM`
- GFA genome assembly graph (`.gfa`), first field is uppercase letter followed by tab
- PDB protein structure file (`.pdb`), starts with `HEADER`

### Supported Compressed Input Formats
Automatically decompress these formats without manual pre-decompression:
- gzip (`.gz`)
- bzip2 (`.bz2`)
- zip (`.zip`)

## Installation
### Build from source (Cross-platform)
```bash
# Clone your local project
git clone https://github.com/fanvanf/any2fasta_rust.git
cd any2fasta
cargo build --release
# Binary located at ./target/release/any2fasta
cp ./target/release/any2fasta /usr/local/bin
```

### Precompiled Binary
You can compile standalone binaries for Windows(x86_64), Linux(x86_64/aarch64), macOS(x86_64/aarch64) via Rust cross-compilation, no additional runtime dependencies required.

### Test Installation
```bash
# Check version
any2fasta -V

# Show full help document
any2fasta -h
```

## Usage
```
NAME
  any2fasta 0.9.0 (Rust Rewrite)
SYNOPSIS
  Convert various biological sequence formats into standard FASTA
  Rust port inspired by https://github.com/tseemann/any2fasta
USAGE
  any2fasta [OPTIONS] FILE.{gb,fa,fq,gff,gfa,clw,sth}[.gz,bz2,zip] [...] > output.fasta
```

### Options
```
  -q, --quiet          No runtime log output, only print error messages
  -k, --skip           Skip corrupted input files instead of terminating the program
  -n, --purify         Replace non-standard nucleotide residues (not A/C/G/T) with N;
                       non-standard amino acid residues will be replaced with X;
                       replacement character case follows the nearest valid base/amino acid
  -l, --lowercase      Convert all output sequences to lowercase
  -u, --uppercase      Convert all output sequences to uppercase
  -g, --include-version Append sequence version number to ID from GenBank / EMBL files
  -s, --strip          Remove description text after sequence ID for FASTA / FASTQ
  -h, --help           Print help information
  -V, --version        Print version information and exit
```

## Usage Examples
```bash
# Convert GenBank to FASTA
any2fasta ref.gbk > ref.fna

# Convert FASTA file (acts like cat)
any2fasta input.fasta > output.fasta

# Parse GFF file with embedded FASTA sequence
any2fasta prokka.gff > prokka.fna

# Read input from standard input
any2fasta - < genome.gbk > genome.fasta

# Auto decompress gzip compressed file
any2fasta genes.gff.gz > genes.ffn

# Batch process multiple files + stdin
any2fasta 1.gb 2.fa.gz 3.gff.bz2 - > combined.fa

# Pipe output to compressed file
any2fasta R1.fq.gz | bzip2 > R1.fa.bz2

# Convert CLUSTAL alignment, keep gap character '-'
any2fasta -q 23S.clw > 23S.aln

# Convert Stockholm alignment, convert '.' gap to standard '-'
any2fasta pfam4321.sth > pfam4321.aln
```

## Parameter Details
1. `-n --purify`: Replace invalid non-ATCG nucleotide characters with `N`; invalid amino acid characters with `X`. Gap character `-` is preserved. The case of replacement characters inherits the nearest valid preceding sequence character.
2. `-l --lowercase`: Force the entire output sequence to lowercase.
3. `-u --uppercase`: Force the entire output sequence to uppercase.
4. `-q --quiet`: Suppress all normal runtime prompt logs, only keep error outputs.
5. `-k --skip`: Print warning messages for corrupted input files and continue processing subsequent files instead of exiting abnormally.
6. `-g --include-version`: Append the sequence version number to the primary accession ID parsed from GenBank/EMBL files.
7. `-s --strip`: Truncate all descriptive content after the first whitespace/tab of the FASTA/FASTQ sequence ID, retaining only the primary identifier.

## Compatibility
This Rust version strictly replicates the parsing rules, parameter logic and output format of the original Perl `any2fasta` tool. All test cases of the original project can be reused for functional verification. The native binary can run independently on Windows, Linux and macOS without installing Perl or any third-party runtime environments.
