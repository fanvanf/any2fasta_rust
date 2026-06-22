use anyhow::{Context, Result};
use clap::Parser;
use std::cell::Cell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

// -----------------------------------------------------------------------------
// Command line arguments with detailed help
// -----------------------------------------------------------------------------
#[derive(Parser)]
#[command(
    version,
    about = "Convert various sequence formats into FASTA",
    long_about = None
)]
struct Args {
    /// No output while running, only errors
    #[arg(short = 'q')]
    quiet: bool,

    /// Skip, don't die, on bad input files
    #[arg(short = 'k')]
    skip: bool,

    /// Replace non-[AGTC] with 'N' (or 'X' if -p)
    #[arg(short = 'n')]
    purify: bool,

    /// Force protein mode (disables auto-detect)
    #[arg(short = 'p')]
    protein: bool,

    /// Lowercase the sequence
    #[arg(short = 'l')]
    lowercase: bool,

    /// Uppercase the sequence
    #[arg(short = 'u')]
    uppercase: bool,

    /// Include VERSION (GENBANK,EMBL)
    #[arg(short = 'g')]
    include_version: bool,

    /// Strip sequence descriptions (FASTA,FASTQ)
    #[arg(short = 's')]
    strip: bool,

    /// Input files (use '-' for stdin)
    #[arg(value_name = "FILE", required = true)]
    files: Vec<String>,
}

// -----------------------------------------------------------------------------
// Global configuration
// -----------------------------------------------------------------------------
struct Config {
    quiet: bool,
    skip: bool,
    purify: bool,
    protein: Cell<bool>,
    lowercase: bool,
    uppercase: bool,
    include_version: bool,
    strip: bool,
}

impl Config {
    fn from_args(args: &Args) -> Self {
        Config {
            quiet: args.quiet,
            skip: args.skip,
            purify: args.purify,
            protein: Cell::new(args.protein),
            lowercase: args.lowercase,
            uppercase: args.uppercase,
            include_version: args.include_version,
            strip: args.strip,
        }
    }
}

// -----------------------------------------------------------------------------
// Messaging helpers
// -----------------------------------------------------------------------------
macro_rules! msg {
    ($config:expr, $($arg:tt)*) => {
        if !$config.quiet {
            eprintln!($($arg)*);
        }
    };
}
macro_rules! warn {
    ($config:expr, $($arg:tt)*) => {
        if !$config.quiet {
            eprintln!("WARNING: {}", format!($($arg)*));
        }
    };
}
macro_rules! error_exit {
    ($($arg:tt)*) => {
        eprintln!("ERROR: {}", format!($($arg)*));
        std::process::exit(1);
    };
}

// -----------------------------------------------------------------------------
// Open input with automatic decompression
// -----------------------------------------------------------------------------
fn open_input(path: &str) -> Result<Box<dyn Read>> {
    if path == "-" {
        return Ok(Box::new(io::stdin()));
    }
    let file = File::open(path).with_context(|| format!("Cannot open '{}'", path))?;
    let path_lower = path.to_lowercase();

    if path_lower.ends_with(".gz") {
        Ok(Box::new(flate2::read::MultiGzDecoder::new(file)))
    } else if path_lower.ends_with(".bz2") {
        Ok(Box::new(bzip2::read::BzDecoder::new(file)))
    } else if path_lower.ends_with(".zst") || path_lower.ends_with(".zstd") {
        Ok(Box::new(zstd::stream::read::Decoder::new(file)?))
    } else if path_lower.ends_with(".zip") {
        let mut archive = zip::ZipArchive::new(file)?;
        let mut entry = archive.by_index(0).context("Zip archive has no entries")?;
        let mut content = Vec::new();
        entry.read_to_end(&mut content)?;
        Ok(Box::new(io::Cursor::new(content)))
    } else {
        Ok(Box::new(file))
    }
}

// -----------------------------------------------------------------------------
// Sequence utilities
// -----------------------------------------------------------------------------
fn is_protein(seq: &str) -> bool {
    let len = seq.len();
    if len == 0 {
        return false;
    }
    let non_atgc = seq.chars().filter(|c| !"ATCGatcg".contains(*c)).count() as f64;
    non_atgc / (len as f64) > 0.666
}
fn purify_seq(config: &Config, seq: &str) -> String {
    let mut s = String::new();
    if config.purify {
        // 记录上一个合法字符的大小写，默认大写
        let mut last_upper = true;
        if config.protein.get() {
            let allowed = "ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy-";
            for ch in seq.chars() {
                if ch.is_whitespace() {
                    s.push(ch);
                    continue;
                }
                if allowed.contains(ch) {
                    s.push(ch);
                    // 更新最近合法字符大小写
                    if ch.is_ascii_uppercase() {
                        last_upper = true;
                    } else {
                        last_upper = false;
                    }
                } else {
                    // 非法字符跟随上一个合法字符大小写
                    if last_upper {
                        s.push('X');
                    } else {
                        s.push('x');
                    }
                }
            }
        } else {
            // DNA 模式
            let allowed = "ACGTacgt-";
            for ch in seq.chars() {
                if ch.is_whitespace() {
                    s.push(ch);
                    continue;
                }
                if allowed.contains(ch) {
                    s.push(ch);
                    if ch.is_ascii_uppercase() {
                        last_upper = true;
                    } else {
                        last_upper = false;
                    }
                } else {
                    if last_upper {
                        s.push('N');
                    } else {
                        s.push('n');
                    }
                }
            }
        }
    } else {
        s = seq.to_string();
    }

    // 互斥大小写转换
    if config.lowercase {
        s = s.to_lowercase();
    } else if config.uppercase {
        s = s.to_uppercase();
    }
    s
}

fn purify_id(config: &Config, id: &str) -> String {
    let mut s = id.to_string();
    if config.strip {
        if let Some(p) = s.find(|c: char| c == ' ' || c == '\t') {
            s.truncate(p);
        }
    }
    s
}

// -----------------------------------------------------------------------------
// Parser functions
// -----------------------------------------------------------------------------
fn parse_fasta(config: &Config, lines: &mut Vec<String>) -> Result<usize> {
    let mut count = 0;
    let mut seq_line_num = 0;
    let mut first_seq = true;
    for line in lines.iter_mut() {
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('>') {
            *line = purify_id(config, line);
            count += 1;
            seq_line_num = 0;
            println!("{}", line);
        } else {
            seq_line_num += 1;
            if count == 1 && seq_line_num == 1 && first_seq {
                let is_prot = is_protein(line);
                config.protein.set(is_prot);
                msg!(config, "Auto-detected alphabet: {}", if is_prot { "PROT" } else { "DNA" });
                first_seq = false;
            }
            let pur = purify_seq(config, line);
            println!("{}", pur);
        }
    }
    Ok(count)
}

fn parse_fastq(config: &Config, lines: &mut Vec<String>) -> Result<usize> {
    let mut count = 0;
    let mut i = 0;
    while i + 3 < lines.len() {
        let seq_line = lines[i + 1].clone();
        let id_line = &mut lines[i];
        *id_line = purify_id(config, id_line);
        println!(">{}", &id_line[1..]);
        let pur = purify_seq(config, &seq_line);
        println!("{}", pur);
        count += 1;
        i += 4;
    }
    Ok(count)
}

fn parse_gff(config: &Config, lines: &mut Vec<String>) -> Result<usize> {
    let mut fasta_lines = Vec::new();
    let mut in_fasta = false;
    for line in lines.iter() {
        if line.starts_with("##FASTA") {
            in_fasta = true;
            continue;
        }
        if in_fasta && (line.starts_with('>') || !line.trim().is_empty()) {
            fasta_lines.push(line.clone());
        }
    }
    if fasta_lines.is_empty() {
        return Ok(0);
    }
    parse_fasta(config, &mut fasta_lines)
}

fn parse_gfa(config: &Config, lines: &mut Vec<String>) -> Result<usize> {
    let mut count = 0;
    for line in lines {
        if line.starts_with('S') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                println!(">{}", parts[1]);
                let pur = purify_seq(config, parts[2]);
                println!("{}", pur);
                count += 1;
            }
        }
    }
    Ok(count)
}

fn parse_genbank(config: &Config, lines: &mut Vec<String>) -> Result<usize> {
    let mut acc = String::new();
    let mut dna = String::new();
    let mut in_seq = false;
    let mut count = 0;

    for line in lines {
        let line = line.trim_end();
        if line.starts_with("//") {
            if !acc.is_empty() && !dna.is_empty() {
                println!(">{}", acc);
                let pur = purify_seq(config, &dna);
                println!("{}", pur);
                count += 1;
            }
            in_seq = false;
            dna.clear();
            acc.clear();
            continue;
        }
        if line.starts_with("ORIGIN") {
            in_seq = true;
            continue;
        }
        if in_seq {
            if line.len() > 10 {
                let s = &line[10..];
                let s = s.replace(' ', "").replace('\t', "");
                dna.push_str(&s);
            }
        } else {
            if line.starts_with("LOCUS") {
                let fields: Vec<&str> = line.split_whitespace().collect();
                if fields.len() >= 2 {
                    acc = fields[1].to_string();
                    if line.contains(" aa ") || line.contains(" AA ") {
                        config.protein.set(true);
                    }
                }
            }
            if config.include_version && line.starts_with("VERSION") {
                let fields: Vec<&str> = line.split_whitespace().collect();
                if fields.len() >= 2 {
                    acc = fields[1].to_string();
                }
            }
        }
    }
    if !acc.is_empty() && !dna.is_empty() {
        println!(">{}", acc);
        let pur = purify_seq(config, &dna);
        print!("{}", pur);
        count += 1;
    }
    Ok(count)
}

fn parse_embl(config: &Config, lines: &mut Vec<String>) -> Result<usize> {
    let mut acc = String::new();
    let mut dna = String::new();
    let mut in_seq = false;
    let mut count = 0;

    for line in lines {
        let line = line.trim_end();
        if line.starts_with("//") {
            if !acc.is_empty() && !dna.is_empty() {
                println!(">{}", acc);
                let pur = purify_seq(config, &dna);
                println!("{}", pur);
                count += 1;
            }
            in_seq = false;
            dna.clear();
            acc.clear();
            continue;
        }
        if line.starts_with("SQ") {
            in_seq = true;
            if line.contains(" AA ") || line.contains(" aa ") {
                config.protein.set(true);
            }
            continue;
        }
        if in_seq {
            let mut s = line.to_string();
            s.retain(|c| !c.is_whitespace() && !c.is_ascii_digit());
            if !s.is_empty() {
                dna.push_str(&s);
            }
        } else {
            if line.starts_with("ID") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    let id_part = parts[1];
                    let id = id_part.trim_end_matches(';');
                    acc = id.to_string();
                    if config.include_version {
                        let fields_semi: Vec<&str> = line.split(';').collect();
                        for seg in fields_semi {
                            let seg_trim = seg.trim();
                            if seg_trim.starts_with("SV ") {
                                let sv_num = seg_trim.strip_prefix("SV ").unwrap().trim();
                                acc.push('.');
                                acc.push_str(sv_num);
                                break;
                            }
                        }
                    }
                }
            }
            if line.starts_with("ID") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    let id_part = parts[1];
                    let id = id_part.trim_end_matches(';');
                    acc = id.to_string();
                    if config.include_version {
                        // 提取 SV 数字
                        let fields_semi: Vec<&str> = line.split(';').collect();
                        for seg in fields_semi {
                            let seg_trim = seg.trim();
                            if seg_trim.starts_with("SV ") {
                                let sv_num = seg_trim.strip_prefix("SV ").unwrap().trim();
                                acc.push('.');
                                acc.push_str(sv_num);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    if !acc.is_empty() && !dna.is_empty() {
        println!(">{}", acc);
        let pur = purify_seq(config, &dna);
        print!("{}", pur);
        count += 1;
    }
    Ok(count)
}

fn parse_clustal(config: &Config, lines: &mut Vec<String>) -> Result<usize> {
    let mut seqs: HashMap<String, String> = HashMap::new();
    let mut order = Vec::new();

    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with("CLUSTAL") || trimmed.starts_with("MUSCLE") {
            continue;
        }
        // 分割所有空白符，取前两个非空字段
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 2 {
            continue;
        }
        let id = fields[0];
        let seq = fields[1];
        // 去除行尾可能的数字（但这里没有，因为已经分割）
        let seq_clean = seq.trim_end_matches(|c: char| c.is_ascii_digit());
        if seq_clean.is_empty() {
            continue;
        }
        // 检查是否只含字母和连字符（粗略检查）
        if !seq_clean.chars().all(|c| c.is_ascii_alphabetic() || c == '-') {
            continue;
        }
        if !seqs.contains_key(id) {
            order.push(id.to_string());
        }
        let entry = seqs.entry(id.to_string()).or_default();
        entry.push_str(seq_clean);
    }

    for id in &order {
        if let Some(seq) = seqs.get(id) {
            println!(">{}", id);
            let pur = purify_seq(config, seq);
            println!("{}", pur);
        }
    }
    Ok(order.len())
}

fn parse_stockholm(config: &Config, lines: &mut Vec<String>) -> Result<usize> {
    let mut seqs: HashMap<String, String> = HashMap::new();
    let mut order = Vec::new();

    for line in lines {
        let trimmed = line.trim();
        if trimmed.starts_with('#') || trimmed.starts_with("//") || trimmed.is_empty() {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 2 {
            continue;
        }
        let id = fields[0];
        let seq = fields[1];
        let seq_clean = seq.trim_end_matches(|c: char| c.is_ascii_digit());
        if seq_clean.is_empty() {
            continue;
        }
        // 将 '.' 替换为 '-' (Stockholm 缺省用 . 表示 gap)
        let seq_replaced = seq_clean.replace('.', "-");
        if !seq_replaced.chars().all(|c| c.is_ascii_alphabetic() || c == '-') {
            continue;
        }
        if !seqs.contains_key(id) {
            order.push(id.to_string());
        }
        let entry = seqs.entry(id.to_string()).or_default();
        entry.push_str(&seq_replaced);
        entry.push('\n');
    }

    for id in &order {
        if let Some(seq) = seqs.get(id) {
            println!(">{}", id);
            let pur = purify_seq(config, seq);
            print!("{}", pur);
        }
    }
    Ok(order.len())
}

fn parse_pdb(config: &Config, lines: &mut Vec<String>) -> Result<usize> {
    let aa_map: HashMap<&str, char> = [
        ("ALA", 'A'), ("ARG", 'R'), ("ASN", 'N'), ("ASP", 'D'),
        ("CYS", 'C'), ("GLU", 'E'), ("GLN", 'Q'), ("GLY", 'G'),
        ("HIS", 'H'), ("ILE", 'I'), ("LEU", 'L'), ("LYS", 'K'),
        ("MET", 'M'), ("PHE", 'F'), ("PRO", 'P'), ("SER", 'S'),
        ("THR", 'T'), ("TRP", 'W'), ("TYR", 'Y'), ("VAL", 'V'),
        ("SEC", 'U'), ("PYL", 'O'), ("ASX", 'B'), ("GLX", 'Z'),
        ("XAA", 'X'), ("TER", '*'),
    ]
    .iter()
    .cloned()
    .collect();

    let mut seqs: HashMap<String, String> = HashMap::new();
    let mut pdb_id = "unknown".to_string();

    for line in &*lines {
        if line.starts_with("HEADER") {
            if let Some(id) = line.split_whitespace().last() {
                pdb_id = id.to_string();
            }
            break;
        }
    }

    for line in &*lines {
        if line.starts_with("SEQRES") {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 4 {
                continue;
            }
            let chain = fields[2];
            let key = format!("{}-{}", pdb_id, chain);
            let seq = seqs.entry(key).or_default();
            for &aa3 in &fields[4..] {
                seq.push(*aa_map.get(aa3).unwrap_or(&'X'));
            }
        }
    }

    for (id, seq) in &seqs {
        println!(">{}", id);
        let pur = purify_seq(config, seq);
        println!("{}", pur);
    }
    Ok(seqs.len())
}

// -----------------------------------------------------------------------------
// Format detection
// -----------------------------------------------------------------------------
fn detect_genbank(line: &str) -> bool {
    line.starts_with("LOCUS ")
}
fn detect_uniprot(line: &str) -> bool {
    line.starts_with("ID") && line.contains("Reviewed")
}
fn detect_embl(line: &str) -> bool {
    line.starts_with("ID ") && !line.contains("Reviewed")
}
fn detect_fasta(line: &str) -> bool {
    line.starts_with('>')
}
fn detect_fastq(line: &str) -> bool {
    line.starts_with('@')
}
fn detect_gff(line: &str) -> bool {
    line.starts_with("##gff")
}
fn detect_clustal(line: &str) -> bool {
    line.starts_with("CLUST") || line.starts_with("MUSCL")
}
fn detect_stockholm(line: &str) -> bool {
    line.starts_with("# STOCKHOLM ")
}
fn detect_gfa(line: &str) -> bool {
    // 支持制表符或空格 (但 Perl 只认制表符，这里兼容)
    if line.len() < 2 {
        return false;
    }
    let first = line.chars().next().unwrap();
    first.is_ascii_uppercase() && (line.chars().nth(1).unwrap() == '\t' || line.chars().nth(1).unwrap() == ' ')
}
fn detect_pdb(line: &str) -> bool {
    line.starts_with("HEADER ")
}

type ParserFn = fn(&Config, &mut Vec<String>) -> Result<usize>;

struct Format {
    name: &'static str,
    detect: fn(&str) -> bool,
    parser: ParserFn,
}

const FORMATS: &[Format] = &[
    Format {
        name: "GENBANK",
        detect: detect_genbank,
        parser: parse_genbank,
    },
    Format {
        name: "UNIPROT",
        detect: detect_uniprot,
        parser: parse_embl,
    },
    Format {
        name: "EMBL",
        detect: detect_embl,
        parser: parse_embl,
    },
    Format {
        name: "FASTA",
        detect: detect_fasta,
        parser: parse_fasta,
    },
    Format {
        name: "FASTQ",
        detect: detect_fastq,
        parser: parse_fastq,
    },
    Format {
        name: "GFF",
        detect: detect_gff,
        parser: parse_gff,
    },
    Format {
        name: "CLUSTAL",
        detect: detect_clustal,
        parser: parse_clustal,
    },
    Format {
        name: "STOCKHOLM",
        detect: detect_stockholm,
        parser: parse_stockholm,
    },
    Format {
        name: "GFA",
        detect: detect_gfa,
        parser: parse_gfa,
    },
    Format {
        name: "PDB",
        detect: detect_pdb,
        parser: parse_pdb,
    },
];

// -----------------------------------------------------------------------------
// Main processing loop
// -----------------------------------------------------------------------------
fn run(args: Args) -> Result<()> {
    let config = Config::from_args(&args);
    msg!(config, "This is any2fasta {}", env!("CARGO_PKG_VERSION"));

    let mut processed = 0;
    for infile in args.files {
        let path = if infile == "-" { "-" } else { &infile };
        msg!(config, "Opening '{}'", path);

        if path != "-" && Path::new(path).is_dir() {
            let msg = format!("'{}' is a directory not a file", path);
            if config.skip {
                warn!(config, "{}", msg);
                continue;
            } else {
                error_exit!("{}", msg);
            }
        }

        let reader = match open_input(path) {
            Ok(r) => r,
            Err(e) => {
                let msg = format!("Cannot open '{}': {}", path, e);
                if config.skip {
                    warn!(config, "{}", msg);
                    continue;
                } else {
                    error_exit!("{}", msg);
                }
            }
        };

        let mut buf_reader = BufReader::new(reader);
        let mut first_line = String::new();
        if buf_reader.read_line(&mut first_line)? == 0 {
            let msg = format!("The input '{}' appears to be empty", path);
            if config.skip {
                warn!(config, "{}", msg);
                continue;
            } else {
                error_exit!("{}", msg);
            }
        }

        // 去除第一行的换行符，与 lines() 保持一致
        if first_line.ends_with('\n') {
            first_line.pop();
            if first_line.ends_with('\r') {
                first_line.pop();
            }
        }

        let mut lines = Vec::new();
        lines.push(first_line);
        for line in buf_reader.lines() {
            lines.push(line?);
        }

        let mut detected = false;
        for fmt in FORMATS {
            if (fmt.detect)(&lines[0]) {
                msg!(config, "Detected {} format", fmt.name);
                let nseq = (fmt.parser)(&config, &mut lines)?;
                msg!(config, "Wrote {} sequences from {} file.", nseq, fmt.name);
                if nseq == 0 {
                    let msg = format!("No sequences found in '{}'", path);
                    if config.skip {
                        warn!(config, "{}", msg);
                    } else {
                        error_exit!("{}", msg);
                    }
                }
                detected = true;
                break;
            }
        }

        if !detected {
            let msg = format!("Unfamiliar format with first line: {}", lines[0].trim());
            if config.skip {
                warn!(config, "{}", msg);
            } else {
                error_exit!("{}", msg);
            }
        }
        processed += 1;
    }

    msg!(config, "Processed {} files.", processed);
    msg!(config, "Done.");
    Ok(())
}

fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("{}", e);
        std::process::exit(1);
    }
}