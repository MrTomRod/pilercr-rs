use pilercr_rs::{
    Config, Sequence, detect, fasta_out,
    report::{self, CoordinateMode},
    trace,
};
use std::{
    collections::HashSet,
    env, fs,
    io::{self, Write},
    path::{Path, PathBuf},
};

const HELP: &str = "Usage: pilercr-rs -in INPUT.fa [-in MORE.fa ...] -out REPORT.txt [options]
  -seq PATH            Consensus FASTA output
  -trimseqs            Trim similar consensus sequences (requires -seq)
  -quiet               Suppress progress messages
  -noinfo              Omit explanatory report text
  --legacy-coords      Use reference gap-shifted text coordinates
  --json PATH          JSON output (corrected coordinates)
  --gff PATH           GFF3 output (corrected coordinates)
  --tsv PATH           TSV output (corrected coordinates)
  --trace-dir PATH     Checkpoint TSVs in a new or empty directory (-trace-dir)
  -minrepeat N         Minimum repeat length (16)
  -maxrepeat N         Maximum repeat length (64)
  -minspacer N         Minimum spacer length (8)
  -maxspacer N         Maximum spacer length (64)
  -minarray N          Minimum repeats per array, at least 3 (3)
  -minrepeatratio F    Minimum repeat-length ratio (0.9)
  -minspacerratio F    Minimum spacer-length ratio (0.75)
  -mincons F           Minimum conservation (0.9)
  -minhitlength N      Minimum DP hit length, greater than 4 (16)
  -minid F             Minimum fractional DP identity (0.94)
  -diameter N          Filter search diameter (512)
  -maxmem BYTES        Parameter-selection memory budget (536870912)
  --help, -help        Show this help
  --version, -version  Show version
Multiple inputs are combined in order as separate contigs. Threads are not supported.
";

#[derive(Debug, Default)]
struct Args {
    inputs: Vec<PathBuf>,
    out: Option<PathBuf>,
    seq: Option<PathBuf>,
    json: Option<PathBuf>,
    gff: Option<PathBuf>,
    tsv: Option<PathBuf>,
    trace_dir: Option<PathBuf>,
    trim: bool,
    quiet: bool,
    no_info: bool,
    legacy: bool,
    config: Config,
}

#[derive(Debug)]
enum Command {
    Run(Box<Args>),
    Help,
    Version,
}

fn number<T: std::str::FromStr>(value: &str, flag: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("invalid value for {flag}: {value}"))
}

fn parse(args: impl IntoIterator<Item = String>) -> Result<Command, String> {
    let mut args = args.into_iter();
    let mut result = Args::default();
    let mut seen = HashSet::new();
    while let Some(flag) = args.next() {
        if flag == "--help" || flag == "-help" || flag == "-h" || flag == "-options" {
            return Ok(Command::Help);
        }
        if flag == "--version" || flag == "-version" {
            return Ok(Command::Version);
        }
        if flag == "--threads" || flag == "-threads" {
            return Err("--threads is not implemented; detection is single-threaded".into());
        }
        if flag != "-in" && !seen.insert(flag.clone()) {
            return Err(format!("duplicate option: {flag}"));
        }
        match flag.as_str() {
            "-quiet" => result.quiet = true,
            "-noinfo" => result.no_info = true,
            "-trimseqs" => result.trim = true,
            "--legacy-coords" => result.legacy = true,
            "-in" | "-out" | "-seq" | "--json" | "--gff" | "--tsv" | "-minrepeat"
            | "-maxrepeat" | "-minspacer" | "-maxspacer" | "-minarray" | "-minrepeatratio"
            | "-minspacerratio" | "-mincons" | "-minhitlength" | "-minid" | "-diameter"
            | "-maxmem" | "-trace-dir" | "--trace-dir" => {
                let value = args
                    .next()
                    .filter(|v| !v.is_empty() && !v.starts_with('-'))
                    .ok_or_else(|| format!("missing value for {flag}"))?;
                match flag.as_str() {
                    "-in" => result.inputs.push(value.into()),
                    "-out" => result.out = Some(value.into()),
                    "-seq" => result.seq = Some(value.into()),
                    "--json" => result.json = Some(value.into()),
                    "--gff" => result.gff = Some(value.into()),
                    "--tsv" => result.tsv = Some(value.into()),
                    "-trace-dir" | "--trace-dir" => {
                        if result.trace_dir.replace(value.into()).is_some() {
                            return Err("duplicate option: --trace-dir".into());
                        }
                    }
                    "-minrepeat" => result.config.array.min_repeat = number(&value, &flag)?,
                    "-maxrepeat" => result.config.array.max_repeat = number(&value, &flag)?,
                    "-minspacer" => result.config.array.min_spacer = number(&value, &flag)?,
                    "-maxspacer" => result.config.array.max_spacer = number(&value, &flag)?,
                    "-minarray" => result.config.array.min_array = number(&value, &flag)?,
                    "-minrepeatratio" => {
                        result.config.array.min_repeat_ratio = number(&value, &flag)?
                    }
                    "-minspacerratio" => {
                        result.config.array.min_spacer_ratio = number(&value, &flag)?
                    }
                    "-mincons" => result.config.array.min_cons = number(&value, &flag)?,
                    "-minhitlength" => result.config.dp.min_hit_length = number(&value, &flag)?,
                    "-minid" => result.config.dp.min_identity = number(&value, &flag)?,
                    "-diameter" => result.config.diameter = number(&value, &flag)?,
                    "-maxmem" => result.config.memory_budget = number(&value, &flag)?,
                    _ => unreachable!(),
                }
            }
            _ => return Err(format!("unknown option: {flag}")),
        }
    }
    if result.inputs.is_empty() || result.out.is_none() {
        return Err("-in and -out are required (use --help for options)".into());
    }
    if result.trim && result.seq.is_none() {
        return Err("-trimseqs requires -seq".into());
    }
    result.config.validate()?;
    Ok(Command::Run(Box::new(result)))
}

fn resolved(path: &Path) -> Result<PathBuf, String> {
    if path.exists() {
        return fs::canonicalize(path).map_err(|e| format!("{}: {e}", path.display()));
    }
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let parent = fs::canonicalize(parent).map_err(|e| format!("{}: {e}", path.display()))?;
    let name = path
        .file_name()
        .ok_or_else(|| format!("invalid output path: {}", path.display()))?;
    Ok(parent.join(name))
}

fn run(args: Args) -> Result<(), String> {
    let outputs: Vec<_> = [&args.out, &args.seq, &args.json, &args.gff, &args.tsv]
        .into_iter()
        .flatten()
        .collect();
    let mut paths = HashSet::new();
    for input in &args.inputs {
        paths.insert(resolved(input)?);
    }
    for output in &outputs {
        if !paths.insert(resolved(output)?) {
            return Err(format!(
                "output path aliases an input or another output: {}",
                output.display()
            ));
        }
    }
    if let Some(directory) = &args.trace_dir {
        let directory = resolved(directory)?;
        if paths
            .iter()
            .any(|path| path.starts_with(&directory) || directory.starts_with(path))
        {
            return Err("trace directory aliases an input or output path".into());
        }
    }
    let mut sequence = Sequence {
        bases: Vec::new(),
        contigs: Vec::new(),
    };
    for input in &args.inputs {
        if !args.quiet {
            eprintln!("Reading {}", input.display());
        }
        let mut next =
            Sequence::from_path(input).map_err(|e| format!("{}: {e}", input.display()))?;
        if !sequence.contigs.is_empty() && !next.contigs.is_empty() {
            let mut padding = 32 - sequence.bases.len() % 32;
            if padding < 50 {
                padding += 32;
            }
            sequence.bases.resize(sequence.bases.len() + padding, b'!');
        }
        let offset = sequence.bases.len();
        for contig in &mut next.contigs {
            contig.start = contig
                .start
                .checked_add(offset)
                .ok_or("contig offset overflow")?;
        }
        sequence.bases.extend(next.bases);
        sequence.contigs.extend(next.contigs);
    }
    let arrays = if let Some(directory) = &args.trace_dir {
        trace::with_directory(directory, || detect(&sequence, &args.config))?
    } else {
        detect(&sequence, &args.config)?
    };
    let input_name = args
        .inputs
        .iter()
        .map(|p| p.display().to_string())
        .collect::<Vec<_>>()
        .join(", ");
    let mode = if args.legacy {
        CoordinateMode::Legacy
    } else {
        CoordinateMode::Corrected
    };
    let mut rendered = vec![(
        args.out.as_ref().unwrap(),
        report::text(&sequence, &arrays, &input_name, mode, args.no_info)?,
    )];
    if let Some(path) = &args.seq {
        rendered.push((
            path,
            fasta_out::consensus_fasta(&sequence, &arrays, args.trim)?,
        ));
    }
    if let Some(path) = &args.json {
        rendered.push((path, report::json(&sequence, &arrays)?));
    }
    if let Some(path) = &args.gff {
        rendered.push((path, report::gff(&sequence, &arrays)?));
    }
    if let Some(path) = &args.tsv {
        rendered.push((path, report::tsv(&sequence, &arrays)?));
    }
    for (path, content) in rendered {
        fs::write(path, content).map_err(|e| format!("{}: {e}", path.display()))?;
    }
    if !args.quiet {
        eprintln!("{} putative CRISPR arrays found.", arrays.len());
    }
    Ok(())
}

fn main() {
    let status = match parse(env::args().skip(1)) {
        Ok(Command::Help) => io::stdout()
            .write_all(HELP.as_bytes())
            .map_err(|e| e.to_string()),
        Ok(Command::Version) => writeln!(io::stdout(), "pilercr-rs {}", env!("CARGO_PKG_VERSION"))
            .map_err(|e| e.to_string()),
        Ok(Command::Run(args)) => run(*args),
        Err(error) => {
            eprintln!("pilercr-rs: {error}");
            std::process::exit(2);
        }
    };
    if let Err(error) = status {
        eprintln!("pilercr-rs: {error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn command(args: &[&str]) -> Result<Command, String> {
        parse(args.iter().map(|s| s.to_string()))
    }

    #[test]
    fn parses_reference_and_structured_options() {
        let Command::Run(args) = command(&[
            "-in",
            "a.fa",
            "-in",
            "b.fa",
            "-out",
            "r.txt",
            "-seq",
            "r.fa",
            "-trimseqs",
            "-quiet",
            "-noinfo",
            "--legacy-coords",
            "--json",
            "r.json",
            "--gff",
            "r.gff",
            "--tsv",
            "r.tsv",
            "-trace-dir",
            "trace",
            "-minarray",
            "4",
            "-minid",
            "0.95",
        ])
        .unwrap() else {
            panic!()
        };
        assert_eq!(args.inputs.len(), 2);
        assert_eq!(args.trace_dir, Some(PathBuf::from("trace")));
        assert!(args.trim && args.quiet && args.no_info && args.legacy);
        assert_eq!(args.config.array.min_array, 4);
        assert_eq!(args.config.dp.min_identity, 0.95);
    }

    #[test]
    fn writes_empty_reports_and_rejects_output_aliases() {
        let directory = env::temp_dir().join(format!(
            "pilercr-rs-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&directory).unwrap();
        let input = directory.join("empty.fa");
        fs::write(&input, b"").unwrap();
        let args = Args {
            inputs: vec![input.clone(), input.clone()],
            out: Some(directory.join("report.txt")),
            seq: Some(directory.join("consensus.fa")),
            json: Some(directory.join("arrays.json")),
            gff: Some(directory.join("arrays.gff")),
            tsv: Some(directory.join("arrays.tsv")),
            quiet: true,
            no_info: true,
            ..Args::default()
        };
        run(args).unwrap();
        assert!(
            fs::read_to_string(directory.join("report.txt"))
                .unwrap()
                .contains("0 putative CRISPR arrays found.")
        );
        assert_eq!(fs::read(directory.join("consensus.fa")).unwrap(), b"");
        assert!(
            fs::read_to_string(directory.join("arrays.json"))
                .unwrap()
                .contains("[]")
        );
        assert!(
            fs::read_to_string(directory.join("arrays.gff"))
                .unwrap()
                .starts_with("##gff-version 3")
        );
        assert!(!fs::read(directory.join("arrays.tsv")).unwrap().is_empty());
        assert!(
            run(Args {
                inputs: vec![input.clone()],
                out: Some(input.clone()),
                quiet: true,
                ..Args::default()
            })
            .is_err()
        );
        assert_eq!(fs::read(&input).unwrap(), b"");
        let report = directory.join("report.txt");
        let original = fs::read(&report).unwrap();
        let trace_dir = directory.join("trace");
        let traced_args = |trace_dir| Args {
            inputs: vec![input.clone(), input.clone()],
            out: Some(report.clone()),
            trace_dir: Some(trace_dir),
            quiet: true,
            no_info: true,
            ..Args::default()
        };
        run(traced_args(trace_dir.clone())).unwrap();
        assert_eq!(fs::read(&report).unwrap(), original);
        assert_eq!(
            fs::read_dir(&trace_dir).unwrap().count(),
            trace::STAGES.len()
        );
        for path in [trace_dir, input.clone(), report.clone(), directory.clone()] {
            assert!(run(traced_args(path)).is_err());
            assert_eq!(fs::read(&report).unwrap(), original);
            assert_eq!(fs::read(&input).unwrap(), b"");
        }
        let collision = directory.join("collision");
        assert!(
            run(Args {
                inputs: vec![input.clone()],
                out: Some(collision.join("01-parameters.tsv")),
                trace_dir: Some(collision.clone()),
                quiet: true,
                ..Args::default()
            })
            .is_err()
        );
        assert!(!collision.exists());
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn rejects_unknown_missing_duplicate_and_invalid_options() {
        for flags in [
            vec![],
            vec!["--threads", "4"],
            vec!["-in"],
            vec!["-in", "a", "-out", "b", "-trimseqs"],
            vec!["-in", "a", "-out", "b", "-mincons", "NaN"],
            vec!["-in", "a", "-out", "b", "-minarray", "2"],
            vec!["-in", "a", "-out", "b", "-quiet", "-quiet"],
            vec!["-in", "a", "-out", "b", "--bogus"],
            vec!["-in", "a", "-out", "b", "--trace-dir"],
            vec![
                "-in",
                "a",
                "-out",
                "b",
                "--trace-dir",
                "x",
                "-trace-dir",
                "y",
            ],
        ] {
            assert!(command(&flags).is_err(), "{flags:?}");
        }
        assert!(matches!(command(&["--help"]), Ok(Command::Help)));
        assert!(matches!(command(&["--version"]), Ok(Command::Version)));
    }
}
