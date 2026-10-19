use pilercr_rs::{
    Array, Config, Sequence, detect, fasta_out,
    report::{self, CoordinateMode},
};
use std::{
    ffi::OsStr,
    fs,
    path::{Path, PathBuf},
    process::Command,
    time::Instant,
};

fn golden_arrays(text: &str, sequence: &Sequence) -> Vec<Array> {
    let lines: Vec<_> = text.lines().collect();
    let mut arrays = Vec::new();
    let mut i = 0;
    while i < lines.len() && lines[i] != "SUMMARY BY SIMILARITY" {
        if let Some(id) = lines[i].strip_prefix("Array ") {
            let label = lines[i + 1].strip_prefix('>').unwrap();
            let contig = sequence
                .contigs
                .iter()
                .position(|c| c.label == label)
                .unwrap();
            i += 5;
            let mut array = Array {
                id: id.parse().unwrap(),
                ..Array::default()
            };
            while !lines[i].starts_with("==========") {
                let line = lines[i];
                if array.repeats.is_empty() {
                    let local = line[..10].trim().parse::<usize>().unwrap() - 1;
                    array.pos = sequence.local_to_global(contig, local).unwrap();
                }
                // Four fixed-width columns (10, 6, 6, 6), each followed by two spaces.
                array.left_flanks.push(line.as_bytes()[36..46].to_vec());
                let mut fields = line[50..].split_whitespace();
                array
                    .repeats
                    .push(fields.next().unwrap().as_bytes().to_vec());
                array
                    .spacers
                    .push(fields.next().unwrap_or("").as_bytes().to_vec());
                i += 1;
            }
            array.aligned_consensus = lines[i + 1]
                .split_whitespace()
                .last()
                .unwrap()
                .as_bytes()
                .to_vec();
            array.consensus = array
                .aligned_consensus
                .iter()
                .copied()
                .filter(|&c| c != b'-')
                .collect();
            for repeat in &mut array.repeats {
                assert_eq!(repeat.len(), array.aligned_consensus.len());
                for (c, &consensus) in repeat.iter_mut().zip(&array.aligned_consensus) {
                    if *c == b'.' {
                        *c = consensus;
                    }
                }
            }
            arrays.push(array);
        }
        i += 1;
    }
    arrays
}

fn array_difference(actual: &[Array], expected: &[Array], sequence: &Sequence) -> Option<String> {
    if actual.len() != expected.len() {
        let positions = |arrays: &[Array]| {
            arrays
                .iter()
                .map(|a| {
                    let (ci, pos) = sequence.global_to_local(a.pos).unwrap();
                    format!("{}:{}", sequence.contigs[ci].label, pos + 1)
                })
                .collect::<Vec<_>>()
                .join(",")
        };
        return Some(format!(
            "arrays {}/{}; positions actual=[{}] expected=[{}]; earliest possible: kmer/filter/DP/pile/graph or array acceptance",
            actual.len(),
            expected.len(),
            positions(actual),
            positions(expected)
        ));
    }
    for (a, e) in actual.iter().zip(expected) {
        if a.id != e.id || a.pos != e.pos {
            return Some(format!(
                "array {} id/position actual={}/{} expected={}/{}; earliest possible: pile/graph selection or array localization",
                e.id, a.id, a.pos, e.id, e.pos
            ));
        }
        if a.repeats.len() != e.repeats.len() {
            return Some(format!(
                "array {} copies {}/{}; earliest possible: pile/graph selection or array bounds",
                e.id,
                a.repeats.len(),
                e.repeats.len()
            ));
        }
        if a.consensus != e.consensus || a.aligned_consensus != e.aligned_consensus {
            return Some(format!(
                "array {} consensus actual={:?} expected={:?}; earliest possible: array MSA/bounds/consensus",
                e.id,
                String::from_utf8_lossy(&a.aligned_consensus),
                String::from_utf8_lossy(&e.aligned_consensus)
            ));
        }
        // Detail reports expose only the last ten flank characters, left-padded
        // with spaces; internal genome-derived flanks may be longer.
        let left_flanks: Vec<_> = a
            .left_flanks
            .iter()
            .map(|flank| {
                format!(
                    "{:>10}",
                    String::from_utf8_lossy(&flank[flank.len().saturating_sub(10)..])
                )
                .into_bytes()
            })
            .collect();
        for (name, actual, expected) in [
            ("repeats", &a.repeats, &e.repeats),
            ("spacers", &a.spacers, &e.spacers),
            ("left flanks", &left_flanks, &e.left_flanks),
        ] {
            if actual != expected {
                let row = actual
                    .iter()
                    .zip(expected)
                    .position(|(a, e)| a != e)
                    .unwrap_or(actual.len().min(expected.len()));
                return Some(format!(
                    "array {} {name} row {}; earliest possible: array localization/MSA/bounds",
                    e.id,
                    row + 1
                ));
            }
        }
    }
    None
}

fn byte_difference(name: &str, actual: &str, expected: &str) -> Option<String> {
    if actual == expected {
        return None;
    }
    let offset = actual
        .bytes()
        .zip(expected.bytes())
        .position(|(a, e)| a != e)
        .unwrap_or(actual.len().min(expected.len()));
    let line = expected.as_bytes()[..offset.min(expected.len())]
        .iter()
        .filter(|&&c| c == b'\n')
        .count();
    let excerpt = |text: &str| {
        text.lines()
            .nth(line)
            .unwrap_or("<EOF>")
            .chars()
            .take(160)
            .collect::<String>()
    };
    Some(format!(
        "{name} byte {offset}, line {} actual={:?} expected={:?}",
        line + 1,
        excerpt(actual),
        excerpt(expected)
    ))
}

struct Scratch(PathBuf);

impl Scratch {
    /// Derive the candidate name for `prefix` and sequence `seq` under `root`.
    /// Pure naming, so a test can pre-create the exact path an allocation will
    /// try first.
    fn candidate(root: &Path, prefix: &str, seq: u64) -> PathBuf {
        root.join(format!("{prefix}-{}-{seq}", std::process::id()))
    }

    /// A unique, empty scratch directory removed on drop. The name combines a
    /// caller prefix, the process id, and a per-process sequence so parallel
    /// tests cannot collide. Creation is atomic (`create_dir`), so only the
    /// process that successfully creates the directory owns it: a stale path
    /// left by a crashed run (or a reused process id) makes `create_dir` fail
    /// with `AlreadyExists`, and the retry moves to a fresh name rather than
    /// deleting a directory this process does not own.
    fn unique(prefix: &str) -> Self {
        Self::unique_in(&std::env::temp_dir(), prefix)
    }

    /// Like `unique`, but allocates under an explicit `root` instead of the
    /// shared temp directory. Tests use this to keep every path they touch
    /// inside a scratch tree they own.
    fn unique_in(root: &Path, prefix: &str) -> Self {
        use std::sync::atomic::AtomicU64;
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        Self::unique_from(root, prefix, &COUNTER)
    }

    /// Like `unique_in`, but draws sequence numbers from the caller's counter.
    /// Tests pass a private counter so they can pre-create the exact candidate
    /// the first attempt will try without racing the shared `unique` counter
    /// used by the rest of the suite.
    fn unique_from(root: &Path, prefix: &str, counter: &std::sync::atomic::AtomicU64) -> Self {
        use std::sync::atomic::Ordering;
        loop {
            let seq = counter.fetch_add(1, Ordering::Relaxed);
            let path = Self::candidate(root, prefix, seq);
            match fs::create_dir(&path) {
                Ok(()) => return Scratch(path),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => {
                    panic!(
                        "scratch directory {} is not creatable: {error}",
                        path.display()
                    )
                }
            }
        }
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// Write a minimal placeholder executable so default-path resolution sees the
/// same file (and, on unix, executable-bit) shape as a real checkout.
fn write_executable(path: &Path) {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).expect("executable parent directory");
    }
    fs::write(path, b"#!/bin/sh\nexit 0\n").expect("placeholder executable");
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        fs::set_permissions(path, fs::Permissions::from_mode(0o755))
            .expect("placeholder executable permissions");
    }
}

/// Whether `dir` is this crate's package root: it holds a `Cargo.toml` that
/// declares the expected package name.
fn is_package_dir(dir: &Path) -> bool {
    let Ok(manifest) = fs::read_to_string(dir.join("Cargo.toml")) else {
        return false;
    };
    manifest
        .lines()
        .any(|line| line.trim() == "name = \"pilercr-rs\"")
}

/// Resolve this package's directory (`crates/pilercr-rs`) from runtime
/// context.
///
/// Cargo bakes the absolute manifest path into the test binary at compile time,
/// so a relocated build tree (a copied `target/` whose artifacts are newer than
/// the sources) can leave a stale path that still exists but points fixture
/// lookups at a different checkout. Rather than trust that compiled-in path,
/// prefer a runtime `CARGO_MANIFEST_DIR` when it names this package, and
/// otherwise fall back to Cargo's package-root working directory, validating
/// the package manifest in both cases.
fn resolve_package_dir(runtime_manifest: Option<&OsStr>, current_dir: &Path) -> PathBuf {
    if let Some(dir) = runtime_manifest {
        let dir = Path::new(dir);
        if is_package_dir(dir) {
            return dir.to_path_buf();
        }
    }
    assert!(
        is_package_dir(current_dir),
        "test working directory {} is not the pilercr-rs package root; run tests via cargo",
        current_dir.display()
    );
    current_dir.to_path_buf()
}

fn package_dir() -> PathBuf {
    resolve_package_dir(
        std::env::var_os("CARGO_MANIFEST_DIR").as_deref(),
        &std::env::current_dir().expect("test working directory is unavailable"),
    )
}

fn oracle_dir(root: &Path, explicit: Option<&OsStr>) -> Result<PathBuf, String> {
    let Some(value) = explicit else {
        return Ok(root.join("test-data/local-oracle"));
    };
    let path = Path::new(value)
        .canonicalize()
        .map_err(|e| format!("invalid PILERCR_ORACLE_DIR {value:?}: {e}"))?;
    if !path.is_dir() {
        return Err(format!(
            "PILERCR_ORACLE_DIR is not a directory: {}",
            path.display()
        ));
    }
    Ok(path)
}

/// The optional C++ reference for consensus FASTA comparison. Only an explicit
/// `PILERCR_REFERENCE_BINARY` selects one; nothing is picked up implicitly.
fn reference_binary(explicit: Option<&OsStr>) -> Result<Option<PathBuf>, String> {
    let Some(value) = explicit else {
        return Ok(None);
    };
    let path = Path::new(value)
        .canonicalize()
        .map_err(|e| format!("invalid PILERCR_REFERENCE_BINARY {value:?}: {e}"))?;
    if !path.is_file() {
        return Err(format!(
            "reference binary is not a file: {}",
            path.display()
        ));
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let metadata = fs::metadata(&path).map_err(|e| e.to_string())?;
        if metadata.permissions().mode() & 0o111 == 0 {
            return Err(format!(
                "reference binary is not executable: {}",
                path.display()
            ));
        }
    }
    Ok(Some(path))
}

/// One parity case: an input FASTA, its reference legacy report, and the
/// reference consensus FASTA when that is committed alongside the report.
struct Case {
    name: String,
    input: PathBuf,
    golden: String,
    consensus: Option<String>,
}

/// Compare one case against its goldens. `reference` is an optional
/// (binary, scratch directory) pair, used only when no committed consensus exists.
fn fixture(case: &Case, reference: Option<(&Path, &Path)>) -> Result<Vec<String>, String> {
    let input = case.input.as_path();
    let golden = case.golden.as_str();
    let (input_name, count) = golden
        .lines()
        .find_map(|line| {
            let (name, suffix) = line.rsplit_once(": ")?;
            Some((name, suffix.strip_suffix(" putative CRISPR arrays found.")?))
        })
        .ok_or("missing golden input/count line")?;
    let sequence = Sequence::from_path(input).map_err(|e| e.to_string())?;
    let expected = golden_arrays(golden, &sequence);
    if count.parse::<usize>().map_err(|e| e.to_string())? != expected.len() {
        return Err("golden header/detail array count differs".into());
    }
    let arrays = detect(&sequence, &Config::default()).map_err(|error| {
        let stage = if error == "filter hit coordinates are outside the sequence" {
            "earliest failing stage: filter::merge_hits validates kmer::filter coordinates"
        } else {
            "native detect failed before output comparison"
        };
        format!("{error}; {stage}; array/text/FASTA parity blocked")
    })?;
    let mut failures = Vec::new();
    let array_difference = array_difference(&arrays, &expected, &sequence);
    if let Some(difference) = &array_difference {
        failures.push(difference.clone());
    }
    match report::text(
        &sequence,
        &arrays,
        input_name,
        CoordinateMode::Legacy,
        false,
    ) {
        Ok(actual) => {
            if let Some(difference) = byte_difference("legacy text", &actual, golden) {
                failures.push(format!(
                    "{difference}; {}",
                    if array_difference.is_some() {
                        "upstream detection differs"
                    } else {
                        "earliest likely: report formatting/clustering"
                    }
                ));
            }
        }
        Err(error) => failures.push(format!("legacy report error: {error}")),
    }
    let expected_consensus = match (&case.consensus, reference) {
        (Some(committed), _) => Some(committed.clone()),
        (None, Some((binary, scratch))) => {
            let output = scratch.join(format!("{}.txt", case.name));
            let consensus = scratch.join(format!("{}.fasta", case.name));
            let result = Command::new(binary)
                .arg("-in")
                .arg(input)
                .arg("-out")
                .arg(&output)
                .arg("-seq")
                .arg(&consensus)
                .arg("-quiet")
                .output()
                .map_err(|e| e.to_string())?;
            if result.status.success() {
                Some(fs::read_to_string(consensus).map_err(|e| e.to_string())?)
            } else {
                failures.push(format!(
                    "reference failed {}: {}",
                    result.status,
                    String::from_utf8_lossy(&result.stderr)
                        .chars()
                        .take(160)
                        .collect::<String>()
                ));
                None
            }
        }
        (None, None) => None,
    };
    if let Some(expected) = expected_consensus {
        match fasta_out::consensus_fasta(&sequence, &arrays, false) {
            Ok(actual) => {
                if let Some(difference) = byte_difference("consensus FASTA", &actual, &expected) {
                    failures.push(format!(
                        "{difference}; {}",
                        if array_difference.is_some() {
                            "upstream detection differs"
                        } else {
                            "earliest likely: consensus FASTA output"
                        }
                    ));
                }
            }
            Err(error) => failures.push(format!("consensus FASTA error: {error}")),
        }
    }
    Ok(failures)
}

fn read_fixture(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// Excerpt ids (`NN` in `excerpt_NN.legacy.txt`), sorted.
fn excerpt_ids(dir: &Path) -> Vec<String> {
    let mut ids: Vec<String> = fs::read_dir(dir)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()))
        .filter_map(|entry| {
            let name = entry.ok()?.file_name().into_string().ok()?;
            name.strip_prefix("excerpt_")?
                .strip_suffix(".legacy.txt")
                .map(str::to_owned)
        })
        .collect();
    ids.sort();
    ids
}

fn excerpt_case(dir: &Path, id: &str) -> Case {
    Case {
        name: id.to_owned(),
        input: dir.join(format!("excerpt_{id}.fa")),
        golden: read_fixture(&dir.join(format!("excerpt_{id}.legacy.txt"))),
        consensus: Some(read_fixture(
            &dir.join(format!("excerpt_{id}.consensus.fa")),
        )),
    }
}

fn oracle_case(root: &Path, oracle: &Path, id: usize) -> Case {
    Case {
        name: id.to_string(),
        input: root.join(format!("test-data/assemblies/{id}.fasta")),
        golden: read_fixture(&oracle.join(format!("{id}.txt"))),
        consensus: None,
    }
}

/// Run every case, recording per-case problems. A panic in one case is a
/// failure for that case, not an abort of the whole run.
fn run_cases(cases: &[Case], reference: Option<(&Path, &Path)>) -> (usize, Vec<String>) {
    let mut passed = 0;
    let mut failures = Vec::new();
    for case in cases {
        let result = std::panic::catch_unwind(|| fixture(case, reference));
        let problems = match result {
            Ok(Ok(problems)) => problems,
            Ok(Err(error)) => vec![format!("pipeline/fixture error: {error}")],
            Err(error) => vec![format!(
                "panic: {}",
                error
                    .downcast_ref::<String>()
                    .map(String::as_str)
                    .or_else(|| error.downcast_ref::<&str>().copied())
                    .unwrap_or("non-string panic")
            )],
        };
        if problems.is_empty() {
            passed += 1;
        } else {
            failures.push(format!("{}: {}", case.name, problems.join(" | ")));
        }
    }
    (passed, failures)
}

#[test]
fn report_flanks_compare_displayed_suffix_without_truncating_right_flank() {
    let sequence =
        Sequence::from_reader(&b">scf1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"[..]).unwrap();
    let mut actual = Array {
        id: 1,
        pos: 10,
        consensus: b"ACGT".to_vec(),
        aligned_consensus: b"ACGT".to_vec(),
        repeats: vec![b"ACGT".to_vec(); 2],
        spacers: vec![b"TT".to_vec(), b"TATATGCGACT".to_vec()],
        left_flanks: vec![b"GENOME_PREFIXTTTGGCAGAG".to_vec(), b"TT".to_vec()],
        ..Array::default()
    };
    let text = report::text(
        &sequence,
        std::slice::from_ref(&actual),
        "fixture.fasta",
        CoordinateMode::Legacy,
        true,
    )
    .unwrap();
    let expected = golden_arrays(&text, &sequence);
    assert_eq!(
        expected[0].left_flanks,
        [b"TTTGGCAGAG", b"        TT"],
        "{text}"
    );
    assert_eq!(expected[0].spacers[1], b"TATATGCGACT");
    assert_eq!(
        array_difference(std::slice::from_ref(&actual), &expected, &sequence),
        None
    );

    actual.left_flanks[0].push(b'C');
    assert!(
        array_difference(std::slice::from_ref(&actual), &expected, &sequence)
            .unwrap()
            .contains("left flanks row 1")
    );
    actual.left_flanks[0].pop();
    actual.spacers[1].pop();
    assert!(
        array_difference(&[actual], &expected, &sequence)
            .unwrap()
            .contains("spacers row 2")
    );
}

#[test]
fn scratch_directories_are_unique_and_only_the_owner_removes_them() {
    let first = Scratch::unique("pilercr-scratch-ownership");
    let second = Scratch::unique("pilercr-scratch-ownership");
    assert_ne!(
        first.0, second.0,
        "parallel scratch directories must not collide"
    );
    assert!(first.0.is_dir() && second.0.is_dir());
    let first_path = first.0.clone();
    drop(first);
    assert!(
        !first_path.exists(),
        "the owning process removes its own scratch directory on drop"
    );
    assert!(
        second.0.is_dir(),
        "a sibling scratch directory is not removed by another owner"
    );
}

#[test]
fn scratch_unique_skips_an_occupied_candidate_and_preserves_its_contents() {
    use std::sync::atomic::AtomicU64;

    // A unique, owned parent keeps the whole collision inside a scratch tree
    // this test owns, so nothing is pre-created (or later removed) directly
    // under the shared temp directory. The parent's Drop cleans everything.
    let parent = Scratch::unique("pilercr-scratch-collision");
    let prefix = "candidate";

    // A private counter makes the first candidate deterministic: it starts at 0
    // and no other test or process shares it, so occupying `(prefix, pid, 0)`
    // forces exactly one collision on the first allocation attempt.
    let counter = AtomicU64::new(0);
    let occupied = Scratch::candidate(&parent.0, prefix, 0);
    fs::create_dir(&occupied).expect("pre-create the first candidate");
    let sentinel = occupied.join("sentinel");
    fs::write(&sentinel, b"owned by the pre-existing directory").expect("write sentinel");

    let scratch = Scratch::unique_from(&parent.0, prefix, &counter);

    assert_ne!(
        scratch.0, occupied,
        "a collision must move the allocation to a fresh candidate"
    );
    assert!(scratch.0.is_dir(), "the chosen candidate is a directory");
    assert_eq!(
        fs::read_dir(&scratch.0).expect("read chosen dir").count(),
        0,
        "the chosen candidate is empty"
    );
    assert!(
        occupied.is_dir(),
        "the pre-existing candidate survives allocation"
    );
    assert_eq!(
        fs::read(&sentinel).expect("sentinel readable after allocation"),
        b"owned by the pre-existing directory",
        "the sentinel survives allocation"
    );

    let owned = scratch.0.clone();
    drop(scratch);
    assert!(
        !owned.exists(),
        "the owning process removes its own scratch directory on drop"
    );
    assert!(
        occupied.is_dir(),
        "the pre-existing candidate survives Drop"
    );
    assert_eq!(
        fs::read(&sentinel).expect("sentinel readable after Drop"),
        b"owned by the pre-existing directory",
        "the sentinel survives Drop"
    );

    // The parent Scratch owns the occupied candidate and the sentinel; dropping
    // it here removes the whole tree, including anything the collision left.
    drop(parent);
    assert!(
        !occupied.exists(),
        "dropping the parent Scratch cleans up the pre-existing candidate"
    );
}

#[test]
fn oracle_resolution_defaults_use_local_oracle_and_no_implicit_reference() {
    // A workspace-shaped root with a locally built reference present: it must
    // not be selected unless PILERCR_REFERENCE_BINARY names it.
    let scratch = Scratch::unique("pilercr-oracle-resolution");
    let root = scratch.0.as_path();

    let oracle = root.join("test-data/local-oracle");
    fs::create_dir_all(&oracle).unwrap();
    assert_eq!(oracle_dir(root, None).unwrap(), oracle);

    let built = root.join("target/reference-trace/pilercr");
    write_executable(&built);
    assert_eq!(
        reference_binary(None).unwrap(),
        None,
        "a built reference must not be selected without PILERCR_REFERENCE_BINARY"
    );
    assert_eq!(
        reference_binary(Some(built.as_os_str())).unwrap(),
        Some(built.canonicalize().unwrap()),
        "an explicit reference is selected"
    );
}

#[test]
fn oracle_resolution_explicit_paths_are_canonical() {
    let root = package_dir();
    assert_eq!(
        oracle_dir(&root, Some(OsStr::new("."))).unwrap(),
        Path::new(".").canonicalize().unwrap()
    );
    let binary = std::env::current_exe().unwrap();
    assert_eq!(
        reference_binary(Some(binary.as_os_str())).unwrap(),
        Some(binary.canonicalize().unwrap())
    );
}

#[test]
fn oracle_resolution_invalid_explicit_paths_fail() {
    let root = package_dir();
    let missing = root.join("__missing_reference_parity_oracle__");
    for value in [OsStr::new(""), missing.as_os_str()] {
        assert!(oracle_dir(&root, Some(value)).is_err());
        assert!(reference_binary(Some(value)).is_err());
    }
    let file = root.join("Cargo.toml");
    assert!(oracle_dir(&root, Some(file.as_os_str())).is_err());
    assert!(reference_binary(Some(root.as_os_str())).is_err());
    #[cfg(unix)]
    assert!(reference_binary(Some(file.as_os_str())).is_err());
}

#[test]
fn committed_excerpts_match_reference_goldens() {
    let dir = package_dir().join("tests/fixtures/excerpts");
    let ids = excerpt_ids(&dir);
    assert!(!ids.is_empty(), "no excerpt goldens in {}", dir.display());
    let cases: Vec<Case> = ids.iter().map(|id| excerpt_case(&dir, id)).collect();
    let (passed, failures) = run_cases(&cases, None);
    let summary = format!(
        "committed excerpt parity: {passed}/{} passed, {} failed",
        cases.len(),
        failures.len()
    );
    eprintln!("{summary}");
    assert!(failures.is_empty(), "{summary}\n{}", failures.join("\n"));
}

#[test]
fn package_dir_prefers_runtime_manifest_and_rejects_stale_existing_path() {
    let package = package_dir();
    assert!(is_package_dir(&package), "{}", package.display());
    // The workspace root exists and holds a `Cargo.toml`, but not this package's
    // manifest: a stale compiled-in path like this must not be trusted.
    let stale = package.join("../..");
    assert!(
        stale.is_dir() && !is_package_dir(&stale),
        "{}",
        stale.display()
    );
    assert_eq!(
        resolve_package_dir(Some(stale.as_os_str()), &package),
        package
    );
    // A runtime manifest naming this package wins over the working directory.
    assert_eq!(
        resolve_package_dir(Some(package.as_os_str()), &stale),
        package
    );
}

#[test]
#[ignore = "requires test-data/; run with --ignored"]
fn all_52_native_detection_matches_local_oracle() {
    let root = package_dir().join("../..");
    let oracle = oracle_dir(&root, std::env::var_os("PILERCR_ORACLE_DIR").as_deref())
        .unwrap_or_else(|error| panic!("{error}"));
    let explicit_reference = std::env::var_os("PILERCR_REFERENCE_BINARY");
    let reference =
        reference_binary(explicit_reference.as_deref()).unwrap_or_else(|error| panic!("{error}"));
    let scratch = Scratch::unique("pilercr-reference-parity");
    if reference.is_none() {
        eprintln!("PILERCR_REFERENCE_BINARY not set; optional consensus FASTA comparison skipped");
    }
    if cfg!(debug_assertions) {
        eprintln!(
            "All52 detection is expensive; prefer cargo test --release -p pilercr-rs -- --ignored"
        );
    }
    let start = Instant::now();
    let cases: Vec<Case> = (1..=52).map(|id| oracle_case(&root, &oracle, id)).collect();
    let reference_pair = reference
        .as_deref()
        .map(|binary| (binary, scratch.0.as_path()));
    let (passed, failures) = run_cases(&cases, reference_pair);
    let summary = format!(
        "all52 native parity [oracle={}]: {passed}/52 passed, {} failed, consensus FASTA {}, elapsed {:.2}s",
        oracle.display(),
        failures.len(),
        if reference.is_some() {
            "enabled"
        } else {
            "skipped"
        },
        start.elapsed().as_secs_f64()
    );
    eprintln!("{summary}");
    assert!(failures.is_empty(), "{summary}\n{}", failures.join("\n"));
}

/// One line of the "SUMMARY BY SIMILARITY" table. These are the fields that
/// make an original-2011 vs modern divergence visible without shipping either
/// large report or assembly as a fixture.
#[derive(Clone, Debug, PartialEq, Eq)]
struct SummarySignature {
    position: usize,
    length: usize,
    copies: usize,
    repeat_len: usize,
    spacer_len: usize,
    consensus: String,
}

fn signature(
    position: usize,
    length: usize,
    copies: usize,
    repeat_len: usize,
    spacer_len: usize,
    consensus: &str,
) -> SummarySignature {
    SummarySignature {
        position,
        length,
        copies,
        repeat_len,
        spacer_len,
        consensus: consensus.to_owned(),
    }
}

/// Signatures in the shipped original 2011 reports (`test-data/ref-output/`).
/// Assembly 9's first row has two trailing alignment gaps in its consensus;
/// that is part of the 2011 report, not a transcription error.
fn original_2011_signatures(assembly: usize) -> Vec<SummarySignature> {
    match assembly {
        9 => vec![
            signature(518070, 1066, 18, 28, 33, "GTTTTCTCCACGTATGTGGAGGTGATCC--"),
            signature(519291, 213, 4, 30, 31, "ATTTTCTCCACGTATGTGGAGGTGATCCTA"),
        ],
        16 => vec![
            signature(
                643021,
                1814,
                28,
                35,
                30,
                "GTTTTTGTACTCTCAAGATTTAAGTAACTATACAA",
            ),
            signature(
                892312,
                1058,
                15,
                36,
                37,
                "GATATAGACCTAATTACCTCGAGAGGGGACGGAAAC",
            ),
            signature(
                1375385,
                564,
                9,
                36,
                30,
                "GTTTTGGAACCATTCGAAACAACACAGCTCTAAAAC",
            ),
        ],
        17 => vec![signature(
            1152331,
            828,
            13,
            36,
            30,
            "GATCTAAACCTTATTGATCTAACAACCATCTAAAAC",
        )],
        other => panic!("no original-2011 signature recorded for assembly {other}"),
    }
}

/// The modern baseline this project targets: local-oracle output and a modern
/// rebuild, which the Rust implementation matches. Assembly 16 arrays 2 and 3
/// are unchanged; only array 1 moves and gains copies.
fn modern_signatures(assembly: usize) -> Vec<SummarySignature> {
    match assembly {
        9 => vec![signature(
            518071,
            1431,
            24,
            27,
            34,
            "TTTTCTCCACGTATGTGGAGGTGATCC",
        )],
        16 => vec![
            signature(
                641966,
                2870,
                44,
                36,
                29,
                "GTTTTTGTACTCTCAAGATTTAAGTAACTATACAAC",
            ),
            signature(
                892312,
                1058,
                15,
                36,
                37,
                "GATATAGACCTAATTACCTCGAGAGGGGACGGAAAC",
            ),
            signature(
                1375385,
                564,
                9,
                36,
                30,
                "GTTTTGGAACCATTCGAAACAACACAGCTCTAAAAC",
            ),
        ],
        17 => vec![signature(
            1152595,
            564,
            9,
            36,
            30,
            "GATCTAAACCTTATTGATCTAACAACCATCTAAAAC",
        )],
        other => panic!("no modern signature recorded for assembly {other}"),
    }
}

fn parse_summary_signatures(report: &str) -> Vec<SummarySignature> {
    let mut in_similarity = false;
    let mut signatures = Vec::new();
    for line in report.lines() {
        if line == "SUMMARY BY SIMILARITY" {
            in_similarity = true;
            continue;
        }
        if line == "SUMMARY BY POSITION" {
            break;
        }
        if !in_similarity {
            continue;
        }
        let fields: Vec<_> = line.split_whitespace().collect();
        // Header, separator, footnote and blank lines do not start with an
        // array number. Cluster footnotes are not part of the signature.
        if fields.len() < 9 || fields[0].parse::<usize>().is_err() {
            continue;
        }
        signatures.push(SummarySignature {
            position: fields[2].parse().unwrap(),
            length: fields[3].parse().unwrap(),
            copies: fields[4].parse().unwrap(),
            repeat_len: fields[5].parse().unwrap(),
            spacer_len: fields[6].parse().unwrap(),
            consensus: fields[8].to_owned(),
        });
    }
    signatures
}

fn report_array_count(report: &str) -> usize {
    report
        .lines()
        .find_map(|line| {
            let (_, suffix) = line.rsplit_once(": ")?;
            suffix
                .strip_suffix(" putative CRISPR arrays found.")?
                .parse()
                .ok()
        })
        .expect("report has no array-count header")
}

#[test]
#[ignore = "requires test-data/; run with --ignored"]
fn original_2011_divergences_on_assemblies_9_16_17_are_pinned() {
    let root = package_dir().join("../..");
    let oracle = oracle_dir(&root, std::env::var_os("PILERCR_ORACLE_DIR").as_deref())
        .unwrap_or_else(|error| panic!("{error}"));

    for assembly in [9, 16, 17] {
        let original_expected = original_2011_signatures(assembly);
        let modern_expected = modern_signatures(assembly);
        let original_report =
            read_fixture(&root.join(format!("test-data/ref-output/{assembly}.txt")));
        let oracle_report = read_fixture(&oracle.join(format!("{assembly}.txt")));

        assert_eq!(
            report_array_count(&original_report),
            original_expected.len(),
            "assembly {assembly}: original 2011 array count"
        );
        assert_eq!(
            parse_summary_signatures(&original_report),
            original_expected,
            "assembly {assembly}: original 2011 summary signatures"
        );
        assert_eq!(
            report_array_count(&oracle_report),
            modern_expected.len(),
            "assembly {assembly}: modern rebuild array count"
        );
        assert_eq!(
            parse_summary_signatures(&oracle_report),
            modern_expected,
            "assembly {assembly}: modern rebuild summary signatures"
        );

        let input = root.join(format!("test-data/assemblies/{assembly}.fasta"));
        let sequence = Sequence::from_path(&input).unwrap_or_else(|error| {
            panic!("assembly {assembly} input {}: {error}", input.display())
        });
        let arrays = detect(&sequence, &Config::default())
            .unwrap_or_else(|error| panic!("assembly {assembly}: {error}"));
        let rust_report = report::text(
            &sequence,
            &arrays,
            "assembly.fasta",
            CoordinateMode::Legacy,
            true,
        )
        .unwrap_or_else(|error| panic!("assembly {assembly}: {error}"));
        assert_eq!(
            report_array_count(&rust_report),
            modern_expected.len(),
            "assembly {assembly}: Rust array count"
        );
        assert_eq!(
            parse_summary_signatures(&rust_report),
            modern_expected,
            "assembly {assembly}: Rust summary signatures must match the modern baseline"
        );
        assert_ne!(
            original_expected, modern_expected,
            "assembly {assembly}: the recorded 2011/modern divergence disappeared"
        );
    }
}
