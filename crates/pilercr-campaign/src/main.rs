use pilercr_campaign::{
    CRASH_EXEMPT, classify_signal, generate_case, mutate_repeats, mutation_step_substitutions,
    validated_array_count,
};
use serde_json::json;
#[cfg(unix)]
use std::os::unix::process::CommandExt;
#[cfg(unix)]
use std::os::unix::process::ExitStatusExt;
use std::{
    env, fs,
    io::Read,
    path::{Path, PathBuf},
    process::{Child, Command, Stdio},
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};

#[derive(Debug)]
struct Opt {
    target: usize,
    max_attempts: usize,
    max_misses: usize,
    max_mutations: usize,
    mutation_steps: usize,
    seed: u64,
    negative_every: usize,
    output: PathBuf,
    tmpdir: PathBuf,
    reference: Option<PathBuf>,
    rust: Option<PathBuf>,
    timeout: Duration,
    build_timeout: Duration,
}
fn seconds(value: &str, flag: &str) -> Result<Duration, String> {
    let n = value.parse::<f64>().map_err(|_| format!("bad {flag}"))?;
    if !n.is_finite() || n <= 0.0 {
        return Err(format!("{flag} must be positive finite seconds"));
    }
    Duration::try_from_secs_f64(n).map_err(|_| format!("{flag} is out of range"))
}
/// Repository-local scratch directory (gitignored). Default home for temporary
/// build/work files and the default campaign output.
fn scratch_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../scratch")
}
fn parse() -> Result<Opt, String> {
    let mut o = Opt {
        target: 1000,
        max_attempts: 0,
        max_misses: 100,
        max_mutations: 0,
        mutation_steps: 32,
        seed: 20261017,
        negative_every: 10,
        output: scratch_dir().join("synthetic-campaign-out"),
        tmpdir: scratch_dir(),
        reference: None,
        rust: None,
        timeout: Duration::from_secs(120),
        build_timeout: Duration::from_secs(600),
    };
    let mut a = env::args().skip(1);
    while let Some(k) = a.next() {
        if k == "--help" || k == "-h" {
            println!(
                "pilercr-campaign [--target N] [--max-attempts N] [--max-misses N] [--max-mutations N] [--mutation-steps N] [--seed N] [--negative-every N] [--output PATH] [--tmpdir PATH] [--reference-binary PATH] [--rust-binary PATH] [--timeout SEC] [--build-timeout SEC]"
            );
            std::process::exit(0)
        }
        let v = a.next().ok_or_else(|| format!("missing value for {k}"))?;
        match k.as_str() {
            "--target" => o.target = v.parse().map_err(|_| "bad target")?,
            "--max-attempts" => o.max_attempts = v.parse().map_err(|_| "bad max-attempts")?,
            "--max-misses" => o.max_misses = v.parse().map_err(|_| "bad max-misses")?,
            "--max-mutations" => o.max_mutations = v.parse().map_err(|_| "bad max-mutations")?,
            "--mutation-steps" => o.mutation_steps = v.parse().map_err(|_| "bad mutation-steps")?,
            "--seed" => o.seed = v.parse().map_err(|_| "bad seed")?,
            "--negative-every" => o.negative_every = v.parse().map_err(|_| "bad negative-every")?,
            "--output" => o.output = v.into(),
            "--tmpdir" => o.tmpdir = v.into(),
            "--reference-binary" => o.reference = Some(v.into()),
            "--rust-binary" => o.rust = Some(v.into()),
            "--timeout" => o.timeout = seconds(&v, "--timeout")?,
            "--build-timeout" => o.build_timeout = seconds(&v, "--build-timeout")?,
            _ => return Err(format!("unknown option {k}")),
        }
    }
    if o.target == 0 {
        return Err("--target must be >= 1".into());
    }
    Ok(o)
}
fn canonical_destination(path: &Path) -> Result<PathBuf, String> {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        env::current_dir().map_err(|e| e.to_string())?.join(path)
    };
    let mut tail = Vec::new();
    let mut ancestor = absolute.as_path();
    while fs::symlink_metadata(ancestor).is_err() {
        let name = ancestor
            .file_name()
            .ok_or_else(|| format!("cannot resolve output path: {}", path.display()))?;
        tail.push(name.to_os_string());
        ancestor = ancestor
            .parent()
            .ok_or_else(|| format!("cannot resolve output path: {}", path.display()))?;
    }
    let mut resolved = fs::canonicalize(ancestor)
        .map_err(|e| format!("cannot resolve output ancestor {}: {e}", ancestor.display()))?;
    for part in tail.into_iter().rev() {
        if part == ".." {
            resolved.pop();
        } else if part != "." {
            resolved.push(part);
        }
    }
    Ok(resolved)
}
fn overlap(a: &Path, b: &Path) -> bool {
    a.starts_with(b) || b.starts_with(a)
}
fn now() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    format!("{secs}")
}
fn normalized(report: &[u8], input: &Path) -> Vec<u8> {
    let input_name = input.file_name().and_then(|x| x.to_str());
    if validated_array_count(report, &input.to_string_lossy()).is_none() || input_name.is_none() {
        return report.to_vec();
    }
    let mut out = Vec::new();
    for line in report.split_inclusive(|b| *b == b'\n') {
        let body = line.strip_suffix(b"\n").unwrap_or(line);
        let body = body.strip_suffix(b"\r").unwrap_or(body);
        let parsed = std::str::from_utf8(body)
            .ok()
            .and_then(|s| s.rsplit_once(": "))
            .and_then(|(path, tail)| {
                tail.strip_suffix(" putative CRISPR arrays found.")
                    .map(|digits| (path, digits))
            });
        let current = parsed.is_some_and(|(path, digits)| {
            !digits.is_empty()
                && digits.bytes().all(|b| b.is_ascii_digit())
                && Path::new(path).file_name().and_then(|x| x.to_str()) == input_name
        });
        if current {
            let delimiter = line.windows(2).rposition(|window| window == b": ");
            if let Some(delimiter) = delimiter {
                out.extend_from_slice(b"<INPUT>: ");
                out.extend_from_slice(&line[delimiter + 2..]);
            } else {
                out.extend_from_slice(line)
            }
        } else {
            out.extend_from_slice(line)
        }
    }
    out
}
fn sha_bytes(bytes: &[u8]) -> Result<String, String> {
    let mut c = Command::new("sha256sum")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|e| e.to_string())?;
    use std::io::Write;
    c.stdin
        .take()
        .unwrap()
        .write_all(bytes)
        .map_err(|e| e.to_string())?;
    let o = c.wait_with_output().map_err(|e| e.to_string())?;
    if !o.status.success() {
        return Err("sha256sum failed".into());
    }
    Ok(String::from_utf8_lossy(&o.stdout)
        .split_whitespace()
        .next()
        .unwrap_or_default()
        .to_owned())
}
fn hash(p: &Path) -> Result<String, String> {
    let out = Command::new("sha256sum")
        .arg(p)
        .output()
        .map_err(|e| format!("sha256sum {}: {e}", p.display()))?;
    if !out.status.success() {
        return Err(format!("sha256sum failed for {}", p.display()));
    }
    String::from_utf8(out.stdout)
        .map_err(|e| e.to_string())?
        .split_whitespace()
        .next()
        .map(str::to_owned)
        .ok_or_else(|| "sha256sum returned no digest".into())
}
#[derive(Debug)]
struct Run {
    code: Option<i32>,
    signal: Option<i32>,
    timed_out: bool,
    stdout: Vec<u8>,
    stderr: Vec<u8>,
    elapsed: f64,
    command: Vec<String>,
}
const MAX_LOG_BYTES: u64 = 2 * 1024 * 1024;
fn finish_child(
    mut child: Child,
    timeout: Duration,
) -> (Option<std::process::ExitStatus>, bool, f64) {
    let start = Instant::now();
    let mut timed_out = false;
    let status = loop {
        match child.try_wait() {
            Ok(Some(s)) => break Some(s),
            Ok(None) if start.elapsed() < timeout => std::thread::sleep(Duration::from_millis(10)),
            _ => {
                timed_out = true;
                let _ = child.kill();
                kill_process_group(child.id());
                break child.wait().ok();
            }
        }
    };
    kill_process_group(child.id());
    (status, timed_out, start.elapsed().as_secs_f64())
}
fn read_log(path: &Path) -> Result<Vec<u8>, String> {
    let file = fs::File::open(path).map_err(|e| format!("read log {}: {e}", path.display()))?;
    let mut bytes = Vec::new();
    file.take(MAX_LOG_BYTES)
        .read_to_end(&mut bytes)
        .map_err(|e| e.to_string())?;
    Ok(bytes)
}
fn sidecar(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_os_string();
    name.push(suffix);
    PathBuf::from(name)
}
fn cleanup_case_artifacts(input: &Path, reports: &[&Path]) -> Result<(), String> {
    let mut paths = vec![input.to_path_buf()];
    for report in reports {
        paths.push(report.to_path_buf());
        paths.push(sidecar(report, ".stdout"));
        paths.push(sidecar(report, ".stderr"));
    }
    for path in paths {
        match fs::remove_file(&path) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(format!("remove temporary artifact {}: {e}", path.display())),
        }
    }
    Ok(())
}
fn spawn_logged(
    mut command: Command,
    stdout_path: &Path,
    stderr_path: &Path,
) -> Result<(Child, Vec<Child>), String> {
    let stdout_file = fs::File::create(stdout_path)
        .map_err(|e| format!("create log {}: {e}", stdout_path.display()))?;
    let stderr_file = fs::File::create(stderr_path)
        .map_err(|e| format!("create log {}: {e}", stderr_path.display()))?;
    #[cfg(unix)]
    command.process_group(0);
    let mut child = command
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| e.to_string())?;
    let out = child.stdout.take().ok_or("missing child stdout")?;
    let err = child.stderr.take().ok_or("missing child stderr")?;
    let output = match capped_logger(out, stdout_file) {
        Ok(x) => x,
        Err(e) => {
            let _ = child.kill();
            kill_process_group(child.id());
            let _ = child.wait();
            return Err(e);
        }
    };
    match capped_logger(err, stderr_file) {
        Ok(error) => Ok((child, vec![output, error])),
        Err(error) => {
            let _ = child.kill();
            kill_process_group(child.id());
            let _ = child.wait();
            let mut output = output;
            let _ = output.kill();
            let _ = output.wait();
            Err(error)
        }
    }
}
fn capped_logger(input: impl Into<Stdio>, file: fs::File) -> Result<Child, String> {
    // Keep consuming after the prefix is saved; otherwise the tool sees EPIPE
    // as soon as head reaches its byte cap.
    let mut cmd = Command::new("sh");
    cmd.args([
        "-c",
        "head -c \"$1\"; cat >/dev/null",
        "campaign-log",
        &MAX_LOG_BYTES.to_string(),
    ])
    .stdin(input)
    .stdout(Stdio::from(file))
    .stderr(Stdio::null());
    #[cfg(unix)]
    cmd.process_group(0);
    cmd.spawn()
        .map_err(|e| format!("start capped log writer: {e}"))
}
fn finish_loggers(mut loggers: Vec<Child>) {
    let grace = Instant::now() + Duration::from_millis(50);
    loop {
        let running = loggers
            .iter_mut()
            .any(|logger| logger.try_wait().ok().flatten().is_none());
        if !running || Instant::now() >= grace {
            break;
        }
        std::thread::sleep(Duration::from_millis(2));
    }
    for logger in &mut loggers {
        if logger.try_wait().ok().flatten().is_none() {
            kill_process_group(logger.id());
        }
        let _ = logger.wait();
    }
}
fn kill_process_group(pid: u32) {
    #[cfg(unix)]
    {
        let _ = Command::new("kill")
            .args(["-KILL", "--", &format!("-{pid}")])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status();
    }
    #[cfg(not(unix))]
    {
        let _ = pid;
    }
}
#[allow(clippy::too_many_arguments)] // one call-site per evidence kind; a struct adds noise
fn persist_case(
    out: &Path,
    id: &str,
    input: &[u8],
    input_path: &Path,
    cpp: &Run,
    rust: &Run,
    cr: &[u8],
    rr: &[u8],
    cn: Option<u64>,
    rn: Option<u64>,
    extra: serde_json::Value,
) -> Result<PathBuf, String> {
    let dir = out.join(id);
    fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
    fs::write(dir.join("input.fa"), input).map_err(|e| e.to_string())?;
    if !cr.is_empty() {
        fs::write(dir.join("cpp.report.txt"), cr).map_err(|e| e.to_string())?
    }
    if !rr.is_empty() {
        fs::write(dir.join("rust.report.txt"), rr).map_err(|e| e.to_string())?
    }
    fs::write(dir.join("cpp.stdout"), &cpp.stdout).map_err(|e| e.to_string())?;
    fs::write(dir.join("cpp.stderr"), &cpp.stderr).map_err(|e| e.to_string())?;
    fs::write(dir.join("rust.stdout"), &rust.stdout).map_err(|e| e.to_string())?;
    fs::write(dir.join("rust.stderr"), &rust.stderr).map_err(|e| e.to_string())?;
    let a = normalized(cr, input_path);
    let b = normalized(rr, input_path);
    let first = first_difference(&a, &b);
    let metadata = json!({"input_sha256":sha_bytes(input)?,"first_difference":first,"cpp":{"command":cpp.command,"returncode":cpp.code,"signal":cpp.signal,"timed_out":cpp.timed_out,"elapsed_seconds":cpp.elapsed,"report_status":report_status(cr,cn),"array_count":cn,"report_sha256":sha_bytes(cr)?,"normalized_sha256":sha_bytes(&a)?},"rust":{"command":rust.command,"returncode":rust.code,"signal":rust.signal,"timed_out":rust.timed_out,"elapsed_seconds":rust.elapsed,"report_status":report_status(rr,rn),"array_count":rn,"report_sha256":sha_bytes(rr)?,"normalized_sha256":sha_bytes(&b)?},"case":extra});
    fs::write(
        dir.join("metadata.json"),
        serde_json::to_vec_pretty(&metadata).map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())?;
    Ok(dir)
}
fn report_status(bytes: &[u8], count: Option<u64>) -> &'static str {
    if count.is_some() {
        "ok"
    } else if bytes.is_empty() {
        "missing"
    } else {
        "malformed"
    }
}
fn first_difference(a: &[u8], b: &[u8]) -> Option<serde_json::Value> {
    let aa: Vec<_> = a.split_inclusive(|x| *x == b'\n').collect();
    let bb: Vec<_> = b.split_inclusive(|x| *x == b'\n').collect();
    for i in 0..aa.len().max(bb.len()) {
        if aa.get(i) != bb.get(i) {
            return Some(
                json!({"line":i+1,"expected":aa.get(i).map(|x|String::from_utf8_lossy(x)),"actual":bb.get(i).map(|x|String::from_utf8_lossy(x))}),
            );
        }
    }
    None
}
fn run(
    binary: &Path,
    input: &Path,
    report: &Path,
    timeout: Duration,
    legacy: bool,
) -> Result<Run, String> {
    let mut c = Command::new(binary);
    c.args(["-in"])
        .arg(input)
        .args(["-out"])
        .arg(report)
        .arg("-quiet");
    if legacy {
        c.arg("--legacy-coords");
    }
    let command = std::iter::once(binary.display().to_string())
        .chain(c.get_args().map(|a| a.to_string_lossy().into_owned()))
        .collect();
    let stdout_path = sidecar(report, ".stdout");
    let stderr_path = sidecar(report, ".stderr");
    match spawn_logged(c, &stdout_path, &stderr_path) {
        Ok((child, loggers)) => {
            let (status, timed_out, elapsed) = finish_child(child, timeout);
            finish_loggers(loggers);
            let code = status.as_ref().and_then(|s| s.code());
            #[cfg(unix)]
            let signal = status.as_ref().and_then(|s| s.signal());
            #[cfg(not(unix))]
            let signal = None;
            Ok(Run {
                code,
                signal,
                timed_out,
                stdout: read_log(&stdout_path)?,
                stderr: read_log(&stderr_path)?,
                elapsed,
                command,
            })
        }
        Err(e) => Ok(Run {
            code: Some(127),
            signal: None,
            timed_out: false,
            stdout: vec![],
            stderr: e.into_bytes(),
            elapsed: 0.0,
            command,
        }),
    }
}
fn main() -> std::process::ExitCode {
    match execute() {
        Ok(c) => {
            if c == 0 {
                std::process::ExitCode::SUCCESS
            } else {
                std::process::ExitCode::from(1)
            }
        }
        Err(e) => {
            eprintln!("error: {e}");
            std::process::ExitCode::from(1)
        }
    }
}
fn bounded_status(
    command: Command,
    timeout: Duration,
    stdout_path: &Path,
    stderr_path: &Path,
) -> Result<(), String> {
    let rendered = format!("{:?}", command);
    let (child, loggers) = spawn_logged(command, stdout_path, stderr_path)
        .map_err(|e| format!("start {rendered}: {e}"))?;
    let (status, timed_out, _) = finish_child(child, timeout);
    finish_loggers(loggers);
    let stderr = read_log(stderr_path)?;
    if timed_out {
        return Err(format!(
            "build timed out after {}s: {rendered}: {}",
            timeout.as_secs_f64(),
            String::from_utf8_lossy(&stderr)
        ));
    }
    if status.as_ref().and_then(|s| s.code()) != Some(0) {
        return Err(format!(
            "build failed: {}",
            String::from_utf8_lossy(&stderr)
        ));
    }
    Ok(())
}
struct WorkDir(PathBuf);
impl Drop for WorkDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
struct SummaryGuard {
    path: PathBuf,
    summary: serde_json::Value,
    error: Option<String>,
    done: bool,
}
impl Drop for SummaryGuard {
    fn drop(&mut self) {
        if !self.done {
            self.summary["status"] = "failed".into();
            self.summary["error"] = self
                .error
                .clone()
                .unwrap_or_else(|| "campaign terminated before completion".to_owned())
                .into();
            self.summary["finished_at"] = now().into();
            if let Ok(bytes) = serde_json::to_vec_pretty(&self.summary) {
                let _ = fs::write(&self.path, bytes);
            }
        }
    }
}
fn execute() -> Result<i32, String> {
    let o = parse()?;
    let root = fs::canonicalize(Path::new(env!("CARGO_MANIFEST_DIR")).join("../.."))
        .map_err(|e| e.to_string())?;
    let out = canonical_destination(&o.output)?;
    for p in [
        root.join("ref/pilercr1.06-src"),
        root.join("ref"),
        root.join("scripts"),
        root.join("test-data"),
        root.join("crates"),
        root.join("target"),
    ] {
        // Protected paths such as ref/ or test-data/ may not exist on a fresh
        // checkout; canonical_destination resolves the existing ancestors instead.
        if overlap(&out, &canonical_destination(&p)?) {
            return Err(format!("protected destination: {}", out.display()));
        }
    }
    if out.exists()
        && (!out.is_dir()
            || fs::read_dir(&out)
                .map_err(|e| e.to_string())?
                .next()
                .is_some())
    {
        return Err("output must be absent or empty".into());
    }
    fs::create_dir_all(&out).map_err(|e| e.to_string())?;
    let mut finalizer = SummaryGuard {
        path: out.join("summary.json"),
        summary: json!({"status":"running","started_at":now(),"output":out,"master_seed":o.seed,"target":o.target,"attempts":0,"accepted":0,"mismatches":0,"reference_crashes":0}),
        error: None,
        done: false,
    };
    let result = (|| -> Result<i32, String> {
        fs::create_dir_all(&o.tmpdir).map_err(|e| e.to_string())?;
        let tmp = fs::canonicalize(&o.tmpdir).map_err(|e| e.to_string())?;
        if !tmp.is_dir() {
            return Err(format!(
                "temporary directory unavailable: {}",
                tmp.display()
            ));
        }
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let work = tmp.join(format!("pilercr-campaign-{unique}"));
        fs::create_dir_all(&work).map_err(|e| e.to_string())?;
        let _work_guard = WorkDir(work.clone());
        let reference = if let Some(p) = o.reference {
            fs::canonicalize(&p).map_err(|e| format!("reference binary {}: {e}", p.display()))?
        } else {
            let bin = work.join("reference/pilercr");
            fs::create_dir_all(bin.parent().unwrap()).map_err(|e| e.to_string())?;
            let src = root.join("ref/pilercr1.06-src");
            let mut cmd = Command::new(env::var("CXX").unwrap_or("c++".into()));
            cmd.args(["-std=c++17", "-O3", "-DNDEBUG=1"])
                .args(
                    fs::read_dir(src)
                        .map_err(|e| e.to_string())?
                        .filter_map(Result::ok)
                        .map(|e| e.path())
                        .filter(|p| p.extension().is_some_and(|e| e == "cpp")),
                )
                .arg("-lm")
                .arg("-o")
                .arg(&bin);
            cmd.current_dir(&root);
            bounded_status(
                cmd,
                o.build_timeout,
                &work.join("reference-build.stdout"),
                &work.join("reference-build.stderr"),
            )?;
            bin
        };
        let rust = if let Some(p) = o.rust {
            fs::canonicalize(&p).map_err(|e| format!("Rust binary {}: {e}", p.display()))?
        } else {
            let target = work.join("cargo-target");
            let mut cmd = Command::new("cargo");
            cmd.args(["build", "--release", "--manifest-path"])
                .arg(root.join("Cargo.toml"))
                .args(["-p", "pilercr-rs"])
                .env("CARGO_TARGET_DIR", &target)
                .current_dir(&root);
            bounded_status(
                cmd,
                o.build_timeout,
                &work.join("cargo-build.stdout"),
                &work.join("cargo-build.stderr"),
            )?;
            target.join("release/pilercr-rs")
        };
        let mut summary = json!({"started_at":now(),"environment":{"os":env::consts::OS,"arch":env::consts::ARCH,"CXX":env::var("CXX").ok()},"master_seed":o.seed,"target":o.target,"max_attempts":if o.max_attempts==0{2000usize.max(o.target*5)}else{o.max_attempts},"max_misses":o.max_misses,"max_mutations":o.max_mutations,"mutation_steps":o.mutation_steps,"negative_every":o.negative_every,"output":out,"reference":{"path":reference,"sha256":hash(&reference)?},"rust":{"path":rust,"sha256":hash(&rust)?},"attempts":0,"original_attempts":0,"mutation_attempts":0,"accepted":0,"positives_generated":0,"negatives_generated":0,"negatives_detected":0,"undetected_positives":0,"mismatches":0,"reference_crashes":0,"miss_streak":0,"max_miss_streak":0,"mutation_ladders":0,"mutation_ladders_completed":0,"mutation_ladders_capped":0,"mutations_detected":0,"mutations_undetected":0,"mismatch_cases":[],"reference_crash_cases":[],"status":"running"});
        finalizer.summary = summary.clone();
        let budget = summary["max_attempts"].as_u64().unwrap() as usize;
        let mut attempt = 0;
        let mut streak = 0;
        while (summary["accepted"].as_u64().unwrap() as usize) < o.target
            && attempt < budget
            && !(o.max_misses > 0 && streak >= o.max_misses)
        {
            attempt += 1;
            let neg = o.negative_every > 0 && attempt % o.negative_every == 0;
            let (seed, fasta, meta) = generate_case(o.seed, attempt, neg);
            if neg {
                summary["negatives_generated"] =
                    (summary["negatives_generated"].as_u64().unwrap() + 1).into()
            } else {
                summary["positives_generated"] =
                    (summary["positives_generated"].as_u64().unwrap() + 1).into()
            }
            let tag = format!("case-{attempt:06}");
            let input = work.join(format!("{tag}.fa"));
            fs::write(&input, &fasta).map_err(|e| e.to_string())?;
            let cp = work.join(format!("{tag}.cpp.txt"));
            let rp = work.join(format!("{tag}.rust.txt"));
            let cpp = run(&reference, &input, &cp, o.timeout, false)?;
            let rs = run(&rust, &input, &rp, o.timeout, true)?;
            let cr = fs::read(&cp).unwrap_or_default();
            let rr = fs::read(&rp).unwrap_or_default();
            let cn = validated_array_count(&cr, input.to_str().unwrap());
            let rn = validated_array_count(&rr, input.to_str().unwrap());
            let detected = cn.is_some_and(|n| n > 0) || rn.is_some_and(|n| n > 0);
            let class = classify_signal(
                if cpp.timed_out { None } else { cpp.signal },
                !rs.timed_out && rs.code == Some(0),
                &rr,
                input.to_str().unwrap(),
            );
            let mismatch = cpp.timed_out
                || rs.timed_out
                || cpp.code != Some(0)
                || rs.code != Some(0)
                || cn.is_none()
                || rn.is_none()
                || normalized(&cr, &input) != normalized(&rr, &input);
            summary["attempts"] = (summary["attempts"].as_u64().unwrap() + 1).into();
            summary["original_attempts"] = (attempt as u64).into();
            if neg {
                if detected {
                    summary["negatives_detected"] =
                        (summary["negatives_detected"].as_u64().unwrap() + 1).into()
                }
            } else if detected {
                summary["accepted"] = (summary["accepted"].as_u64().unwrap() + 1).into()
            } else {
                summary["undetected_positives"] =
                    (summary["undetected_positives"].as_u64().unwrap() + 1).into()
            }
            streak = if detected { 0 } else { streak + 1 };
            summary["miss_streak"] = (streak as u64).into();
            summary["max_miss_streak"] = summary["max_miss_streak"]
                .as_u64()
                .unwrap()
                .max(streak as u64)
                .into();
            finalizer.summary = summary.clone();
            if mismatch || class == CRASH_EXEMPT {
                let id = format!("attempt-{attempt:06}-seed-{seed:016x}");
                let case_out = if class == CRASH_EXEMPT {
                    summary["reference_crashes"] =
                        (summary["reference_crashes"].as_u64().unwrap() + 1).into();
                    out.join("reference-crashes")
                } else {
                    summary["mismatches"] = (summary["mismatches"].as_u64().unwrap() + 1).into();
                    out.clone()
                };
                let path = persist_case(
                    &case_out,
                    &id,
                    &fasta,
                    &input,
                    &cpp,
                    &rs,
                    &cr,
                    &rr,
                    cn,
                    rn,
                    json!({"attempt":attempt,"seed":seed,"negative":neg,"classification":class,"generator":{"positive":meta.positive,"kind":meta.kind,"repeat_spans":meta.repeat_spans}}),
                )?;
                let index = if class == CRASH_EXEMPT {
                    "reference_crash_cases"
                } else {
                    "mismatch_cases"
                };
                summary[index]
                    .as_array_mut()
                    .unwrap()
                    .push(json!({"case_id":id,"directory":path}));
                finalizer.summary = summary.clone();
            }
            cleanup_case_artifacts(&input, &[&cp, &rp])?;
            if !neg
                && detected
                && !meta.repeat_spans.is_empty()
                && o.mutation_steps > 0
                && (o.max_mutations == 0
                    || summary["mutation_attempts"].as_u64().unwrap() < (o.max_mutations as u64))
            {
                summary["mutation_ladders"] =
                    (summary["mutation_ladders"].as_u64().unwrap() + 1).into();
                let seq = String::from_utf8_lossy(&fasta)
                    .lines()
                    .skip(1)
                    .collect::<String>();
                let n: usize = meta.repeat_spans.iter().map(|(a, b)| b - a).sum();
                let mut complete = false;
                for (step, count) in mutation_step_substitutions(n, o.mutation_steps)
                    .into_iter()
                    .enumerate()
                {
                    if o.max_mutations > 0
                        && summary["mutation_attempts"].as_u64().unwrap() >= o.max_mutations as u64
                    {
                        break;
                    }
                    let m = mutate_repeats(
                        &seq,
                        &meta.repeat_spans,
                        count,
                        seed.wrapping_add(step as u64),
                    );
                    let mf = format!(">mutation_{attempt:06}_step{:02}\n{m}\n", step + 1);
                    let mi = work.join(format!("m-{attempt}-{step}.fa"));
                    let mr = work.join(format!("m-{attempt}-{step}.out"));
                    let mf = mf.into_bytes();
                    fs::write(&mi, &mf).map_err(|e| e.to_string())?;
                    let cpp_mut = run(&reference, &mi, &mr, o.timeout, false)?;
                    let cr = fs::read(&mr).unwrap_or_default();
                    let rustout = work.join(format!("m-{attempt}-{step}.rust"));
                    let rust_mut = run(&rust, &mi, &rustout, o.timeout, true)?;
                    let rr = fs::read(&rustout).unwrap_or_default();
                    let cn = validated_array_count(&cr, mi.to_str().unwrap());
                    let rn = validated_array_count(&rr, mi.to_str().unwrap());
                    let det = cn.is_some_and(|x| x > 0) || rn.is_some_and(|x| x > 0);
                    let class = classify_signal(
                        if cpp_mut.timed_out {
                            None
                        } else {
                            cpp_mut.signal
                        },
                        !rust_mut.timed_out && rust_mut.code == Some(0),
                        &rr,
                        mi.to_str().unwrap(),
                    );
                    let mismatch = cpp_mut.timed_out
                        || rust_mut.timed_out
                        || cpp_mut.code != Some(0)
                        || rust_mut.code != Some(0)
                        || cn.is_none()
                        || rn.is_none()
                        || normalized(&cr, &mi) != normalized(&rr, &mi);
                    summary["attempts"] = (summary["attempts"].as_u64().unwrap() + 1).into();
                    summary["mutation_attempts"] =
                        (summary["mutation_attempts"].as_u64().unwrap() + 1).into();
                    let key = if det {
                        "mutations_detected"
                    } else {
                        "mutations_undetected"
                    };
                    summary[key] = (summary[key].as_u64().unwrap() + 1).into();
                    finalizer.summary = summary.clone();
                    if mismatch || class == CRASH_EXEMPT {
                        let mseed = seed.wrapping_add(step as u64);
                        let id = format!(
                            "mutation-{attempt:06}-step-{:02}-seed-{mseed:016x}",
                            step + 1
                        );
                        let base = if class == CRASH_EXEMPT {
                            summary["reference_crashes"] =
                                (summary["reference_crashes"].as_u64().unwrap() + 1).into();
                            out.join("reference-crashes")
                        } else {
                            summary["mismatches"] =
                                (summary["mismatches"].as_u64().unwrap() + 1).into();
                            out.clone()
                        };
                        let path = persist_case(
                            &base,
                            &id,
                            &mf,
                            &mi,
                            &cpp_mut,
                            &rust_mut,
                            &cr,
                            &rr,
                            cn,
                            rn,
                            json!({"classification":class,"generator":{"mutation":true,"parent_attempt":attempt,"parent_seed":seed,"mutation_seed":mseed,"step":step+1,"substitutions":count,"repeat_positions":n}}),
                        )?;
                        let ix = if class == CRASH_EXEMPT {
                            "reference_crash_cases"
                        } else {
                            "mismatch_cases"
                        };
                        summary[ix]
                            .as_array_mut()
                            .unwrap()
                            .push(json!({"case_id":id,"directory":path}));
                        finalizer.summary = summary.clone();
                    }
                    cleanup_case_artifacts(&mi, &[&mr, &rustout])?;
                    if !det {
                        complete = true;
                        break;
                    }
                }
                let key = if complete {
                    "mutation_ladders_completed"
                } else {
                    "mutation_ladders_capped"
                };
                summary[key] = (summary[key].as_u64().unwrap() + 1).into();
                finalizer.summary = summary.clone();
            }
        }
        summary["status"] = if summary["accepted"].as_u64().unwrap() as usize >= o.target {
            "complete"
        } else {
            "target-unmet"
        }
        .into();
        summary["stop_reason"] = if summary["status"] == "complete" {
            "target"
        } else if streak >= o.max_misses && o.max_misses > 0 {
            "max-misses"
        } else {
            "max-attempts"
        }
        .into();
        summary["finished_at"] = now().into();
        finalizer.summary = summary.clone();
        let text = serde_json::to_vec_pretty(&summary).map_err(|e| e.to_string())?;
        fs::write(out.join("summary.json"), text).map_err(|e| e.to_string())?;
        finalizer.done = true;
        println!(
            "status={} accepted={}/{} attempts={} mismatches={} reference_crashes={} output={}",
            summary["status"],
            summary["accepted"],
            o.target,
            summary["attempts"],
            summary["mismatches"],
            summary["reference_crashes"],
            out.display()
        );
        Ok(if summary["status"] == "complete" {
            0
        } else {
            1
        })
    })();
    if let Err(error) = &result {
        finalizer.error = Some(error.clone());
    }
    result
}

#[cfg(test)]
mod tests {
    use super::{Run, normalized, persist_case};
    use serde_json::Value;
    use std::{
        fs,
        path::Path,
        time::{SystemTime, UNIX_EPOCH},
    };
    #[test]
    fn normalization_changes_only_the_valid_current_input_header() {
        let input = Path::new("/tmp/current.fa");
        let report=b"/old/current.fa: 0 putative CRISPR arrays found.\nbody marker: not a header putative CRISPR arrays found.\n";
        assert_eq!(normalized(report,input),b"<INPUT>: 0 putative CRISPR arrays found.\nbody marker: not a header putative CRISPR arrays found.\n");
        let other = b"/old/other.fa: 0 putative CRISPR arrays found.\n";
        assert_eq!(normalized(other, input), other);
    }
    #[test]
    fn persisted_normalized_hashes_use_actual_input_filename() {
        let root = std::env::temp_dir().join(format!(
            "campaign-normalized-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&root).unwrap();
        let input = root.join("same.fa");
        let report_a = b"/first/path/same.fa: 0 putative CRISPR arrays found.\nidentical body\n";
        let report_b = b"/elsewhere/same.fa: 0 putative CRISPR arrays found.\nidentical body\n";
        let run = |command: Vec<String>| Run {
            code: Some(0),
            signal: None,
            timed_out: false,
            stdout: vec![],
            stderr: vec![],
            elapsed: 0.01,
            command,
        };
        let dir = persist_case(
            &root,
            "case",
            b">same\nACGT\n",
            &input,
            &run(vec!["cpp".into()]),
            &run(vec!["rust".into()]),
            report_a,
            report_b,
            Some(0),
            Some(0),
            Value::Null,
        )
        .unwrap();
        let metadata: Value =
            serde_json::from_slice(&fs::read(dir.join("metadata.json")).unwrap()).unwrap();
        assert_eq!(
            metadata["cpp"]["normalized_sha256"],
            metadata["rust"]["normalized_sha256"]
        );
        assert!(metadata["first_difference"].is_null());
        fs::remove_dir_all(root).unwrap();
    }
}
