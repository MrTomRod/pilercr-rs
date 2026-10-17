use serde_json::Value;
use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
    time::{Duration, SystemTime, UNIX_EPOCH},
};

struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "campaign-test-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&p).unwrap();
        Self(p)
    }
    fn path(&self) -> &Path {
        &self.0
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
struct Cleanup(PathBuf);
impl Drop for Cleanup {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
struct EscapedCleanup(PathBuf);
impl Drop for EscapedCleanup {
    fn drop(&mut self) {
        if let Ok(pid) = fs::read_to_string(&self.0)
            && let Ok(pid) = pid.trim().parse::<u32>()
        {
            let _ = Command::new("kill")
                .args(["-KILL", "--", &format!("-{pid}")])
                .stderr(std::process::Stdio::null())
                .status();
            let _ = Command::new("kill")
                .args(["-KILL", &pid.to_string()])
                .stderr(std::process::Stdio::null())
                .status();
        }
    }
}
fn tool(dir: &Path, name: &str, body: &str) -> PathBuf {
    let p = dir.join(name);
    fs::write(&p, format!("#!/bin/sh\n{body}\n")).unwrap();
    use std::os::unix::fs::PermissionsExt;
    fs::set_permissions(&p, fs::Permissions::from_mode(0o755)).unwrap();
    p
}
fn invoke(
    tmp: &Temp,
    out: &Path,
    reference: &Path,
    rust: &Path,
    more: &[&str],
) -> std::process::Output {
    let mut c = Command::new(env!("CARGO_BIN_EXE_pilercr-campaign"));
    c.args([
        "--target",
        "1",
        "--max-attempts",
        "1",
        "--negative-every",
        "0",
        "--mutation-steps",
        "0",
        "--output",
    ])
    .arg(out)
    .args(["--tmpdir"])
    .arg(tmp.path())
    .args(["--reference-binary"])
    .arg(reference)
    .args(["--rust-binary"])
    .arg(rust)
    .args(more);
    c.output().unwrap()
}
fn scripts(dir: &Path, cpp: &str, rust: &str) -> (PathBuf, PathBuf) {
    let body = |count: &str, tail: &str| {
        format!(
            "src=; out=; while [ $# -gt 0 ]; do case $1 in -in) src=$2; shift 2;; -out) out=$2; shift 2;; *) shift;; esac; done; printf '%s: {count} putative CRISPR arrays found.\\n' \"$src\" > \"$out\"; {tail}"
        )
    };
    (
        tool(dir, "ref", &body(cpp, "")),
        tool(dir, "rust", &body(rust, "")),
    )
}

#[test]
fn real_signal_is_exempt_but_exit_one_is_a_mismatch() {
    let t = Temp::new();
    let (r, u) = scripts(t.path(), "0", "0");
    let (signal, ordinary) = (
        tool(t.path(), "signal", "kill -ABRT $$"),
        tool(t.path(), "exit-one", "exit 1"),
    );
    for (bin, expected) in [(signal, "reference-crash-exempt"), (ordinary, "comparable")] {
        let out = t.path().join(expected);
        let got = invoke(&t, &out, &bin, &u, &[]);
        assert_eq!(got.status.code(), Some(1));
        let summary: Value =
            serde_json::from_slice(&fs::read(out.join("summary.json")).unwrap()).unwrap();
        if expected == "reference-crash-exempt" {
            assert_eq!(summary["reference_crashes"], 1);
            assert_eq!(
                summary["reference_crash_cases"].as_array().unwrap().len(),
                1
            )
        } else {
            assert_eq!(summary["mismatches"], 1)
        }
        let _ = &r;
    }
}

#[test]
fn mutation_mismatch_persists_provenance_and_index() {
    let t = Temp::new();
    let reference = tool(
        t.path(),
        "ref",
        "src=; out=; while [ $# -gt 0 ]; do case $1 in -in) src=$2; shift 2;; -out) out=$2; shift 2;; *) shift;; esac; done; printf '%s: 1 putative CRISPR arrays found.\\n' \"$src\" > \"$out\"",
    );
    let rust = tool(
        t.path(),
        "rust",
        "src=; out=; while [ $# -gt 0 ]; do case $1 in -in) src=$2; shift 2;; -out) out=$2; shift 2;; *) shift;; esac; done; printf '%s: 2 putative CRISPR arrays found.\\n' \"$src\" > \"$out\"",
    );
    let out = t.path().join("out");
    let result = Command::new(env!("CARGO_BIN_EXE_pilercr-campaign"))
        .args([
            "--target",
            "1",
            "--max-attempts",
            "1",
            "--negative-every",
            "0",
            "--mutation-steps",
            "1",
            "--output",
        ])
        .arg(&out)
        .args(["--tmpdir"])
        .arg(t.path())
        .args(["--reference-binary"])
        .arg(reference)
        .args(["--rust-binary"])
        .arg(rust)
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(0));
    let s: Value = serde_json::from_slice(&fs::read(out.join("summary.json")).unwrap()).unwrap();
    assert!(s["mismatches"].as_u64().unwrap() >= 2);
    let entry = s["mismatch_cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|x| x["case_id"].as_str().unwrap().starts_with("mutation-"))
        .unwrap();
    let meta: Value = serde_json::from_slice(
        &fs::read(Path::new(entry["directory"].as_str().unwrap()).join("metadata.json")).unwrap(),
    )
    .unwrap();
    assert!(meta["case"]["generator"]["parent_attempt"].is_number());
    assert!(meta["case"]["generator"]["mutation_seed"].is_number());
    assert!(meta["cpp"]["input_sha256"].is_null());
    assert!(meta["input_sha256"].is_string());
}

#[test]
fn symlink_output_to_protected_tree_is_rejected() {
    let t = Temp::new();
    let alias = t.path().join("alias");
    std::os::unix::fs::symlink(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target"),
        &alias,
    )
    .unwrap();
    let (r, u) = scripts(t.path(), "0", "0");
    let result = invoke(&t, &alias.join("nested"), &r, &u, &[]);
    assert_eq!(result.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&result.stderr).contains("protected destination"));
    let dotdot = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/../target/child");
    let result = invoke(&t, &dotdot, &r, &u, &[]);
    assert_eq!(result.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&result.stderr).contains("protected destination"));
}

#[test]
fn rejected_outputs_are_never_rewritten() {
    let t = Temp::new();
    let (r, u) = scripts(t.path(), "0", "0");
    let nonempty = t.path().join("occupied");
    fs::create_dir(&nonempty).unwrap();
    let original = b"keep-this-summary-byte-for-byte\n";
    fs::write(nonempty.join("summary.json"), original).unwrap();
    let result = invoke(&t, &nonempty, &r, &u, &[]);
    assert_eq!(result.status.code(), Some(1));
    assert_eq!(fs::read(nonempty.join("summary.json")).unwrap(), original);
    let protected = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../target")
        .join(format!(
            "campaign-protected-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
    fs::create_dir_all(&protected).unwrap();
    let _cleanup = Cleanup(protected.clone());
    let summary = protected.join("summary.json");
    fs::write(&summary, original).unwrap();
    let result = invoke(&t, &protected, &r, &u, &[]);
    assert_eq!(result.status.code(), Some(1));
    assert_eq!(fs::read(summary).unwrap(), original);
}

#[test]
fn failed_summary_keeps_progress_after_a_persisted_case() {
    let t = Temp::new();
    let (r, u) = scripts(t.path(), "1", "2");
    let out = t.path().join("failed-after-case");
    let state = t.path().join("sha-count");
    let _sha = tool(
        t.path(),
        "sha256sum",
        &format!(
            "n=$(cat {} 2>/dev/null || echo 0); n=$((n+1)); printf '%s' \"$n\" > {}; if [ \"$n\" -ge 8 ]; then echo injected-hash-failure >&2; exit 1; fi; exec /usr/bin/sha256sum \"$@\"",
            state.display(),
            state.display()
        ),
    );
    let result = Command::new(env!("CARGO_BIN_EXE_pilercr-campaign"))
        .args([
            "--target",
            "2",
            "--max-attempts",
            "2",
            "--negative-every",
            "0",
            "--mutation-steps",
            "0",
            "--seed",
            "17679207994573498142",
            "--output",
        ])
        .arg(&out)
        .args(["--tmpdir"])
        .arg(t.path())
        .args(["--reference-binary"])
        .arg(r)
        .args(["--rust-binary"])
        .arg(u)
        .env("PATH", format!("{}:/usr/bin:/bin", t.path().display()))
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(1));
    let s: Value = serde_json::from_slice(&fs::read(out.join("summary.json")).unwrap()).unwrap();
    assert_eq!(s["status"], "failed");
    assert_eq!(s["attempts"], 2);
    assert_eq!(s["mismatches"], 1);
    assert_eq!(s["mismatch_cases"].as_array().unwrap().len(), 1);
    assert_eq!(s["max_attempts"], 2);
    assert!(s["error"].as_str().unwrap().contains("sha256sum failed"));
}

#[test]
fn escaped_pipe_descendant_cannot_hold_campaign_past_deadline() {
    let t = Temp::new();
    let (r, _) = scripts(t.path(), "0", "0");
    let pidfile = t.path().join("escaped.pid");
    let _escaped_cleanup = EscapedCleanup(pidfile.clone());
    let rust = tool(
        t.path(),
        "rust",
        &format!(
            "printf 'partial-evidence\\n' >&2; setsid sh -c 'echo $$ > {} ; sleep 10' & sleep 10",
            pidfile.display()
        ),
    );
    let out = t.path().join("escaped-output");
    let started = std::time::Instant::now();
    let mut c = Command::new("timeout");
    c.args([
        "4",
        env!("CARGO_BIN_EXE_pilercr-campaign"),
        "--target",
        "1",
        "--max-attempts",
        "1",
        "--negative-every",
        "0",
        "--mutation-steps",
        "0",
        "--output",
    ])
    .arg(&out)
    .args(["--tmpdir"])
    .arg(t.path())
    .args(["--reference-binary"])
    .arg(&r)
    .args(["--rust-binary"])
    .arg(&rust)
    .args(["--timeout", "0.5"]);
    let result = c.output().unwrap();
    let elapsed = started.elapsed();
    assert!(
        elapsed < Duration::from_secs(4),
        "campaign blocked for {elapsed:?}"
    );
    assert_ne!(
        result.status.code(),
        Some(124),
        "external watchdog had to kill campaign"
    );
    assert_eq!(result.status.code(), Some(1));
    let s: Value = serde_json::from_slice(&fs::read(out.join("summary.json")).unwrap()).unwrap();
    let index = s["mismatch_cases"].as_array().unwrap();
    let dir = Path::new(index[0]["directory"].as_str().unwrap());
    let metadata: Value =
        serde_json::from_slice(&fs::read(dir.join("metadata.json")).unwrap()).unwrap();
    assert_eq!(metadata["rust"]["timed_out"], true);
    assert!(dir.join("rust.stderr").is_file());
    assert!(dir.join("rust.stdout").is_file());
}

#[test]
fn cli_generation_terminates_for_zero_mapped_seed() {
    let t = Temp::new();
    let (r, u) = scripts(t.path(), "1", "1");
    let out = t.path().join("zero-state");
    let result = Command::new(env!("CARGO_BIN_EXE_pilercr-campaign"))
        .args([
            "--target",
            "1",
            "--max-attempts",
            "1",
            "--negative-every",
            "0",
            "--mutation-steps",
            "0",
            "--seed",
            "17679207994573498142",
            "--output",
        ])
        .arg(&out)
        .args(["--tmpdir"])
        .arg(t.path())
        .args(["--reference-binary"])
        .arg(r)
        .args(["--rust-binary"])
        .arg(u)
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(0));
}

#[test]
fn unrepresentable_cli_timeouts_are_rejected_before_spawn() {
    let t = Temp::new();
    for (flag, label) in [("--timeout", "tool"), ("--build-timeout", "build")] {
        let marker = t.path().join(format!("{label}-started"));
        let cxx = tool(
            t.path(),
            &format!("cxx-{label}"),
            &format!("touch {}", marker.display()),
        );
        let output = t.path().join(format!("output-{label}"));
        let result = Command::new(env!("CARGO_BIN_EXE_pilercr-campaign"))
            .args([flag, "1e100", "--output"])
            .arg(&output)
            .args(["--tmpdir"])
            .arg(t.path())
            .env("CXX", cxx)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(1));
        assert!(!marker.exists(), "child unexpectedly started for {flag}");
        assert!(
            !output.exists(),
            "output should not be owned for invalid {flag}"
        );
        assert!(String::from_utf8_lossy(&result.stderr).contains("out of range"));
    }
}

#[test]
fn build_timeout_writes_failure_and_cleans_workdir() {
    let t = Temp::new();
    let cxx = tool(t.path(), "slow-cxx", "sleep 10");
    let out = t.path().join("out");
    let result = Command::new(env!("CARGO_BIN_EXE_pilercr-campaign"))
        .args(["--target", "1", "--output"])
        .arg(&out)
        .args(["--tmpdir"])
        .arg(t.path())
        .args(["--build-timeout", "1"])
        .env("CXX", cxx)
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(1));
    let s: Value = serde_json::from_slice(&fs::read(out.join("summary.json")).unwrap()).unwrap();
    assert_eq!(s["status"], "failed");
    assert!(s["error"].as_str().unwrap().contains("timed out"));
    assert_eq!(
        fs::read_dir(t.path())
            .unwrap()
            .filter_map(Result::ok)
            .filter(|e| e
                .file_name()
                .to_string_lossy()
                .starts_with("pilercr-campaign-"))
            .count(),
        0
    );
}

#[test]
fn cargo_build_timeout_is_enforced() {
    let t = Temp::new();
    let cargo = tool(t.path(), "cargo", "sleep 10");
    let reference = tool(t.path(), "reference", "exit 0");
    let out = t.path().join("cargo-out");
    let result = Command::new(env!("CARGO_BIN_EXE_pilercr-campaign"))
        .args(["--target", "1", "--output"])
        .arg(&out)
        .args(["--tmpdir"])
        .arg(t.path())
        .args(["--reference-binary"])
        .arg(reference)
        .args(["--build-timeout", "1"])
        .env("PATH", format!("{}:/usr/bin:/bin", t.path().display()))
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(1));
    let summary: Value =
        serde_json::from_slice(&fs::read(out.join("summary.json")).unwrap()).unwrap();
    assert_eq!(summary["status"], "failed");
    assert!(summary["error"].as_str().unwrap().contains("timed out"));
    let _ = cargo;
}

#[test]
fn log_cap_does_not_break_the_tool_pipe() {
    let t = Temp::new();
    let r = tool(
        t.path(),
        "ref",
        "src=; out=; while [ $# -gt 0 ]; do case $1 in -in) src=$2; shift 2;; -out) out=$2; shift 2;; *) shift;; esac; done; head -c 3145729 /dev/zero >&2 || exit 3; printf '%s: 1 putative CRISPR arrays found.\\n' \"$src\" > \"$out\"",
    );
    let u = tool(
        t.path(),
        "rust",
        "src=; out=; while [ $# -gt 0 ]; do case $1 in -in) src=$2; shift 2;; -out) out=$2; shift 2;; *) shift;; esac; done; printf '%s: 2 putative CRISPR arrays found.\\n' \"$src\" > \"$out\"",
    );
    let out = t.path().join("out");
    let result = invoke(&t, &out, &r, &u, &[]);
    assert_eq!(result.status.code(), Some(0));
    let s: Value = serde_json::from_slice(&fs::read(out.join("summary.json")).unwrap()).unwrap();
    assert_eq!(s["mismatches"], 1);
    let dir = Path::new(s["mismatch_cases"][0]["directory"].as_str().unwrap());
    let metadata: Value =
        serde_json::from_slice(&fs::read(dir.join("metadata.json")).unwrap()).unwrap();
    assert_eq!(metadata["cpp"]["returncode"], 0);
    assert!(metadata["cpp"]["signal"].is_null());
    assert_eq!(metadata["cpp"]["timed_out"], false);
    assert_eq!(metadata["cpp"]["report_status"], "ok");
    let stderr_len = fs::metadata(dir.join("cpp.stderr")).unwrap().len();
    assert!(stderr_len > 0);
    assert!(stderr_len <= 2 * 1024 * 1024);
}

#[test]
fn work_files_are_removed_per_case_and_evidence_is_kept() {
    let t = Temp::new();
    let log = t.path().join("ls-log");
    let r = tool(
        t.path(),
        "ref",
        &format!(
            "src=; out=; while [ $# -gt 0 ]; do case $1 in -in) src=$2; shift 2;; -out) out=$2; shift 2;; *) shift;; esac; done; d=${{out%/*}}; ls -A \"$d\" | wc -l >> {}; printf '%s: 1 putative CRISPR arrays found.\\n' \"$src\" > \"$out\"",
            log.display()
        ),
    );
    let u = tool(
        t.path(),
        "rust",
        "src=; out=; while [ $# -gt 0 ]; do case $1 in -in) src=$2; shift 2;; -out) out=$2; shift 2;; *) shift;; esac; done; printf '%s: 2 putative CRISPR arrays found.\\n' \"$src\" > \"$out\"",
    );
    let out = t.path().join("out");
    let result = Command::new(env!("CARGO_BIN_EXE_pilercr-campaign"))
        .args([
            "--target",
            "3",
            "--max-attempts",
            "3",
            "--negative-every",
            "0",
            "--mutation-steps",
            "0",
            "--output",
        ])
        .arg(&out)
        .args(["--tmpdir"])
        .arg(t.path())
        .args(["--reference-binary"])
        .arg(&r)
        .args(["--rust-binary"])
        .arg(&u)
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(0));
    let s: Value = serde_json::from_slice(&fs::read(out.join("summary.json")).unwrap()).unwrap();
    assert_eq!(s["attempts"], 3);
    assert_eq!(s["mismatches"], 3);
    for entry in s["mismatch_cases"].as_array().unwrap() {
        let dir = Path::new(entry["directory"].as_str().unwrap());
        assert!(dir.join("input.fa").is_file());
        assert!(dir.join("cpp.stdout").is_file());
        assert!(dir.join("metadata.json").is_file());
    }
    let counts: Vec<String> = fs::read_to_string(&log)
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect();
    assert_eq!(counts.len(), 3);
    assert!(counts.iter().all(|c| c == "3"), "work dir grew: {counts:?}");
    let leftovers = fs::read_dir(t.path())
        .unwrap()
        .filter_map(Result::ok)
        .any(|e| {
            e.file_name()
                .to_string_lossy()
                .starts_with("pilercr-campaign-")
        });
    assert!(!leftovers);
}

#[test]
fn timeout_and_large_stderr_are_bounded_and_retained() {
    let t = Temp::new();
    let r = tool(
        t.path(),
        "ref",
        "src=; out=; while [ $# -gt 0 ]; do case $1 in -in) src=$2; shift 2;; -out) out=$2; shift 2;; *) shift;; esac; done; yes diagnostic | head -c 1000000 >&2; printf '%s: 0 putative CRISPR arrays found.\\n' \"$src\" > \"$out\"",
    );
    let u = tool(t.path(), "rust", "printf timeout-evidence >&2; sleep 5");
    let out = t.path().join("out");
    let result = invoke(&t, &out, &r, &u, &["--timeout", "1"]);
    assert_eq!(result.status.code(), Some(1));
    let s: Value = serde_json::from_slice(&fs::read(out.join("summary.json")).unwrap()).unwrap();
    let entry = &s["mismatch_cases"][0];
    let dir = Path::new(entry["directory"].as_str().unwrap());
    let stderr_len = fs::metadata(dir.join("cpp.stderr")).unwrap().len();
    assert!(stderr_len > 500_000);
    assert!(stderr_len <= 2 * 1024 * 1024);
    let metadata: Value =
        serde_json::from_slice(&fs::read(dir.join("metadata.json")).unwrap()).unwrap();
    assert_eq!(metadata["rust"]["timed_out"], true);
    assert!(
        String::from_utf8_lossy(&fs::read(dir.join("rust.stderr")).unwrap())
            .contains("timeout-evidence")
    );
}
