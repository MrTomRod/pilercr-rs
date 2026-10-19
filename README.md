# pilercr-rs

A Rust reimplementation of [PILER-CR](https://www.drive5.com/pilercr/), the
CRISPR repeat finder by Robert C. Edgar. The goal is to reproduce the
reference program's output exactly, and to add structured outputs and a
reusable library API on top.

PILER-CR is described in: Edgar, R.C. (2007) *PILER-CR: fast and accurate
identification of CRISPR repeats*, BMC Bioinformatics 8:18.

## Status

- The detection pipeline is complete. The independent per-trapezoid DP work
  runs in parallel with Rayon; results, coverage and report ordering are
  deterministic across thread-pool sizes.
- By default, text reports and consensus FASTA are byte-identical to the
  modern local C++ rebuild of PILER-CR v1.06 on all 52 assemblies in the
  verification set (see [Verification](#verification)).
- `--legacy` (library: `Config::sort_mode = SortMode::Legacy`) switches the
  two tie-sensitive sorts to the reverse-engineered GCC 4.4.3 `std::sort`
  order used by the shipped 2011 executable. This is opt-in and does not
  change the default. It is not universal parity: this implementation
  reproduces 51 of the 52 locally available 2011 text reports byte-for-byte;
  assembly 20 still differs (see [Known differences](#known-differences-from-the-reference)).
- JSON, GFF3 and TSV exports are implemented. They use corrected
  coordinates.
- The DP stage is parallel; the remaining stages are single-threaded.
  See [Performance](#performance).
- Results are checked on one platform (Linux x86-64). Bitwise equality on
  other platforms or compilers is not claimed.

## How this code was written

pilercr-rs was written by autonomous LLM coding agents, directed and reviewed
by a human maintainer (Thomas Roder). The agents wrote the Rust code, the
tests, the differential campaign harness and the documentation. Each change
was checked by running the reference comparisons and the test suite, not by
reading the code alone. The code has not had an independent line-by-line
human audit. Treat the correctness claims below as the results of those tests,
and report any divergence from the reference as a bug.

## Reference and provenance

- **Upstream project:** <https://www.drive5.com/pilercr/>
- **Reference source:** `ref/pilercr1.06-src/` is the unmodified PILER-CR
  1.06 C++ source from
  <https://www.drive5.com/pilercr/pilercr1.06.tar.gz> (SHA256
  `50175f7aa171674cda5ba255631f340f9cc7f80e8cc25135a4cb857147d91068`).
  `ref/README.md` gives the provenance and `ref/SHA256SUMS` lists every file.
  Check them with `bash scripts/verify_reference_source.sh`.
- **Trace hooks:** the trace tooling needs extra hooks. They are kept in
  `patches/pilercr-trace.patch` and applied only to a build copy, by
  `scripts/build_reference_trace.sh`. The committed reference tree is pristine.
- **Upstream license:** public domain (`ref/pilercr1.06-src/LICENSE.md`; the
  drive5.com PILER-CR page also describes it as public domain).
- **Excerpt fixtures:** `crates/pilercr-rs/tests/fixtures/excerpts/` holds
  15 anonymized windows from proprietary verification assemblies. Each window
  is one or more arrays plus 200 bp of flank on each side, merged when they
  overlap. Source names, labels and coordinates are not included. The reports
  and consensus FASTA are goldens from the pristine reference.
- **This project:** MIT licensed, see [LICENSE](LICENSE).
- **`pilercr-parser`:** the MIT-licensed crate from crates.io (1.1.0) for
  parsing PILER-CR text reports. It is used only as a dev-dependency to check
  coordinate recovery.

## Workspace layout

| Path | Purpose |
| --- | --- |
| `crates/pilercr-rs` | Library (`pilercr_rs`) and command-line tool `pilercr-rs`: FASTA loading, detection, reports, JSON/GFF3/TSV, consensus FASTA. Tests and committed excerpt fixtures are in `tests/`. |
| `crates/pilercr-campaign` | Differential test harness: generates synthetic inputs and compares this tool with the C++ reference. Dev-only, not published. |
| `ref/` | Pristine PILER-CR 1.06 C++ source, with `README.md` and `SHA256SUMS`. |
| `patches/` | `pilercr-trace.patch`: trace hooks, applied to a build copy only. |
| `test-data/` | Full verification assemblies and oracle reports. Local only, not in the repository. |
| `scripts/` | Build the trace-enabled reference, verify the reference source checksums, fetch the original 2011 executable, and compare trace files. |
| `scratch/` | Local temporary and work files (ignored by git). |

## Building

Requirements:

- A recent stable Rust toolchain (edition 2024).
- A C++17 compiler (`c++` or `$CXX`) only if you want to build the reference
  for comparisons (`scripts/build_reference_trace.sh`) or run the campaign.

```sh
cargo build --release
```

The binary is `target/release/pilercr-rs`.

## Command-line usage

```sh
# Text report and consensus FASTA
pilercr-rs -in genome.fasta -out report.txt -seq consensus.fasta

# Structured outputs (corrected coordinates)
pilercr-rs -in genome.fasta -out report.txt \
  --json arrays.json --gff arrays.gff3 --tsv arrays.tsv

# Several inputs are treated as separate contigs, in order
pilercr-rs -in a.fasta -in b.fasta -out report.txt
```

Options (abridged; run `pilercr-rs --help` for the full list):

| Option | Meaning |
| --- | --- |
| `-in PATH` | FASTA input. Repeat to combine several inputs. |
| `-out PATH` | Text report. |
| `-seq PATH` | Consensus FASTA. |
| `-trimseqs` | Trim similar consensus sequences (requires `-seq`). |
| `--json`, `--gff`, `--tsv PATH` | Structured output with corrected coordinates. |
| `--legacy` | Use the GCC 4.4.3 `std::sort` tie order of the shipped 2011 executable for hit and image sorting. Opt-in; the default remains the modern sort. |
| `--legacy-coords` | Reproduce the reference's gap-shifted coordinates in the text report. Independent of `--legacy`. |
| `--trace-dir DIR` | Write checkpoint TSVs for debugging, into a new or empty directory. |
| `-quiet`, `-noinfo` | Suppress progress messages or the explanatory text. |
| `-minrepeat`, `-maxrepeat`, `-minspacer`, `-maxspacer`, `-minarray` | Repeat, spacer and array size limits (defaults 16, 64, 8, 64, 3). |
| `-minid`, `-mincons`, `-minrepeatratio`, `-minspacerratio`, `-minhitlength` | Similarity and length thresholds. |
| `-diameter`, `-maxmem` | Filter search diameter and parameter-selection memory budget. |

### Coordinates

The reference text report advances each repeat's position by the length of
the gapped alignment row, so coordinates drift after any repeat that
contains a gap. This affects only the text report. The bug was identified
and documented by Jim Rybarski in
[pilercr-parser](https://github.com/jimrybarski/pilercr-parser), which
recovers the true coordinates from legacy reports.

This project does not port that code. Corrected coordinates are computed
directly during report generation, from ungapped lengths. The parser is
used in the test suite to check that the two agree.

- By default, text and structured outputs use **corrected** coordinates,
  which count only real bases.
- `--legacy-coords` reproduces the reference's gap-shifted text
  coordinates. It is separate from `--legacy`, which changes detection
  tie-breaking, not coordinate formatting.

## Library usage

```rust
use pilercr_rs::{
    Config, Sequence, detect,
    report::{self, CoordinateMode},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sequence = Sequence::from_path("genome.fasta")?;
    let arrays = detect(&sequence, &Config::default())?;
    let text = report::text(
        &sequence,
        &arrays,
        "genome.fasta",
        CoordinateMode::Corrected,
        false,
    )?;
    print!("{text}");
    Ok(())
}
```

`detect` returns structured `Array` values. Rendering is separate, so you can
produce any output format from the same results. `Config::sort_mode` selects
`SortMode::Modern` (the default) or `SortMode::Legacy`, the opt-in GCC 4.4.3
sort order used by the CLI's `--legacy` flag.

## Verification

Correctness is checked against the reference, not only against the
author's published reports.

- **Two baselines.** The original 2011 executable is not in the repository.
  `bash scripts/fetch_reference_binary.sh` downloads the upstream tarball,
  verifies its published SHA256, and extracts the prebuilt executable into
  `target/reference-binary/pilercr` for optional diagnostic runs. `cargo test`
  never downloads it.
  Its shipped reports (local `test-data/ref-output/`) differ from a modern
  rebuild on assemblies 9, 16 and 17, because the reference's `std::sort` tie
  order depends on the compiler and standard library. The acceptance target is
  the modern rebuild of the pristine source, which matches the local oracle on
  all 52 assemblies.
- **Sort modes.** The default detection pipeline uses the modern-sort order
  and matches the modern rebuild/local oracle on all 52 assemblies.
  `--legacy` selects a port of GCC 4.4.3's libstdc++ `std::sort` for both
  tie-sensitive sites (DP hits, then image endpoints); the algorithm and its
  equal-key permutations are covered by unit tests. This reproduces the
  shipped 2011 executable's tie order, but not all of its output: measured
  against the locally available 2011 reports, legacy mode is byte-identical on 51/52
  assemblies. The ignored test
  `all_52_legacy_text_matches_original_2011_except_known_divergences` pins
  the one remaining report (20).
- **Committed parity tests (always run).** `reference_parity.rs` compares
  arrays, text reports and consensus FASTA for each excerpt against its
  committed golden. `report.rs` round-trips the same goldens, and
  `parser_integration.rs` checks coordinates from the excerpt reports with
  `pilercr-parser`.
- **Full local verification (`#[ignore]`).** The 52-assembly parity and parser
  tests read `test-data/` and are ignored by default. They run only with
  `--ignored`.
- **Differential campaign.** `pilercr-campaign` generates seeded synthetic
  arrays and negatives, runs both tools, and records mismatches and reference
  crashes. Example:

  ```sh
  cargo run --release -p pilercr-campaign -- \
    --target 50 --mutation-steps 4 \
    --reference-binary target/reference-trace/pilercr \
    --rust-binary target/release/pilercr-rs
  ```

  By default it writes to `scratch/synthetic-campaign-out` and uses
  `scratch/` for temporary files. Run `--help` for all options. Rust's
  generator produces a different random stream than the earlier Python
  harness, so results are not comparable across the two.

### Running the tests

```sh
# Default: committed excerpts and goldens only. No test-data/ or C++ build needed.
cargo test --workspace

# Full local verification: needs test-data/ (assemblies and local-oracle/).
# The reference build is optional and enables the consensus FASTA comparison.
bash scripts/build_reference_trace.sh
PILERCR_ORACLE_DIR="$PWD/test-data/local-oracle" \
PILERCR_REFERENCE_BINARY="$PWD/target/reference-trace/pilercr" \
cargo test --release -p pilercr-rs -- --ignored
```

The ignored set includes a differential test against `g++` (`filter.rs`).
The reference build uses the pristine source plus `patches/pilercr-trace.patch`.
Output is unchanged; the trace hooks are only needed for trace tooling.

## Performance

`dp::align_traps` parallelizes the independent per-trapezoid DP work with
Rayon. Each candidate computes its own hits and coverage marks; a sequential
pass then applies the marks in sorted trapezoid order and appends hits in the
original order. This preserves the single-threaded coverage and ordering
semantics exactly.

Measured on a shared 16-thread Linux x86-64 host (background load was present)
with `cargo build --release`, using median wall time of 15 alternating runs per
binary (process startup and output to `/dev/null` included):

| Assembly | Before (ms) | After (ms) | Speedup |
| --- | ---: | ---: | ---: |
| 9 | 624 | 374 | 1.67x |
| 16 | 491 | 331 | 1.48x |
| 17 | 554 | 345 | 1.61x |
| 20 | 655 | 409 | 1.60x |
| 50 | 918 | 454 | 2.02x |

With `RAYON_NUM_THREADS=1`, median times were within 0.3–2.5% of the previous
single-threaded binary (the DP refactor still allocates per-trapezoid
candidates). The multi-threaded improvement is material; the single-threaded
difference is not.

## Known differences from the reference

- The default output targets the modern rebuild. The original 2011 binary's
  output differs on assemblies 9, 16 and 17 because of its GCC 4.4.3 sort
  tie order; `--legacy` restores that order.
- `--legacy` is not full 2011 parity. On this code base, 51/52 locally available
  2011 text reports are byte-identical (with `--legacy-coords`, as those
  goldens use the reference's gap-shifted coordinates). The remaining case is
  **20**: the local 2011 report has one 52-copy array, while legacy mode produces
  two adjacent arrays with 24 and 28 copies. This is pinned by an ignored
  full-data test. It is not caused by the sorter (its equal-key permutations
  are pinned by unit tests); its precise cause is unresolved. Floating-point
  differences are one possibility, but have not been demonstrated.
- With default settings, corrected coordinates differ from the reference's
  text output for arrays containing gaps. `--legacy-coords` restores the
  reference's coordinates; it is independent of `--legacy`.
- Reference crashes on very short inputs (assertion failures) are not
  reproduced. The campaign records them separately as safety differences.

## Acknowledgements

- **Robert C. Edgar** wrote PILER-CR and its algorithm. This project is an
  independent reimplementation based on his public-domain source.
- **Jim Rybarski** ([pilercr-parser](https://github.com/jimrybarski/pilercr-parser),
  jim@rybarski.com, MIT) identified the gap-shift coordinate bug in PILER-CR
  v1.06 and wrote the parser used to verify coordinates here.

## License

MIT. See [LICENSE](LICENSE). PILER-CR itself is public domain.
