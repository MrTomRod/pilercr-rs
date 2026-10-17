# Reference source: PILER-CR 1.06

`pilercr1.06-src/` is the unmodified upstream PILER-CR 1.06 C++ source.

- Upstream: https://www.drive5.com/pilercr/pilercr1.06.tar.gz
- Tarball SHA256: `50175f7aa171674cda5ba255631f340f9cc7f80e8cc25135a4cb857147d91068`
- Extracted from the tarball's `pilercr1.06/` directory, flattened. The tarball's prebuilt
  `pilercr` and `pilercr.exe` binaries are not included.
- `pilercr1.06-src/LICENSE.md` is not part of the tarball. It was copied from the project's
  earlier local copy and states that PILER, including source and documentation, is donated to
  the public domain. The drive5.com PILER-CR page also describes PILERCR as public domain software.
- Snapshot date: 2026-10-17.

Integrity: `SHA256SUMS` lists every file. Check with `bash scripts/verify_reference_source.sh`
from the repository root.

Trace hooks are not part of this tree. They live in `patches/pilercr-trace.patch`, which
`scripts/build_reference_trace.sh` applies to a copy under `target/reference-trace/`.
