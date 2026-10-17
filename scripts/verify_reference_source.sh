#!/usr/bin/env bash
set -euo pipefail
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
expected=$(awk '{print $2}' ref/SHA256SUMS | sort)
actual=$(find ref/pilercr1.06-src -type f | sort)
if [ "$expected" != "$actual" ]; then
  echo "file list differs from ref/SHA256SUMS" >&2
  diff <(printf '%s\n' "$expected") <(printf '%s\n' "$actual") >&2 || true
  exit 1
fi
sha256sum -c ref/SHA256SUMS
