#!/usr/bin/env bash
set -euo pipefail
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
work="$root/target/reference-trace/src"
rm -rf "$work"
mkdir -p "$work"
cp -p "$root"/ref/pilercr1.06-src/*.cpp "$root"/ref/pilercr1.06-src/*.h "$work"/
(cd "$work" && patch -s -p1 < "$root/patches/pilercr-trace.patch")
"${CXX:-c++}" -std=c++17 -O3 -DNDEBUG=1 \
  "$work"/*.cpp -lm \
  -o "$root/target/reference-trace/pilercr"
