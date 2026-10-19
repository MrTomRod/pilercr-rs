#!/usr/bin/env bash
# Fetch the original 2011 PILER-CR 1.06 binary (not committed) into target/reference-binary/.
# Diagnostic only: its reports differ from a modern rebuild on assemblies 9, 16 and 17.
set -euo pipefail

root=$(CDPATH='' cd -- "$(dirname -- "$0")/.." && pwd)
url=https://www.drive5.com/pilercr/pilercr1.06.tar.gz
sha=50175f7aa171674cda5ba255631f340f9cc7f80e8cc25135a4cb857147d91068
dir="$root/target/reference-binary"

mkdir -p "$dir"
curl -fsSL -o "$dir/pilercr1.06.tar.gz" "$url"
[ "$(sha256sum "$dir/pilercr1.06.tar.gz" | cut -d' ' -f1)" = "$sha" ] ||
  { echo "checksum mismatch for $url" >&2; exit 1; }
tar -xzf "$dir/pilercr1.06.tar.gz" -C "$dir" --strip-components=2 ./pilercr1.06/pilercr
chmod +x "$dir/pilercr"

echo "$dir/pilercr"
