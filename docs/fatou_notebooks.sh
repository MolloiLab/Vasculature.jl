#!/usr/bin/env bash
# Format (or check) and lint the Pluto notebooks in docs/notebooks/ with Fatou.
#
#   docs/fatou_notebooks.sh          format the notebooks' cells in place, then lint them
#   docs/fatou_notebooks.sh --check  fail if a cell is unformatted or a lint rule fires (CI)
#
# fatou.toml excludes docs/notebooks/ because Pluto rewrites its own header (the `@bind` mock
# above the first cell) on every save, and Fatou has no directive to leave a region alone.
# This script hands Fatou everything from the first cell marker on, so the cells follow
# Fatou's defaults and the header stays exactly as Pluto writes it.
set -euo pipefail

FATOU="${FATOU:-fatou}"
check=false
[[ "${1:-}" == "--check" ]] && check=true

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

notebooks=("$root"/docs/notebooks/*.jl)
for nb in "${notebooks[@]}"; do
    first=$(grep -n -m1 '^# ╔═╡ ' "$nb" | cut -d: -f1)
    head -n $((first - 1)) "$nb" >"$tmp/$(basename "$nb").header"
    tail -n +"$first" "$nb" >"$tmp/$(basename "$nb")"
done

# --no-config: the repository's fatou.toml would exclude these copies' names; the defaults apply
if $check; then
    (cd "$tmp" && "$FATOU" format --no-config --check ./*.jl)
else
    (cd "$tmp" && "$FATOU" format --no-config ./*.jl)
    for nb in "${notebooks[@]}"; do
        cat "$tmp/$(basename "$nb").header" "$tmp/$(basename "$nb")" >"$nb"
    done
fi
(cd "$tmp" && "$FATOU" lint --no-config ./*.jl)
