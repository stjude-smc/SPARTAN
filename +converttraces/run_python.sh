#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONVERTER="$SCRIPT_DIR/tracesMat2Npz.py"

usage() {
  echo "Usage: $0 <python_executable> <input.mat> <output.npz> [--var NAME]" >&2
  exit 1
}

[[ $# -ge 3 ]] || usage

PY="$1"
[[ -x "$PY" ]] || { echo "error: not executable: $PY" >&2; exit 1; }

exec "$PY" "$CONVERTER" "$2" "$3" "${@:4}"