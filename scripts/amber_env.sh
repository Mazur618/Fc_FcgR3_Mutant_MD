#!/usr/bin/env bash
# Optional setup script; otherwise use tools already available on PATH.
if [ -n "${AMBER_ENV:-}" ]; then
  [ -r "$AMBER_ENV" ] || { echo "Cannot read AMBER_ENV: $AMBER_ENV" >&2; exit 1; }
  source "$AMBER_ENV"
fi
command -v cpptraj >/dev/null || { echo "Load AmberTools or set AMBER_ENV to its amber.sh" >&2; exit 1; }
