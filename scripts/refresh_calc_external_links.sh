#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HOST="${LO_UNO_HOST:-127.0.0.1}"
PORT="${LO_UNO_PORT:-2002}"
PROFILE_DIR="${LO_UNO_PROFILE_DIR:-/tmp/lo_uno_profile}"
TMP_DIR="${LO_UNO_TMPDIR:-/tmp/lo_uno_tmp}"
CSV_GLOB="${1:-package_tests/comprehensive_tests_results_nc_1_*.csv}"

listener_pattern="soffice.bin.*socket,host=${HOST},port=${PORT};urp"

ensure_listener() {
  if pgrep -af "${listener_pattern}" >/dev/null 2>&1; then
    return
  fi

  mkdir -p "${PROFILE_DIR}" "${TMP_DIR}"
  (
    cd "${ROOT_DIR}"
    export HOME=/tmp
    export TMPDIR="${TMP_DIR}"
    exec libreoffice \
      --headless \
      --nologo \
      --nodefault \
      --norestore \
      "-env:UserInstallation=file://${PROFILE_DIR}" \
      "--accept=socket,host=${HOST},port=${PORT};urp" \
      >/tmp/lo_uno_listener.log 2>&1
  ) &

  for _ in $(seq 1 30); do
    sleep 1
    if pgrep -af "${listener_pattern}" >/dev/null 2>&1; then
      return
    fi
  done

  echo "LibreOffice UNO listener failed to start." >&2
  echo "See /tmp/lo_uno_listener.log for details." >&2
  exit 1
}

cd "${ROOT_DIR}"
ensure_listener
python3 scripts/build_calc_external_links_uno.py \
  --host "${HOST}" \
  --port "${PORT}" \
  "${CSV_GLOB}"
