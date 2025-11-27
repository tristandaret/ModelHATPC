#!/usr/bin/env bash
set -euo pipefail
python - <<'PY'
import os,sys
print('ModelHATPC container ready. Python', sys.version)
print('Top-level files:', sorted(os.listdir('.')))
PY
