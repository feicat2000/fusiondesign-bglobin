#!/usr/bin/env bash
# =============================================================================
# mRNA Neoantigen Vaccine Pipeline — Environment Setup
# =============================================================================
# Creates a Python virtual environment, installs all dependencies,
# and downloads required model data (mhcflurry weights, codon tables).
#
# Usage:
#   chmod +x setup.sh
#   ./setup.sh
#
# After setup:
#   source venv/bin/activate
#   python pipeline.py --help
# =============================================================================

set -euo pipefail

VENV_DIR="venv"
PYTHON_MIN="3.9"

# ── helpers ──────────────────────────────────────────────────────────────────
log()  { echo -e "\033[1;34m[SETUP]\033[0m $*"; }
ok()   { echo -e "\033[1;32m[  OK ]\033[0m $*"; }
warn() { echo -e "\033[1;33m[ WARN]\033[0m $*"; }
fail() { echo -e "\033[1;31m[FAIL ]\033[0m $*"; exit 1; }

# ── check Python version ─────────────────────────────────────────────────────
log "Checking Python version …"
PYTHON=$(command -v python3 || fail "python3 not found. Install Python >= ${PYTHON_MIN}")

PY_VER=$("$PYTHON" -c "import sys; print('%d.%d' % sys.version_info[:2])")
PY_MAJOR=$("$PYTHON" -c "import sys; print(sys.version_info[0])")
PY_MINOR=$("$PYTHON" -c "import sys; print(sys.version_info[1])")

if [[ "$PY_MAJOR" -lt 3 || ( "$PY_MAJOR" -eq 3 && "$PY_MINOR" -lt 9 ) ]]; then
    fail "Python >= 3.9 required, found $PY_VER"
fi
ok "Python $PY_VER"

# ── create virtual environment ────────────────────────────────────────────────
if [[ -d "$VENV_DIR" ]]; then
    warn "Virtual environment '$VENV_DIR' already exists — skipping creation."
    warn "To rebuild from scratch: rm -rf $VENV_DIR && ./setup.sh"
else
    log "Creating virtual environment in ./$VENV_DIR …"
    "$PYTHON" -m venv "$VENV_DIR"
    ok "Virtual environment created."
fi

# activate
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"
ok "Virtual environment activated."

# ── upgrade pip ───────────────────────────────────────────────────────────────
log "Upgrading pip …"
pip install --quiet --upgrade pip wheel
# Install setuptools < 74: mhcflurry depends on pkg_resources which was
# separated from setuptools in v74+. Pinning below 74 keeps it bundled.
pip install --quiet "setuptools<74"
ok "pip upgraded."

# ── install ViennaRNA (system package preferred, Python binding via pip) ───────
log "Checking ViennaRNA Python bindings …"
if python -c "import RNA" 2>/dev/null; then
    ok "ViennaRNA Python bindings already available."
else
    warn "ViennaRNA 'RNA' module not found in venv."
    # Try the unofficial PyPI wheel first (works on most Linux/Mac)
    if pip install --quiet ViennaRNA 2>/dev/null; then
        ok "ViennaRNA installed via PyPI."
    else
        warn "PyPI ViennaRNA failed. Falling back to system package instructions:"
        warn "  Ubuntu/Debian : sudo apt install vienna-rna python3-rna"
        warn "  Conda         : conda install -c bioconda viennarna"
        warn "  Source        : https://www.tbi.univie.ac.at/RNA/#download"
        warn "Continuing — stability analysis will use fallback GC/MFE estimates."
    fi
fi

# ── install Python dependencies ───────────────────────────────────────────────
log "Installing Python dependencies from requirements.txt …"
pip install --quiet -r requirements.txt
ok "Dependencies installed."

# ── download mhcflurry models ─────────────────────────────────────────────────
log "Downloading mhcflurry pan-allele MHC-I models (~500 MB, cached) …"
if mhcflurry-downloads info 2>/dev/null | grep -q "available"; then
    ok "mhcflurry models already downloaded."
else
    mhcflurry-downloads fetch || warn "mhcflurry model download failed — check internet connection."
fi

# ── create output directories ─────────────────────────────────────────────────
log "Creating output directories …"
mkdir -p output/{constructs,mrna,analysis,reports}
ok "Output directories ready."

# ── pepsickle (optional — proteasomal cleavage) ───────────────────────────────
log "Checking pepsickle (proteasomal cleavage predictor) …"
if python -c "import pepsickle" 2>/dev/null; then
    ok "pepsickle available."
else
    warn "pepsickle not installed. Attempting …"
    pip install --quiet pepsickle 2>/dev/null \
        && ok "pepsickle installed." \
        || warn "pepsickle unavailable — will use NetChop web API fallback."
fi

# ── netMHCpan (optional external binary) ─────────────────────────────────────
if command -v netMHCpan &>/dev/null; then
    ok "netMHCpan binary found in PATH."
else
    warn "netMHCpan not in PATH (optional). Download from:"
    warn "  https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/"
    warn "  Set NETMHCPAN_PATH env var to its location if installed locally."
fi

# ── verify core imports ───────────────────────────────────────────────────────
log "Verifying core imports …"
python - <<'PYCHECK'
import sys, importlib
required = [
    ("Bio",             "biopython"),
    ("mhcflurry",       "mhcflurry"),
    ("pandas",          "pandas"),
    ("numpy",           "numpy"),
    ("matplotlib",      "matplotlib"),
    ("python_codon_tables", "python-codon-tables"),
    ("dnachisel",       "dnachisel"),
    ("requests",        "requests"),
    ("jinja2",          "jinja2"),
]
failed = []
for mod, pkg in required:
    try:
        importlib.import_module(mod)
        print(f"  \033[32m✓\033[0m {pkg}")
    except ImportError:
        print(f"  \033[31m✗\033[0m {pkg}  ← MISSING")
        failed.append(pkg)

optional = [
    ("RNA",             "ViennaRNA"),
    ("pepsickle",       "pepsickle"),
]
for mod, pkg in optional:
    try:
        importlib.import_module(mod)
        print(f"  \033[32m✓\033[0m {pkg} (optional)")
    except ImportError:
        print(f"  \033[33m~\033[0m {pkg} (optional — fallback active)")

if failed:
    print(f"\n\033[31mSetup incomplete — missing: {failed}\033[0m")
    sys.exit(1)
PYCHECK

ok "All core imports verified."

# ── done ──────────────────────────────────────────────────────────────────────
echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║  Setup complete!                                             ║"
echo "╠══════════════════════════════════════════════════════════════╣"
echo "║  Activate environment : source venv/bin/activate             ║"
echo "║  Run example          : python pipeline.py \\                 ║"
echo "║       --fasta data/example_neoantigens.fasta \\               ║"
echo "║       --hla HLA-A*02:01,HLA-B*07:02                         ║"
echo "║  Help                 : python pipeline.py --help            ║"
echo "╚══════════════════════════════════════════════════════════════╝"
