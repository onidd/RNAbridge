import os

from pathlib import Path

# Base project directory
BASE_DIR = Path(__file__).parent.absolute()

# Data Directories
DATASET = os.environ.get("RNABRIDGE_DATASET", "pdb").lower()
DATA_DIR = BASE_DIR / "rnasolo" if DATASET == "rnasolo" else BASE_DIR

CIF_DIR = Path(os.environ.get("RNABRIDGE_CIF_DIR", DATA_DIR / "cif"))
JSON_DIR = Path(os.environ.get("RNABRIDGE_JSON_DIR", DATA_DIR / "json"))
RESULTS_DIR = Path(os.environ.get("RNABRIDGE_RESULTS_DIR", DATA_DIR / "analyze"))
INPUT_DIR = RESULTS_DIR
OUTPUT_DIR = Path(os.environ.get("RNABRIDGE_OUTPUT_DIR", DATA_DIR / "result"))

# Ensure directories exist
for d in [CIF_DIR, JSON_DIR, RESULTS_DIR, OUTPUT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Analysis Parameters
MAX_BEND_ANGLE = float(os.environ.get("RNABRIDGE_MAX_ANGLE", 50.0))

# Database Configuration
DATABASE_URL = os.environ.get("DATABASE_URL", "sqlite:///rnabridge.db")

# External Tools Paths
VARNA_CLI_PATH = os.environ.get("VARNA_CLI_PATH", "cli2rest-bio")

# X3DNA Configuration
# Automatically setup X3DNA environment if directory exists in project root
X3DNA_DIR = Path(os.environ.get("X3DNA", BASE_DIR / "x3dna-v2.4"))
if X3DNA_DIR.exists():
    os.environ["X3DNA"] = str(X3DNA_DIR)
    x3dna_bin = X3DNA_DIR / "bin"
    if str(x3dna_bin) not in os.environ.get("PATH", ""):
        os.environ["PATH"] = str(x3dna_bin) + os.pathsep + os.environ.get("PATH", "")
