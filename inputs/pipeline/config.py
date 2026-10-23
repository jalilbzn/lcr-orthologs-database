import os
from pathlib import Path


def _env_int(name: str, default: int) -> int:
    """Read integer env var with a safe fallback."""
    value = os.getenv(name)
    if value is None:
        return default
    try:
        return int(value)
    except ValueError:
        return default

# ── Database ──
DB_CONFIG = {
    "host": os.getenv("LCR_DB_HOST", os.getenv("PGHOST", "127.0.0.1")),
    "port": _env_int("LCR_DB_PORT", _env_int("PGPORT", 3000)),
    "dbname": os.getenv("LCR_DB_NAME", os.getenv("PGDATABASE", "lcr_matrix")),
    "user": os.getenv("LCR_DB_USER", os.getenv("PGUSER", "postgres")),
    "password": os.getenv("LCR_DB_PASSWORD", os.getenv("PGPASSWORD", "")),
}

# ── Project paths ──
PROJECT_ROOT = Path(__file__).resolve().parent.parent
ARTIFACTS_DIR = PROJECT_ROOT / "artifacts"
REPORTS_DIR = ARTIFACTS_DIR / "reports"
LOGS_DIR = ARTIFACTS_DIR / "logs"
VISUALS_DIR = ARTIFACTS_DIR / "visuals"
HISTORY_DIR = ARTIFACTS_DIR / "history"
DATA_DIR = PROJECT_ROOT / "data"
OMA_DIR = PROJECT_ROOT / "oma"

for _dir in (ARTIFACTS_DIR, REPORTS_DIR, LOGS_DIR, VISUALS_DIR, HISTORY_DIR, DATA_DIR, OMA_DIR):
    _dir.mkdir(parents=True, exist_ok=True)

# ── LCRAnnotationsDB API ──
API_BASE = "https://lcrannotdb.lcr-lab.org/api"
API_ENDPOINTS = {
    "proteins": f"{API_BASE}/proteins/",
    "lcrs":     f"{API_BASE}/lcrs/",
    "annotations": f"{API_BASE}/annotations/",
    "categories":  f"{API_BASE}/categories/",
}

# ── Rate limiting ──
REQUEST_DELAY_SEC = 1.5  # seconds between API requests
