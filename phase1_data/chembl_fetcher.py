import requests
import json
from pathlib import Path
from utils import get_logger
from config import CHEMBL_DIR, CASP9_CHEMBL_ID, CHEMBL_LIMIT

log = get_logger(__name__)

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"

def fetch_casp9_inhibitors() -> list[dict]:
    out = CHEMBL_DIR / "casp9_inhibitors.json"
    if out.exists():
        log.info("ChEMBL data already cached")
        return json.loads(out.read_text())

    log.info("Fetching Caspase-9 inhibitors from ChEMBL...")
    params = {
        "target_chembl_id": CASP9_CHEMBL_ID,
        "standard_type": "Ki",
        "limit": CHEMBL_LIMIT,
        "format": "json",
    }
    r = requests.get(f"{CHEMBL_BASE}/activity", params=params, timeout=30)
    r.raise_for_status()
    activities = r.json().get("activities", [])

    cleaned = [
        {
            "chembl_id": a.get("molecule_chembl_id"),
            "smiles": a.get("canonical_smiles"),
            "ki_nm": float(a["standard_value"]) if a.get("standard_value") else None,
        }
        for a in activities
        if a.get("canonical_smiles") and a.get("standard_value")
    ]
    out.write_text(json.dumps(cleaned, indent=2))
    log.info(f"Saved {len(cleaned)} inhibitors → {out}")
    return cleaned

def run():
    fetch_casp9_inhibitors()


# allow env override
import os as _os
CHEMBL_LIMIT = int(_os.getenv("CHEMBL_LIMIT", str(CHEMBL_LIMIT)))