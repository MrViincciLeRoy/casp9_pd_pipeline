import requests
import json
from pathlib import Path
from utils import get_logger
from config import ZINC_DIR

log = get_logger(__name__)

ZINC_API = "https://zinc.docking.org/substances/subsets/lead-like.json"

def fetch_zinc_subset(limit: int = 10000) -> list[str]:
    out = ZINC_DIR / f"zinc_subset_{limit}.json"
    if out.exists():
        log.info("ZINC subset already cached")
        return json.loads(out.read_text())

    log.info(f"Fetching {limit} ZINC lead-like compounds...")
    smiles_list = []
    page = 1
    while len(smiles_list) < limit:
        r = requests.get(ZINC_API, params={"page": page, "count": 100}, timeout=30)
        if r.status_code != 200:
            break
        batch = r.json()
        if not batch:
            break
        smiles_list += [c["smiles"] for c in batch if c.get("smiles")]
        page += 1

    smiles_list = smiles_list[:limit]
    out.write_text(json.dumps(smiles_list, indent=2))
    log.info(f"Saved {len(smiles_list)} ZINC SMILES → {out}")
    return smiles_list

def run():
    fetch_zinc_subset()
