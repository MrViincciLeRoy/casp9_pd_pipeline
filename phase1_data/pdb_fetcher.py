import requests
from pathlib import Path
from utils import get_logger
from config import STRUCTURES_DIR

log = get_logger(__name__)

def fetch_pdb(pdb_id: str) -> Path:
    out = STRUCTURES_DIR / f"{pdb_id}.pdb"
    if out.exists():
        log.info(f"{pdb_id}.pdb already downloaded")
        return out
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    log.info(f"Downloading {pdb_id} from RCSB...")
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    out.write_text(r.text)
    log.info(f"Saved → {out}")
    return out

def run():
    from config import CASP9_PDB_ID, ALPHASYN_PDB_ID
    fetch_pdb(CASP9_PDB_ID)
    fetch_pdb(ALPHASYN_PDB_ID)
