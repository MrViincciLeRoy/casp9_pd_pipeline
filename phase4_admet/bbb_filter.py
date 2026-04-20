import requests
from utils import get_logger
from config import PKCSM_API_URL

log = get_logger(__name__)

def bbb_check(smiles: str) -> bool:
    try:
        r = requests.post(PKCSM_API_URL, json={"smiles": smiles}, timeout=15)
        r.raise_for_status()
        return r.json().get("BBB", False)
    except Exception as e:
        log.warning(f"BBB check failed for {smiles[:20]}...: {e}")
        return False

def run_admet_filter(candidates: list[dict]) -> list[dict]:
    from .lipinski import lipinski_filter
    import json
    from config import RESULTS_DIR

    log.info(f"ADMET filtering {len(candidates)} candidates...")
    passed = []
    for c in candidates:
        smiles = c.get("smiles", "")
        if not lipinski_filter(smiles):
            continue
        if not bbb_check(smiles):
            continue
        passed.append({**c, "lipinski": True, "bbb_permeable": True})

    out = RESULTS_DIR / "admet_passed.json"
    out.write_text(json.dumps(passed, indent=2))
    log.info(f"{len(passed)} candidates passed ADMET → {out}")
    return passed
