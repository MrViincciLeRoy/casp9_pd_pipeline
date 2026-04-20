import json
from pathlib import Path
from typing import Optional
from utils import get_logger
from config import STRUCTURES_DIR, RESULTS_DIR, DOCKING_CENTER, DOCKING_BOX_SIZE, DOCKING_EXHAUSTIVENESS, DOCKING_N_POSES, DOCKING_SCORE_CUTOFF, CASP9_PDB_ID

log = get_logger(__name__)

def dock_compound(ligand_pdbqt: Path, receptor_pdbqt: Optional[Path] = None) -> Optional[float]:
    try:
        from vina import Vina
        if receptor_pdbqt is None:
            receptor_pdbqt = STRUCTURES_DIR / f"{CASP9_PDB_ID}_prepared.pdbqt"
        v = Vina(sf_name="vina", verbosity=0)
        v.set_receptor(str(receptor_pdbqt))
        v.compute_vina_maps(center=DOCKING_CENTER, box_size=DOCKING_BOX_SIZE)
        v.set_ligand_from_file(str(ligand_pdbqt))
        v.dock(exhaustiveness=DOCKING_EXHAUSTIVENESS, n_poses=DOCKING_N_POSES)
        return v.score()[0]
    except Exception as e:
        log.warning(f"Docking failed for {ligand_pdbqt.name}: {e}")
        return None

def run_vina_screen(compound_pdbqts: list[Path]) -> list[dict]:
    log.info(f"Vina screening {len(compound_pdbqts)} compounds...")
    results = []
    for pdbqt in compound_pdbqts:
        score = dock_compound(pdbqt)
        if score and score < DOCKING_SCORE_CUTOFF:
            results.append({"file": str(pdbqt), "vina_score": score})
    results.sort(key=lambda x: x["vina_score"])
    out = RESULTS_DIR / "vina_hits.json"
    out.write_text(json.dumps(results, indent=2))
    log.info(f"{len(results)} hits saved → {out}")
    return results
