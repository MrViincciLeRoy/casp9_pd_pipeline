import subprocess
import json
from pathlib import Path
from utils import get_logger
from config import STRUCTURES_DIR, RESULTS_DIR, DOCKING_CENTER, DOCKING_BOX_SIZE, CASP9_PDB_ID

log = get_logger(__name__)

def dock_with_gnina(ligand_pdbqt: Path, receptor_pdbqt: Path | None = None) -> dict | None:
    if receptor_pdbqt is None:
        receptor_pdbqt = STRUCTURES_DIR / f"{CASP9_PDB_ID}_prepared.pdbqt"
    out_path = ligand_pdbqt.with_suffix(".gnina.pdbqt")
    cx, cy, cz = DOCKING_CENTER
    sx, sy, sz = DOCKING_BOX_SIZE
    cmd = [
        "gnina",
        "-r", str(receptor_pdbqt),
        "-l", str(ligand_pdbqt),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "-o", str(out_path),
        "--cnn_scoring", "rescore",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        score_line = [l for l in result.stdout.splitlines() if "Affinity" in l]
        score = float(score_line[0].split()[1]) if score_line else None
        return {"file": str(ligand_pdbqt), "gnina_score": score}
    except Exception as e:
        log.warning(f"GNINA failed for {ligand_pdbqt.name}: {e}")
        return None

def run_gnina_validation(top_hits: list[dict]) -> list[dict]:
    log.info(f"GNINA validating top {len(top_hits)} Vina hits...")
    validated = []
    for hit in top_hits:
        result = dock_with_gnina(Path(hit["file"]))
        if result:
            validated.append({**hit, **result})
    validated.sort(key=lambda x: x.get("gnina_score") or 0)
    out = RESULTS_DIR / "gnina_validated.json"
    out.write_text(json.dumps(validated, indent=2))
    log.info(f"GNINA validation done → {out}")
    return validated

