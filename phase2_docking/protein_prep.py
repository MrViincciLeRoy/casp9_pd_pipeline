import subprocess
from pathlib import Path
from utils import get_logger
from config import STRUCTURES_DIR, CASP9_PDB_ID

log = get_logger(__name__)

def prepare_receptor(pdb_id: str = CASP9_PDB_ID) -> Path:
    pdb = STRUCTURES_DIR / f"{pdb_id}.pdb"
    pdbqt = STRUCTURES_DIR / f"{pdb_id}_prepared.pdbqt"
    if pdbqt.exists():
        log.info(f"{pdbqt.name} already prepared")
        return pdbqt
    # requires: conda install -c conda-forge autodock-utilities
    log.info(f"Preparing receptor {pdb_id} → pdbqt...")
    cmd = ["prepare_receptor", "-r", str(pdb), "-o", str(pdbqt), "-A", "hydrogens"]
    subprocess.run(cmd, check=True)
    log.info(f"Saved → {pdbqt}")
    return pdbqt

def prepare_ligand(smiles: str, out_path: Path) -> Path | None:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)
        Chem.MolToMolFile(mol, str(out_path.with_suffix(".mol")))
        subprocess.run(["mk_prepare_ligand.py", "-i", str(out_path.with_suffix(".mol")), "-o", str(out_path)], check=True)
        return out_path
    except Exception as e:
        log.warning(f"Ligand prep failed: {e}")
        return None
