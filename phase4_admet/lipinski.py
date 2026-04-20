from utils import get_logger

log = get_logger(__name__)

def lipinski_filter(smiles: str) -> bool:
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        from config import LIPINSKI_MW_MAX, LIPINSKI_LOGP_MAX, LIPINSKI_HBD_MAX, LIPINSKI_HBA_MAX
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        return (
            Descriptors.MolWt(mol) <= LIPINSKI_MW_MAX
            and Descriptors.MolLogP(mol) <= LIPINSKI_LOGP_MAX
            and rdMolDescriptors.CalcNumHBD(mol) <= LIPINSKI_HBD_MAX
            and rdMolDescriptors.CalcNumHBA(mol) <= LIPINSKI_HBA_MAX
        )
    except Exception as e:
        log.warning(f"Lipinski check error: {e}")
        return False
