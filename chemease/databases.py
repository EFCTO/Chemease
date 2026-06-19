from typing import Dict, List, Optional

from chembl_webresource_client.new_client import new_client
from rdkit import Chem

from chemease.exceptions import FetchError, InvalidMoleculeError
from chemease.logger import logger


def get_target_by_name(target_name: str) -> List[Dict]:
    logger.info(f"Searching ChEMBL target: {target_name}")
    target = new_client.target
    res = target.filter(pref_name__icontains=target_name).only(
        ["target_chembl_id", "pref_name", "organism"]
    )
    return list(res)


def get_bioactivities(
    target_chembl_id: str, activity_type: str = "IC50", limit: int = 100
) -> List[Dict]:
    logger.info(
        f"Fetching bioactivities for target {target_chembl_id} (type: {activity_type})"
    )
    activity = new_client.activity
    res = activity.filter(target_chembl_id=target_chembl_id).filter(
        standard_type=activity_type
    )
    return list(res[:limit])


def get_mol_from_chembl(chembl_id: str) -> Optional[Chem.Mol]:
    logger.info(f"Fetching molecule from ChEMBL: {chembl_id}")
    molecule = new_client.molecule
    try:
        res = molecule.filter(molecule_chembl_id=chembl_id).only(
            ["molecule_structures"]
        )
    except Exception as e:
        raise FetchError(f"Failed to fetch {chembl_id} from ChEMBL: {e}")

    if res and res[0].get("molecule_structures"):
        smiles = res[0]["molecule_structures"].get("canonical_smiles")
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol.SetProp("_Name", chembl_id)
                mol.SetProp("ChEMBL_ID", chembl_id)
                return mol
            else:
                logger.warning(f"Failed to parse SMILES for {chembl_id}")
                raise InvalidMoleculeError(f"Invalid SMILES for {chembl_id}: {smiles}")

    logger.warning(f"No structure found for ChEMBL ID: {chembl_id}")
    return None
