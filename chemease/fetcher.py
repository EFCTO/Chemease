"""PubChem data fetching module."""

from typing import Optional

import pubchempy as pcp
from rdkit import Chem

from chemease.exceptions import FetchError, InvalidMoleculeError
from chemease.logger import logger


def get_mol_from_name(name: str) -> Chem.Mol:
    logger.info(f"Fetching compound by name: {name}")
    compounds = pcp.get_compounds(name, "name")
    if not compounds:
        logger.error(f"No compound found for name: {name}")
        raise FetchError(f"No compound found for name: {name}")

    c = compounds[0]
    smiles = getattr(c, "isomeric_smiles", getattr(c, "canonical_smiles", None))
    if not smiles:
        raise FetchError(f"No SMILES found for name: {name}")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidMoleculeError(f"Failed to parse SMILES for {name}: {smiles}")

    mol.SetProp("_Name", name)
    mol.SetProp("CID", str(c.cid))
    return mol


def get_mol_from_cid(cid: int) -> Chem.Mol:
    logger.info(f"Fetching compound by CID: {cid}")
    try:
        c = pcp.Compound.from_cid(cid)
    except Exception as e:
        raise FetchError(f"Failed to fetch CID {cid}: {e}")

    smiles = getattr(c, "isomeric_smiles", getattr(c, "canonical_smiles", None))
    if not smiles:
        raise FetchError(f"No SMILES found for CID: {cid}")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidMoleculeError(f"Failed to parse SMILES for CID {cid}: {smiles}")

    if c.iupac_name:
        mol.SetProp("_Name", c.iupac_name)
    mol.SetProp("CID", str(c.cid))
    return mol
