import os

from rdkit import Chem


def read_sdf(filepath: str, sanitize: bool = True) -> list[Chem.Mol]:
    """read molecules from an SDF file."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    suppl = Chem.SDMolSupplier(filepath, sanitize=sanitize)
    return [mol for mol in suppl if mol is not None]


def write_sdf(mols: list[Chem.Mol], filepath: str):
    """Write molecules to an SDF file."""
    writer = Chem.SDWriter(filepath)
    for mol in mols:
        if mol is not None:
            writer.write(mol)
    writer.close()


def read_smiles(
    filepath: str, delimiter: str = " ", smilesColumn: int = 0, nameColumn: int = 1
) -> list[Chem.Mol]:
    """read molecules from a SMILES (.smi/.csv) file"""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    suppl = Chem.SmilesMolSupplier(
        filepath,
        delimiter=delimiter,
        smilesColumn=smilesColumn,
        nameColumn=nameColumn,
        titleLine=False,
    )
    return [mol for mol in suppl if mol is not None]
