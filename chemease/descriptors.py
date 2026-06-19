import pandas as pd
from rdkit import Chem
from rdkit.Chem import MACCSkeys

from chemease.properties import calculate_properties


def get_maccs_keys(mol: Chem.Mol) -> list[int]:
    if mol is None:
        return [0] * 167
    fp = MACCSkeys.GenMACCSKeys(mol)
    return list(fp)


def calculate_advanced_descriptors_df(
    mols: list[Chem.Mol], names: list[str] = None
) -> pd.DataFrame:
    data = []
    for i, mol in enumerate(mols):
        if mol is None:
            continue
        props = calculate_properties(mol)
        name = (
            names[i]
            if names and i < len(names)
            else mol.GetProp("_Name") if mol.HasProp("_Name") else f"Mol_{i}"
        )
        props["Name"] = name

        # Add MACCS keys
        maccs = get_maccs_keys(mol)
        for j in range(1, 167):
            props[f"MACCS_{j}"] = maccs[j]

        data.append(props)

    df = pd.DataFrame(data)
    if not df.empty:
        cols = ["Name"] + [c for c in df.columns if c != "Name"]
        df = df[cols]
    return df
