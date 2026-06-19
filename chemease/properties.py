import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


def calculate_properties(mol: Chem.Mol) -> dict:
    if mol is None:
        return {}

    return {
        "MW": Descriptors.ExactMolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "NumHDonors": rdMolDescriptors.CalcNumHBD(mol),
        "NumHAcceptors": rdMolDescriptors.CalcNumHBA(mol),
        "NumRotatableBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
    }


def calculate_properties_df(
    mols: list[Chem.Mol], names: list[str] = None
) -> pd.DataFrame:
    data = []
    for i, mol in enumerate(mols):
        props = calculate_properties(mol)
        if not props:
            continue
        name = None
        if names and i < len(names):
            name = names[i]
        elif mol.HasProp("_Name"):
            name = mol.GetProp("_Name")
        else:
            name = f"Mol_{i}"

        props["Name"] = name
        data.append(props)

    df = pd.DataFrame(data)
    # Reorder columns to put Name first
    if not df.empty:
        cols = ["Name"] + [c for c in df.columns if c != "Name"]
        df = df[cols]
    return df
