from rdkit import Chem

from chemease import calculate_properties, calculate_properties_df


def test_calculate_properties():
    mol = Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin
    props = calculate_properties(mol)
    assert "MW" in props
    assert "LogP" in props
    assert "TPSA" in props
    assert props["MW"] > 100


def test_calculate_properties_df():
    mols = [Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O"), Chem.MolFromSmiles("CCO")]
    df = calculate_properties_df(mols, names=["Aspirin", "Ethanol"])
    assert len(df) == 2
    assert "Name" in df.columns
    assert df.iloc[0]["Name"] == "Aspirin"
    assert df.iloc[1]["Name"] == "Ethanol"
    assert "MW" in df.columns
