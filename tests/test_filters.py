from rdkit import Chem

from chemease import filter_by_substructure, passes_lipinski


def test_passes_lipinski():
    mol_pass = Chem.MolFromSmiles("CCO")  # Ethanol
    assert passes_lipinski(mol_pass)

    mol_fail = Chem.MolFromSmiles("C" * 100)  # MW > 500, LogP > 5
    assert not passes_lipinski(mol_fail)


def test_filter_by_substructure():
    mols = [
        Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O"),  # Aspirin
        Chem.MolFromSmiles("CCO"),  # Ethanol
        Chem.MolFromSmiles("c1ccccc1"),  # Benzene
    ]
    filtered = filter_by_substructure(mols, "c1ccccc1")
    assert len(filtered) == 2
