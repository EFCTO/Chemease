from rdkit import Chem

from chemease import calculate_similarity, find_similar


def test_calculate_similarity():
    mol1 = Chem.MolFromSmiles("CCO")
    mol2 = Chem.MolFromSmiles("CCO")
    mol3 = Chem.MolFromSmiles("CC")

    sim_identical = calculate_similarity(mol1, mol2)
    assert sim_identical == 1.0

    sim_different = calculate_similarity(mol1, mol3)
    assert sim_different < 1.0


def test_find_similar():
    query = Chem.MolFromSmiles("CCO")
    targets = [
        Chem.MolFromSmiles("CCO"),
        Chem.MolFromSmiles("CCCO"),
        Chem.MolFromSmiles("c1ccccc1"),
    ]

    results = find_similar(query, targets, threshold=0.1)
    assert len(results) > 0
    assert results[0][1] == 1.0

# hell