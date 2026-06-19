import pytest
from unittest.mock import patch, MagicMock
from rdkit import Chem

from chemease import (calculate_advanced_descriptors_df, cluster_molecules,
                      generate_conformers, get_conformer_energies,
                      get_maccs_keys, get_mol_from_chembl, get_murcko_scaffold,
                      get_target_by_name, perform_pca, run_virtual_reaction,
                      synthesize_amides)

@patch('chemease.databases.new_client')
def test_databases(mock_new_client):
    # Mock target
    mock_target_filter = MagicMock()
    mock_target_filter.only.return_value = [{'target_chembl_id': 'CHEMBL203', 'pref_name': 'EGFR', 'organism': 'Homo sapiens'}]
    mock_new_client.target.filter.return_value = mock_target_filter
    
    # Mock molecule
    mock_mol_filter = MagicMock()
    mock_mol_filter.only.return_value = [{'molecule_structures': {'canonical_smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O'}}]
    mock_new_client.molecule.filter.return_value = mock_mol_filter

    targets = get_target_by_name("EGFR")
    assert len(targets) > 0
    
    mol = get_mol_from_chembl("CHEMBL25")
    assert mol is not None
    assert mol.GetNumAtoms() > 0

def test_descriptors():
    mol = Chem.MolFromSmiles("CCO")
    maccs = get_maccs_keys(mol)
    assert len(maccs) == 167

    df = calculate_advanced_descriptors_df([mol])
    assert "MACCS_1" in df.columns
    assert len(df) == 1

def test_conformers():
    mol = Chem.MolFromSmiles("CCO")
    mol_3d = generate_conformers(mol, num_confs=3, optimize=False)
    assert mol_3d.GetNumConformers() == 3

    energies = get_conformer_energies(mol_3d)
    assert len(energies) == 3

def test_analysis():
    mols = [
        Chem.MolFromSmiles("c1ccccc1CC"),
        Chem.MolFromSmiles("c1ccccc1C"),
        Chem.MolFromSmiles("CCO"),
    ]
    scaffold = get_murcko_scaffold(mols[0])
    assert scaffold != ""

    clusters = cluster_molecules(mols, cutoff=0.5)
    assert len(clusters) > 0

    df_pca = perform_pca(mols)
    assert len(df_pca) == 3
    assert "PC1" in df_pca.columns

def test_reactions():
    acid = Chem.MolFromSmiles("CC(=O)O")
    amine = Chem.MolFromSmiles("NCC")
    amides = synthesize_amides([acid], [amine])
    assert len(amides) > 0
