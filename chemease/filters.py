from rdkit import Chem

from chemease.properties import calculate_properties


def passes_lipinski(mol: Chem.Mol) -> bool:
    if mol is None:
        return False

    props = calculate_properties(mol)
    if not props:
        return False

    violations = 0
    if props.get("MW", 0) > 500:
        violations += 1
    if props.get("LogP", 0) > 5:
        violations += 1
    if props.get("NumHDonors", 0) > 5:
        violations += 1
    if props.get("NumHAcceptors", 0) > 10:
        violations += 1

    return violations <= 1


def filter_by_substructure(mols: list[Chem.Mol], smarts_pattern: str) -> list[Chem.Mol]:
    query = Chem.MolFromSmarts(smarts_pattern)
    if query is None:
        raise ValueError(f"Invalid SMARTS pattern: {smarts_pattern}")

    return [mol for mol in mols if mol is not None and mol.HasSubstructMatch(query)]
