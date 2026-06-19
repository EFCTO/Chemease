from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs


def get_morgan_fp(mol: Chem.Mol, radius: int = 2, nBits: int = 2048):
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)


def calculate_similarity(mol1: Chem.Mol, mol2: Chem.Mol) -> float:
    fp1 = get_morgan_fp(mol1)
    fp2 = get_morgan_fp(mol2)

    if fp1 is None or fp2 is None:
        return 0.0

    return DataStructs.TanimotoSimilarity(fp1, fp2)


def find_similar(
    query_mol: Chem.Mol, target_mols: list[Chem.Mol], threshold: float = 0.7
) -> list[tuple[Chem.Mol, float]]:
    query_fp = get_morgan_fp(query_mol)
    if query_fp is None:
        return []

    results = []
    for mol in target_mols:
        if mol is None:
            continue
        fp = get_morgan_fp(mol)
        if fp is not None:
            sim = DataStructs.TanimotoSimilarity(query_fp, fp)
            if sim >= threshold:
                results.append((mol, sim))

    # Sort by similarity descending
    results.sort(key=lambda x: x[1], reverse=True)
    return results

# going to sleep now
