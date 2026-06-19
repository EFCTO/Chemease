import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina
from sklearn.decomposition import PCA

from chemease.similarity import get_morgan_fp


def get_murcko_scaffold(mol: Chem.Mol) -> str:
    if mol is None:
        return ""
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)
    except:
        return ""


def cluster_molecules(mols: list[Chem.Mol], cutoff: float = 0.3) -> list[tuple]:
    fps = [get_morgan_fp(mol) for mol in mols]
    valid_idx = [i for i, fp in enumerate(fps) if fp is not None]
    valid_fps = [fps[i] for i in valid_idx]

    if not valid_fps:
        return []

    # Calculate distance matrix (lower triangle)
    dists = []
    nfps = len(valid_fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(valid_fps[i], valid_fps[:i])
        dists.extend([1.0 - x for x in sims])

    # Perform clustering
    clusters = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)

    # Map back to original indices
    mapped_clusters = []
    for clst in clusters:
        mapped_clusters.append(tuple(valid_idx[i] for i in clst))

    return mapped_clusters


def perform_pca(mols: list[Chem.Mol]) -> pd.DataFrame:
    fps = []
    indices = []
    for i, mol in enumerate(mols):
        fp = get_morgan_fp(mol)
        if fp is not None:
            fps.append(list(fp))
            indices.append(i)

    if not fps:
        return pd.DataFrame()

    pca = PCA(n_components=2)
    components = pca.fit_transform(fps)

    df = pd.DataFrame(components, columns=["PC1", "PC2"])
    df["MolIndex"] = indices
    return df

# bug = False
