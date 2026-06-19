# ChemEase 🧪

[![PyPI version](https://badge.fury.io/py/chemease.svg)](https://badge.fury.io/py/chemease)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**ChemEase** is a multi-purpose cheminformatics library built on top of RDKit, PubChem, and ChEMBL.

[🇰🇷 한국어 문서(Korean Readme) 읽기](./README.ko.md)

## Features
- **Database Integration:** Fetch molecules, targets, and bioactivities easily from PubChem and ChEMBL.
- **Advanced Descriptors:** Calculate physicochemical properties and advanced fingerprints (e.g., MACCS Keys) for machine learning.
- **3D Structure & Optimization:** Generate 3D conformers using ETKDG and optimize them with MMFF94 force fields.
- **Clustering & Analysis:** Perform Butina clustering and PCA on large chemical libraries.
- **Virtual Synthesis:** Simulate chemical reactions (e.g., Amide coupling) using SMIRKS/SMARTS patterns.
- **CLI Support:** Run major workflows directly from your terminal.

## Installation

```bash
pip install chemease
```

## Quick Start (CLI)

```bash
# Fetch a molecule and save it as SDF
chemease fetch --name "aspirin" -o aspirin.sdf

# Calculate advanced descriptors
chemease descriptors library.sdf -o features.csv

# Cluster a library
chemease cluster library.sdf --cutoff 0.6
```

## Quick Start (Python)

```python
from chemease import get_mol_from_chembl, calculate_advanced_descriptors_df, generate_conformers

# Load from ChEMBL
mol = get_mol_from_chembl("CHEMBL25") # Aspirin

# Calculate ML Features
df_features = calculate_advanced_descriptors_df([mol])

# Generate 3D Conformers
mol_3d = generate_conformers(mol, num_confs=5, optimize=True)
```

## Contributing
Contributions are welcome! Please run `black`, `isort`, `flake8`, and `pytest` before submitting a pull request.
