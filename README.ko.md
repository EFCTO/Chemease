# ChemEase 🧪

[![PyPI version](https://badge.fury.io/py/chemease.svg)](https://badge.fury.io/py/chemease)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**ChemEase**는 RDKit, PubChem, ChEMBL을 기반으로 구축된 다용도 화학정보 라이브러리입니다.

[🇬🇧 English Documentation](./README.md)

## 주요 기능
- **데이터베이스 연동:** PubChem 및 ChEMBL에서 화합물 구조, 타겟 정보, 생물학적 활성 데이터를 손쉽게 수집합니다.
- **고급 Descriptor:** 기계 학습(ML) 모델링을 위해 물리화학적 특성 및 MACCS Keys와 같은 지문을 일괄 계산합니다.
- **3D 구조 최적화:** ETKDG 알고리즘으로 3D Conformer를 생성하고 MMFF94 Force Field로 에너지를 최적화합니다.
- **클러스터링 및 분석:** 대규모 화합물 라이브러리에 대해 Butina 클러스터링과 PCA 분석을 지원합니다.
- **가상 합성:** SMIRKS 패턴을 이용하여 아미드 커플링 등의 가상 화학 반응을 시뮬레이션합니다.
- **CLI 지원:** 파이썬 코드 없이 터미널 명령어만으로 핵심 작업들을 수행할 수 있습니다.

## 설치 방법

```bash
pip install chemease
```

## 빠른 시작 (CLI)

```bash
# 이름으로 화합물 검색 후 SDF 파일로 저장
chemease fetch --name "aspirin" -o aspirin.sdf

# 분자 특성(Descriptors) 계산 후 CSV 저장
chemease descriptors library.sdf -o features.csv

# 유사도 기반 화합물 클러스터링
chemease cluster library.sdf --cutoff 0.6
```

## 빠른 시작 (Python API)

```python
from chemease import get_mol_from_chembl, calculate_advanced_descriptors_df, generate_conformers

# ChEMBL에서 가져오기
mol = get_mol_from_chembl("CHEMBL25") # Aspirin

# 머신러닝 피처(Feature) 생성
df_features = calculate_advanced_descriptors_df([mol])

# 3D 구조 생성 및 최적화
mol_3d = generate_conformers(mol, num_confs=5, optimize=True)
```
