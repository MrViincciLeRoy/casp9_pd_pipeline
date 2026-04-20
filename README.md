# Caspase-9 × Parkinson's Disease — Drug Discovery Pipeline

Computational pipeline to identify novel Caspase-9 inhibitors as potential neuroprotectants for Parkinson's Disease.

## Pipeline

```
Phase 1 → Phase 2 → Phase 3 → Phase 4 → Phase 5
Data      Docking    ML         ADMET      Report
```

| Phase | What it does | Key tools |
|-------|-------------|-----------|
| 1 | Downloads 2AR9.pdb, ChEMBL Ki data, ZINC subset | BioPython, requests |
| 2 | Virtual docking of compounds against Caspase-9 | AutoDock Vina, GNINA |
| 3 | Fine-tunes ChemBERTa-2 on Ki data, screens ZINC | HuggingFace, PyTorch |
| 4 | Lipinski + BBB permeability filtering | RDKit, pkCSM API |
| 5 | Generates ranked candidate report (Markdown) | — |

## Install

```bash
pip install -r requirements.txt
# AutoDock Vina: conda install -c conda-forge vina
# GNINA: https://github.com/gnina/gnina/releases
```

## Run

```bash
# Full pipeline (skips Phase 2 — docking needs manual ligand prep)
python run_pipeline.py

# Specific phases
python run_pipeline.py --phases 1 3 4 5

# Phase 2 (docking) — manual, needs pdbqt files prepared first
python run_pipeline.py --phases 2
```

## Project Structure

```
casp9_pd_pipeline/
├── run_pipeline.py        # orchestrator
├── config.py              # all settings
├── phase1_data/           # PDB, ChEMBL, ZINC fetchers
├── phase2_docking/        # Vina + GNINA docking
├── phase3_ml/             # ChemBERTa model + screener
├── phase4_admet/          # Lipinski + BBB filter
├── phase5_report/         # Markdown report generator
├── utils/                 # logger
└── data/                  # downloaded/generated data
    ├── structures/         # .pdb and .pdbqt files
    ├── chembl/             # inhibitor Ki data
    ├── zinc/               # compound library
    └── results/            # hits, model, reports
```

## Output

`data/results/report_TIMESTAMP.md` — ranked candidates with Vina scores, predicted Ki, and BBB status.
