## Installation

### Using pip

```bash
cd biopilot
pip install -e .
```

### Using conda

```bash
conda env create -f biopilot/env/environment.yml
conda activate biopilot
```

### Install optional dependencies

```bash
pip install biopilot[dev]  # Development tools
pip install biopilot[r]    # R integration via rpy2
```

## Folder Structure

```
biopilot/
├── data/
│   ├── raw/           # Raw data files
│   ├── processed/     # Processed data
│   └── metadata/      # Metadata and annotations
├── results/
│   ├── figures/       # Generated plots
│   ├── tables/        # Result tables
│   ├── reports/       # Analysis reports
│   └── pipelines/     # Pipeline outputs
├── notebooks/         # Jupyter notebooks
├── scripts/           # Analysis scripts
├── logs/              # Log files and command history
├── env/               # Environment files (requirements.txt, environment.yml)
└── src/
    ├── fetcher/       # Dataset fetcher module
    ├── annotation/    # Annotation database module
    ├── pipeline/      # Pipeline manager module
    ├── analyzer/      # Data analyzer module
    ├── dashboard/     # Streamlit dashboard
    └── reproducibility/  # Reproducibility logger
```
