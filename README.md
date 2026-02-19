# BioPilot

Bioinformatics assistant for personal use on a laptop, to act as:

- **Dataset fetcher**: query GEO, NCBI SRA, and ENA databases; download raw 
data and metadata
- **Annotation module**: store sample metadata in SQLite database with 
species, tissue, and sequencing type tags
- **Pipeline manager**: wrap Snakemake/Nextflow pipelines with organized 
folder structure
- **Data analyzer**: normalization (CPM/RPKM/TPM), PCA, differential 
expression analysis
- **Interactive dashboard**: streamlit-based web interface for exploration 
and visualization
- **Reproducibility module**: log commands, library versions, timestamps; 
save environment snapshots

## Installation

[DEV INSTALL document](./INSTALL.md)

Requirements

- Python 3.10+
- 8-16 GB RAM (for typical RNA-seq datasets)
- Optional: Snakemake or Nextflow for pipeline execution
- Optional: SRA Toolkit for SRA downloads

Development

```bash
pip install -e ".[dev]"
pytest biopilot/tests
black biopilot/
ruff check biopilot/
mypy biopilot/
```

## Quick start

[QUICK START document](./QUICKSTART.md)

## Supported Analyses

[SUPPORT document](./SUPPORT.md)

## License

MIT License
