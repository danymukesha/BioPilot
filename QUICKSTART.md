## Quick Start

### 1. Search for datasets

```bash
biopilot search "RNA-seq human cancer" --database geo --max-results 10
biopilot search "Homo sapiens" --database ena
biopilot search "RNA-seq" --database sra
```

### 2. Download data

```bash
biopilot download SRR123456 --source sra
biopilot download GSE12345 --source geo
biopilot download ERR123456 --source ena --url "ftp://..."
```

### 3. Manage samples

```bash
biopilot add-sample --sample-id sample_001 --species "Homo sapiens" --tissue liver --seq-type RNA-seq --condition treated
biopilot list --species "Homo sapiens"
biopilot stats
```

### 4. Run analysis

```bash
biopilot analyze normalize -i expression_matrix.tsv --method cpm
biopilot analyze pca -i normalized.tsv --components 3
biopilot analyze de -i normalized.tsv --group1 sample1,sample2 --group2 control1,control2
```

### 5. Run pipelines

```bash
biopilot pipeline rnaseq --pipeline-type snakemake --samples sample1,sample2 --cores 4
biopilot pipeline rnaseq --pipeline-type snakemake --dry-run
```

### 6. Launch dashboard

```bash
biopilot dashboard --port 8501
```

### 7. Capture environment

```bash
biopilot env capture
biopilot env export --format json -o analysis_log.json
```
## Python API Usage

```python
from biopilot.src.fetcher import DatasetFetcher
from biopilot.src.annotation import AnnotationDB, Sample
from biopilot.src.analyzer import DataAnalyzer
from biopilot.src.reproducibility import ReproducibilityLogger

fetcher = DatasetFetcher(data_dir="data")
results = fetcher.search_geo("RNA-seq human", max_results=5)

db = AnnotationDB(db_path="data/metadata/annotations.db")
sample = Sample(
    sample_id="sample_001",
    accession="GSE12345",
    species="Homo sapiens",
    tissue_type="liver",
    sequencing_type="RNA-seq",
    condition="treated"
)
db.add_sample(sample)

analyzer = DataAnalyzer(results_dir="results")
data = analyzer.load_expression_matrix("expression.tsv")
normalized, params = analyzer.normalize(data, method="cpm")
pca_result, pca_params = analyzer.pca(normalized, n_components=2)

logger = ReproducibilityLogger(logs_dir="logs")
logger.capture_environment()
```

## Example Workflow

```bash
biopilot search "GSE121212" --database geo
biopilot download GSE121212 --source geo
biopilot add-sample --sample-id GSM123456 --species "Homo sapiens" --tissue liver --seq-type RNA-seq
biopilot analyze normalize -i data/raw/expression.tsv --method cpm
biopilot analyze pca -i results/tables/analysis_xxx.csv
biopilot analyze de -i results/tables/analysis_xxx.csv --group1 treated1,treated2 --group2 control1,control2
biopilot dashboard
```

CLI commands

```bash
biopilot search "RNA-seq human cancer" --database geo
biopilot download SRR123456 --source sra
biopilot add-sample --sample-id s1 --species "Homo sapiens" --tissue liver
biopilot list --species "Homo sapiens"
biopilot analyze normalize -i expression.tsv --method cpm
biopilot analyze pca -i normalized.tsv
biopilot analyze de -i normalized.tsv --group1 treated1,treated2 --group2 ctrl1,ctrl2
biopilot dashboard --port 8501
biopilot env capture
```

