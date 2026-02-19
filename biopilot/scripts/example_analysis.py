#!/usr/bin/env python3
"""Example analysis script demonstrating BioPilot usage."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np

from biopilot.src.annotation import AnnotationDB, Sample
from biopilot.src.analyzer import DataAnalyzer
from biopilot.src.reproducibility import ReproducibilityLogger


def main():
    print("=" * 60)
    print("BioPilot Example Analysis Script")
    print("=" * 60)
    
    db = AnnotationDB(db_path="data/metadata/annotations.db")
    analyzer = DataAnalyzer(results_dir="results", logs_dir="logs")
    logger = ReproducibilityLogger(logs_dir="logs", env_dir="env")
    
    print("\n1. Adding samples to database...")
    samples = [
        Sample(
            sample_id=f"sample_{i}",
            accession=f"ACC{i:03d}",
            species="Homo sapiens",
            tissue_type="liver" if i < 3 else "brain",
            sequencing_type="RNA-seq",
            platform="Illumina",
            condition="treated" if i % 2 == 0 else "control",
        )
        for i in range(6)
    ]
    
    for sample in samples:
        db.add_sample(sample)
    
    stats = db.get_statistics()
    print(f"   Total samples: {stats['total_samples']}")
    
    print("\n2. Generating example expression data...")
    np.random.seed(42)
    
    n_genes = 500
    sample_ids = [s.sample_id for s in samples]
    
    data = pd.DataFrame(
        np.random.negative_binomial(20, 0.5, size=(n_genes, len(sample_ids))),
        columns=sample_ids,
        index=[f"gene_{i}" for i in range(n_genes)],
    )
    
    print(f"   Expression matrix: {data.shape[0]} genes x {data.shape[1]} samples")
    
    print("\n3. Normalizing data (CPM + log2)...")
    normalized, norm_params = analyzer.normalize(data, method="cpm", log_transform=True)
    print(f"   Normalization complete")
    
    print("\n4. Running PCA...")
    pca_result, pca_params = analyzer.pca(normalized, n_components=2)
    print(f"   PC1: {pca_params['explained_variance'][0]*100:.1f}% variance")
    print(f"   PC2: {pca_params['explained_variance'][1]*100:.1f}% variance")
    
    print("\n5. Performing differential expression analysis...")
    treated = [s.sample_id for s in samples if s.condition == "treated"]
    control = [s.sample_id for s in samples if s.condition == "control"]
    
    de_result, de_params = analyzer.differential_expression(
        normalized,
        {"treated": treated, "control": control},
        fold_change_threshold=1.5,
        pvalue_threshold=0.05,
    )
    print(f"   Significant genes: {de_params['significant_genes']}")
    
    print("\n6. Saving results...")
    output_file, analysis = analyzer.save_results(
        de_result, "differential_expression", de_params
    )
    print(f"   Results saved to: {output_file}")
    
    print("\n7. Capturing environment...")
    snapshot = logger.capture_environment()
    print(f"   Environment hash: {snapshot.hash}")
    
    logger.log_command(
        "python scripts/example_analysis.py",
        exit_code=0,
        duration=5.0,
    )
    
    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
