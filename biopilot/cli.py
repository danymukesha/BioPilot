#!/usr/bin/env python3
"""BioPilot CLI - Command line interface for bioinformatics workflows."""

import argparse
import json
import sys
from pathlib import Path

from biopilot.src.fetcher import DatasetFetcher
from biopilot.src.annotation import AnnotationDB
from biopilot.src.annotation.database import Sample
from biopilot.src.pipeline import PipelineManager
from biopilot.src.analyzer import DataAnalyzer
from biopilot.src.reproducibility import ReproducibilityLogger


def cmd_search(args):
    fetcher = DatasetFetcher(data_dir=args.data_dir)
    
    results = {}
    
    if args.database == "all":
        results = fetcher.search_all(args.query, args.max_results)
    elif args.database == "geo":
        results["geo"] = fetcher.search_geo(args.query, args.max_results)
    elif args.database == "sra":
        results["sra"] = fetcher.search_sra(args.query, args.max_results)
    elif args.database == "ena":
        results["ena"] = fetcher.search_ena(args.query, args.max_results)
    
    for db_name, datasets in results.items():
        print(f"\n=== {db_name.upper()} Results ===")
        for ds in datasets:
            print(f"\n{ds.accession}: {ds.title[:60]}...")
            print(f"  Organism: {ds.organism}")
            print(f"  Platform: {ds.platform}")
            print(f"  Type: {ds.study_type}")


def cmd_download(args):
    fetcher = DatasetFetcher(data_dir=args.data_dir)
    
    if args.source == "sra":
        result = fetcher.download_sra(args.accession)
        if result:
            print(f"Downloaded to: {result}")
    elif args.source == "ena":
        if args.url:
            result = fetcher.download_ena_fastq(args.url)
            if result:
                print(f"Downloaded to: {result}")
        else:
            print("Error: --url required for ENA downloads")
    elif args.source == "geo":
        ds = fetcher.get_geo_dataset(args.accession)
        if ds:
            fetcher.save_metadata(ds)
            print(f"Saved metadata for {args.accession}")


def cmd_add_sample(args):
    db = AnnotationDB(db_path=f"{args.data_dir}/metadata/annotations.db")
    
    sample = Sample(
        sample_id=args.sample_id,
        accession=args.accession or "",
        species=args.species,
        tissue_type=args.tissue or "",
        sequencing_type=args.seq_type,
        platform=args.platform or "",
        condition=args.condition or "",
        replicate=args.replicate,
        metadata=json.loads(args.metadata) if args.metadata else None,
    )
    
    db.add_sample(sample)
    print(f"Added sample: {args.sample_id}")


def cmd_list_samples(args):
    db = AnnotationDB(db_path=f"{args.data_dir}/metadata/annotations.db")
    
    filters = {}
    if args.species:
        filters["species"] = args.species
    if args.tissue:
        filters["tissue_type"] = args.tissue
    if args.seq_type:
        filters["sequencing_type"] = args.seq_type
    
    samples = db.query_samples(**filters)
    
    print(f"\nFound {len(samples)} samples:\n")
    print(f"{'Sample ID':<20} {'Species':<20} {'Tissue':<15} {'Condition':<15}")
    print("-" * 70)
    
    for s in samples:
        print(f"{s.sample_id:<20} {s.species:<20} {s.tissue_type:<15} {s.condition:<15}")


def cmd_stats(args):
    db = AnnotationDB(db_path=f"{args.data_dir}/metadata/annotations.db")
    stats = db.get_statistics()
    
    print("\n=== Database Statistics ===\n")
    print(f"Total samples: {stats.get('total_samples', 0)}")
    print(f"Species count: {stats.get('species_count', 0)}")
    
    if stats.get("by_species"):
        print("\nSamples by species:")
        for species, count in stats["by_species"].items():
            print(f"  {species}: {count}")
    
    if stats.get("by_sequencing_type"):
        print("\nSamples by sequencing type:")
        for seq_type, count in stats["by_sequencing_type"].items():
            print(f"  {seq_type}: {count}")


def cmd_run_pipeline(args):
    manager = PipelineManager(
        data_dir=args.data_dir,
        results_dir=args.results_dir,
        logs_dir=args.logs_dir,
    )
    
    if args.pipeline_type == "snakemake":
        config = {
            "data_dir": args.data_dir,
            "output_dir": f"{args.results_dir}/{args.name}",
            "samples": args.samples.split(",") if args.samples else [],
        }
        
        manager.create_snakemake_pipeline(args.name, config)
        
        run = manager.run_snakemake(
            args.name,
            samples=config["samples"],
            cores=args.cores,
            dry_run=args.dry_run,
        )
        
        print(f"Pipeline run: {run.run_id}")
        print(f"Status: {run.status}")
    
    elif args.pipeline_type == "nextflow":
        config = {
            "data_dir": args.data_dir,
            "output_dir": f"{args.results_dir}/{args.name}",
        }
        
        manager.create_nextflow_pipeline(args.name, config)
        
        run = manager.run_nextflow(args.name)
        
        print(f"Pipeline run: {run.run_id}")
        print(f"Status: {run.status}")


def cmd_analyze(args):
    analyzer = DataAnalyzer(results_dir=args.results_dir, logs_dir=args.logs_dir)
    
    if args.analysis == "normalize":
        data = analyzer.load_expression_matrix(args.input_file, sep=args.sep)
        normalized, params = analyzer.normalize(data, method=args.method)
        
        output_file, _ = analyzer.save_results(normalized, "normalization", params, args.input_file)
        print(f"Normalized data saved to: {output_file}")
    
    elif args.analysis == "pca":
        data = analyzer.load_expression_matrix(args.input_file, sep=args.sep)
        result, params = analyzer.pca(data, n_components=args.components)
        
        output_file, _ = analyzer.save_results(result, "pca", params, args.input_file)
        print(f"PCA results saved to: {output_file}")
        
        if not args.no_plot:
            fig_path = analyzer.plot_pca(result, params["explained_variance"][:2])
            if fig_path:
                print(f"PCA plot saved to: {fig_path}")
    
    elif args.analysis == "de":
        data = analyzer.load_expression_matrix(args.input_file, sep=args.sep)
        
        group1 = args.group1.split(",")
        group2 = args.group2.split(",")
        
        result, params = analyzer.differential_expression(
            data,
            {"group1": group1, "group2": group2},
            fold_change_threshold=args.fc_threshold,
            pvalue_threshold=args.pval_threshold,
        )
        
        output_file, _ = analyzer.save_results(result, "differential_expression", params, args.input_file)
        print(f"DE results saved to: {output_file}")
        print(f"Significant genes: {params['significant_genes']}")
        
        if not args.no_plot:
            fig_path = analyzer.plot_volcano(result)
            if fig_path:
                print(f"Volcano plot saved to: {fig_path}")


def cmd_dashboard(args):
    import subprocess
    
    dashboard_path = Path(__file__).parent / "src" / "dashboard" / "app.py"
    
    cmd = [
        sys.executable, "-m", "streamlit", "run",
        str(dashboard_path),
        "--server.port", str(args.port),
    ]
    
    if args.headless:
        cmd.extend(["--server.headless", "true"])
    
    print(f"Starting dashboard on port {args.port}...")
    print("Running command:", cmd)
    subprocess.run(cmd)


def cmd_env_capture(args):
    logger = ReproducibilityLogger(logs_dir=args.logs_dir, env_dir=args.env_dir)
    
    snapshot = logger.capture_environment(
        save_requirements=not args.no_requirements,
        save_conda_env=not args.no_conda,
    )
    
    print(f"\nEnvironment snapshot captured: {snapshot.hash}")
    print(f"Python version: {snapshot.python_version}")
    print(f"Platform: {snapshot.platform_info.get('system', 'Unknown')}")
    print(f"Total packages: {len(snapshot.packages)}")
    
    if not args.no_requirements:
        print(f"Requirements saved to: {args.env_dir}/requirements.txt")
    if not args.no_conda:
        print(f"Conda env saved to: {args.env_dir}/environment.yml")


def cmd_log_export(args):
    logger = ReproducibilityLogger(logs_dir=args.logs_dir, env_dir=args.env_dir)
    
    output_path = logger.export_log(args.output, format=args.format)
    print(f"Log exported to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        prog="biopilot",
        description="BioPilot - A lightweight bioinformatics assistant"
    )
    
    parser.add_argument("--data-dir", default="biopilot/data", help="Data directory")
    parser.add_argument("--results-dir", default="biopilot/results", help="Results directory")
    parser.add_argument("--logs-dir", default="biopilot/logs", help="Logs directory")
    parser.add_argument("--env-dir", default="biopilot/env", help="Environment directory")
    
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    search_parser = subparsers.add_parser("search", help="Search databases")
    search_parser.add_argument("query", help="Search query")
    search_parser.add_argument("--database", "-d", choices=["geo", "sra", "ena", "all"], default="all")
    search_parser.add_argument("--max-results", "-n", type=int, default=10)
    search_parser.set_defaults(func=cmd_search)
    
    download_parser = subparsers.add_parser("download", help="Download data")
    download_parser.add_argument("accession", help="Accession number")
    download_parser.add_argument("--source", "-s", choices=["sra", "ena", "geo"], default="sra")
    download_parser.add_argument("--url", help="Direct download URL")
    download_parser.set_defaults(func=cmd_download)
    
    add_sample_parser = subparsers.add_parser("add-sample", help="Add sample to database")
    add_sample_parser.add_argument("--sample-id", required=True)
    add_sample_parser.add_argument("--accession", default="")
    add_sample_parser.add_argument("--species", required=True)
    add_sample_parser.add_argument("--tissue", default="")
    add_sample_parser.add_argument("--seq-type", default="RNA-seq")
    add_sample_parser.add_argument("--platform", default="")
    add_sample_parser.add_argument("--condition", default="")
    add_sample_parser.add_argument("--replicate", type=int, default=None)
    add_sample_parser.add_argument("--metadata", help="JSON metadata")
    add_sample_parser.set_defaults(func=cmd_add_sample)
    
    list_parser = subparsers.add_parser("list", help="List samples")
    list_parser.add_argument("--species")
    list_parser.add_argument("--tissue")
    list_parser.add_argument("--seq-type")
    list_parser.set_defaults(func=cmd_list_samples)
    
    subparsers.add_parser("stats", help="Show database statistics").set_defaults(func=cmd_stats)
    
    pipeline_parser = subparsers.add_parser("pipeline", help="Run pipeline")
    pipeline_parser.add_argument("name", help="Pipeline name")
    pipeline_parser.add_argument("--pipeline-type", "-t", choices=["snakemake", "nextflow"], default="snakemake")
    pipeline_parser.add_argument("--samples", help="Comma-separated sample list")
    pipeline_parser.add_argument("--cores", "-j", type=int, default=4)
    pipeline_parser.add_argument("--dry-run", action="store_true")
    pipeline_parser.set_defaults(func=cmd_run_pipeline)
    
    analyze_parser = subparsers.add_parser("analyze", help="Run analysis")
    analyze_parser.add_argument("analysis", choices=["normalize", "pca", "de"])
    analyze_parser.add_argument("--input-file", "-i", required=True)
    analyze_parser.add_argument("--sep", default="\t")
    analyze_parser.add_argument("--method", default="cpm", help="Normalization method")
    analyze_parser.add_argument("--components", type=int, default=2)
    analyze_parser.add_argument("--group1", help="Group 1 samples (comma-separated)")
    analyze_parser.add_argument("--group2", help="Group 2 samples (comma-separated)")
    analyze_parser.add_argument("--fc-threshold", type=float, default=2.0)
    analyze_parser.add_argument("--pval-threshold", type=float, default=0.05)
    analyze_parser.add_argument("--no-plot", action="store_true")
    analyze_parser.set_defaults(func=cmd_analyze)
    
    dashboard_parser = subparsers.add_parser("dashboard", help="Launch web dashboard")
    dashboard_parser.add_argument("--port", "-p", type=int, default=8501)
    dashboard_parser.add_argument("--headless", action="store_true")
    dashboard_parser.set_defaults(func=cmd_dashboard)
    
    env_parser = subparsers.add_parser("env", help="Environment management")
    env_parser.add_argument("action", choices=["capture", "export"])
    env_parser.add_argument("--no-requirements", action="store_true")
    env_parser.add_argument("--no-conda", action="store_true")
    env_parser.add_argument("--output", "-o", default="biopilot/logs/export.json")
    env_parser.add_argument("--format", "-f", choices=["json", "txt"], default="json")
    env_parser.set_defaults(func=lambda args: cmd_env_capture(args) if args.action == "capture" else cmd_log_export(args))
    
    args = parser.parse_args()
    
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
