"""Pipeline manager for Snakemake and Nextflow workflows."""

import os
import json
import logging
import subprocess
import shutil
from pathlib import Path
from typing import Optional, List, Dict, Any
from datetime import datetime
from dataclasses import dataclass, asdict

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class PipelineRun:
    run_id: str
    pipeline_name: str
    status: str
    start_time: str
    end_time: Optional[str] = None
    config: Optional[Dict[str, Any]] = None
    samples: Optional[List[str]] = None
    output_dir: Optional[str] = None
    log_file: Optional[str] = None
    error_message: Optional[str] = None


class PipelineManager:
    """Manage Snakemake and Nextflow pipelines."""
    
    def __init__(
        self,
        data_dir: str = "data",
        results_dir: str = "results",
        logs_dir: str = "logs",
    ):
        self.data_dir = Path(data_dir)
        self.results_dir = Path(results_dir)
        self.logs_dir = Path(logs_dir)
        self.pipelines_dir = self.results_dir / "pipelines"
        self.intermediate_dir = self.results_dir / "intermediate"
        
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        self.pipelines_dir.mkdir(parents=True, exist_ok=True)
        self.intermediate_dir.mkdir(parents=True, exist_ok=True)
        
        self.runs_file = self.logs_dir / "pipeline_runs.json"
        self.runs: List[PipelineRun] = self._load_runs()
        
        self._check_tools()
    
    def _check_tools(self):
        self.tools = {}
        
        for tool in ["snakemake", "nextflow"]:
            result = subprocess.run(
                ["which", tool],
                capture_output=True,
                text=True,
            )
            self.tools[tool] = result.returncode == 0
        
        logger.info(f"Available tools: {self.tools}")
    
    def _load_runs(self) -> List[PipelineRun]:
        if not self.runs_file.exists():
            return []
        
        with open(self.runs_file, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        return [PipelineRun(**run) for run in data]
    
    def _save_runs(self):
        with open(self.runs_file, "w", encoding="utf-8") as f:
            json.dump([asdict(run) for run in self.runs], f, indent=2, default=str)
    
    def _generate_run_id(self) -> str:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return f"run_{timestamp}"
    
    def create_snakemake_pipeline(
        self,
        name: str,
        config: Dict[str, Any],
        snakefile_content: Optional[str] = None,
    ) -> Path:
        pipeline_dir = self.pipelines_dir / name
        pipeline_dir.mkdir(parents=True, exist_ok=True)
        
        config_file = pipeline_dir / "config.yaml"
        import yaml
        with open(config_file, "w", encoding="utf-8") as f:
            yaml.dump(config, f, default_flow_style=False)
        
        if snakefile_content:
            snakefile = pipeline_dir / "Snakefile"
            snakefile.write_text(snakefile_content, encoding="utf-8")
        else:
            self._create_default_snakefile(pipeline_dir, config)
        
        logger.info(f"Created Snakemake pipeline: {name}")
        return pipeline_dir
    
    def _create_default_snakefile(self, pipeline_dir: Path, config: Dict[str, Any]):
        snakefile_content = '''
configfile: "config.yaml"

rule all:
    input:
        expand(config["output_dir"] + "/{sample}_processed.fastq", sample=config["samples"])

rule quality_check:
    input:
        config["data_dir"] + "/raw/{sample}.fastq"
    output:
        config["output_dir"] + "/qc/{sample}_fastqc.html"
    shell:
        "fastqc {input} -o {wildcards.output_dir}/qc/"

rule trim:
    input:
        config["data_dir"] + "/raw/{sample}.fastq"
    output:
        config["output_dir"] + "/{sample}_processed.fastq"
    params:
        adapter=config.get("adapter", "AGATCGGAAGAGC")
    shell:
        "cutadapt -a {params.adapter} -o {output} {input}"
'''
        snakefile = pipeline_dir / "Snakefile"
        snakefile.write_text(snakefile_content.strip(), encoding="utf-8")
    
    def create_nextflow_pipeline(
        self,
        name: str,
        config: Dict[str, Any],
        main_nf_content: Optional[str] = None,
    ) -> Path:
        pipeline_dir = self.pipelines_dir / name
        pipeline_dir.mkdir(parents=True, exist_ok=True)
        
        config_file = pipeline_dir / "nextflow.config"
        config_content = f'''
params {{
    data_dir = "{config.get('data_dir', 'data/raw')}"
    output_dir = "{config.get('output_dir', 'results')}"
    samples = {config.get('samples', [])}
}}
'''
        config_file.write_text(config_content.strip(), encoding="utf-8")
        
        if main_nf_content:
            main_nf = pipeline_dir / "main.nf"
            main_nf.write_text(main_nf_content, encoding="utf-8")
        else:
            self._create_default_nextflow(pipeline_dir)
        
        logger.info(f"Created Nextflow pipeline: {name}")
        return pipeline_dir
    
    def _create_default_nextflow(self, pipeline_dir: Path):
        main_nf_content = '''
#!/usr/bin/env nextflow

params.samples = []
params.data_dir = "data/raw"
params.output_dir = "results"

process FASTQC {
    input:
    path fastq
    
    output:
    path "*_fastqc.html"
    
    script:
    """
    fastqc $fastq
    """
}

process TRIM {
    input:
    path fastq
    
    output:
    path "*_trimmed.fastq"
    
    script:
    """
    cutadapt -a AGATCGGAAGAGC -o \$(basename $fastq .fastq)_trimmed.fastq $fastq
    """
}

workflow {
    samples = Channel.fromPath(params.data_dir + "/*.fastq")
    
    FASTQC(samples)
    TRIM(samples)
}
'''
        main_nf = pipeline_dir / "main.nf"
        main_nf.write_text(main_nf_content.strip(), encoding="utf-8")
    
    def run_snakemake(
        self,
        pipeline_name: str,
        samples: Optional[List[str]] = None,
        cores: int = 4,
        dry_run: bool = False,
        config_override: Optional[Dict[str, Any]] = None,
    ) -> PipelineRun:
        pipeline_dir = self.pipelines_dir / pipeline_name
        snakefile = pipeline_dir / "Snakefile"
        
        if not snakefile.exists():
            raise FileNotFoundError(f"Snakefile not found: {snakefile}")
        
        run_id = self._generate_run_id()
        output_dir = self.results_dir / pipeline_name / run_id
        output_dir.mkdir(parents=True, exist_ok=True)
        
        log_file = self.logs_dir / f"{run_id}.log"
        
        run = PipelineRun(
            run_id=run_id,
            pipeline_name=pipeline_name,
            status="running",
            start_time=datetime.now().isoformat(),
            samples=samples,
            output_dir=str(output_dir),
            log_file=str(log_file),
        )
        
        if config_override:
            config_file = pipeline_dir / "config.yaml"
            import yaml
            with open(config_file, "r", encoding="utf-8") as f:
                config = yaml.safe_load(f)
            config.update(config_override)
            config["output_dir"] = str(output_dir)
            if samples:
                config["samples"] = samples
            with open(config_file, "w", encoding="utf-8") as f:
                yaml.dump(config, f, default_flow_style=False)
        
        cmd = [
            "snakemake",
            "-s", str(snakefile),
            "-j", str(cores),
            "--directory", str(pipeline_dir),
        ]
        
        if dry_run:
            cmd.append("--dryrun")
        
        self.runs.append(run)
        self._save_runs()
        
        try:
            with open(log_file, "w", encoding="utf-8") as log:
                result = subprocess.run(
                    cmd,
                    cwd=str(pipeline_dir),
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
            
            run.status = "completed" if result.returncode == 0 else "failed"
            run.end_time = datetime.now().isoformat()
            
            if result.returncode != 0:
                run.error_message = f"Exit code: {result.returncode}"
        
        except Exception as e:
            run.status = "error"
            run.end_time = datetime.now().isoformat()
            run.error_message = str(e)
        
        self._save_runs()
        logger.info(f"Pipeline run {run_id}: {run.status}")
        return run
    
    def run_nextflow(
        self,
        pipeline_name: str,
        profiles: Optional[List[str]] = None,
        config_override: Optional[Dict[str, Any]] = None,
    ) -> PipelineRun:
        pipeline_dir = self.pipelines_dir / pipeline_name
        main_nf = pipeline_dir / "main.nf"
        
        if not main_nf.exists():
            raise FileNotFoundError(f"main.nf not found: {main_nf}")
        
        run_id = self._generate_run_id()
        output_dir = self.results_dir / pipeline_name / run_id
        output_dir.mkdir(parents=True, exist_ok=True)
        
        log_file = self.logs_dir / f"{run_id}.log"
        
        run = PipelineRun(
            run_id=run_id,
            pipeline_name=pipeline_name,
            status="running",
            start_time=datetime.now().isoformat(),
            output_dir=str(output_dir),
            log_file=str(log_file),
        )
        
        self.runs.append(run)
        self._save_runs()
        
        cmd = ["nextflow", "run", str(pipeline_dir)]
        
        if profiles:
            cmd.extend(["-profile", ",".join(profiles)])
        
        if config_override:
            for key, value in config_override.items():
                cmd.append(f"--{key}={value}")
        
        try:
            with open(log_file, "w", encoding="utf-8") as log:
                result = subprocess.run(
                    cmd,
                    cwd=str(pipeline_dir),
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
            
            run.status = "completed" if result.returncode == 0 else "failed"
            run.end_time = datetime.now().isoformat()
            
            if result.returncode != 0:
                run.error_message = f"Exit code: {result.returncode}"
        
        except Exception as e:
            run.status = "error"
            run.end_time = datetime.now().isoformat()
            run.error_message = str(e)
        
        self._save_runs()
        logger.info(f"Pipeline run {run_id}: {run.status}")
        return run
    
    def list_pipelines(self) -> List[str]:
        if not self.pipelines_dir.exists():
            return []
        
        pipelines = []
        for item in self.pipelines_dir.iterdir():
            if item.is_dir():
                if (item / "Snakefile").exists() or (item / "main.nf").exists():
                    pipelines.append(item.name)
        
        return pipelines
    
    def get_run(self, run_id: str) -> Optional[PipelineRun]:
        for run in self.runs:
            if run.run_id == run_id:
                return run
        return None
    
    def list_runs(self, pipeline_name: Optional[str] = None) -> List[PipelineRun]:
        if pipeline_name:
            return [r for r in self.runs if r.pipeline_name == pipeline_name]
        return self.runs
    
    def get_intermediate_files(self, run_id: str) -> List[Path]:
        run = self.get_run(run_id)
        if not run or not run.output_dir:
            return []
        
        output_dir = Path(run.output_dir)
        files = []
        
        for pattern in ["**/*.bam", "**/*.sam", "**/*.fastq", "**/*.fq", "**/*.vcf"]:
            files.extend(output_dir.glob(pattern))
        
        return files
    
    def cleanup_intermediate(self, run_id: str, keep_final: bool = True) -> int:
        run = self.get_run(run_id)
        if not run or not run.output_dir:
            return 0
        
        output_dir = Path(run.output_dir)
        removed = 0
        
        intermediate_patterns = ["**/*.bam", "**/*.sam", "**/tmp_*"]
        
        for pattern in intermediate_patterns:
            for file_path in output_dir.glob(pattern):
                if keep_final and file_path.name.startswith("final_"):
                    continue
                if file_path.is_file():
                    file_path.unlink()
                    removed += 1
        
        logger.info(f"Cleaned up {removed} intermediate files for {run_id}")
        return removed
    
    def get_pipeline_status(self, pipeline_name: str) -> Dict[str, Any]:
        runs = [r for r in self.runs if r.pipeline_name == pipeline_name]
        
        if not runs:
            return {"status": "never_run", "runs": 0}
        
        latest = runs[-1]
        
        return {
            "status": latest.status,
            "latest_run": latest.run_id,
            "total_runs": len(runs),
            "completed": sum(1 for r in runs if r.status == "completed"),
            "failed": sum(1 for r in runs if r.status in ["failed", "error"]),
        }
    
    def export_pipeline(self, pipeline_name: str, output_path: str) -> Path:
        pipeline_dir = self.pipelines_dir / pipeline_name
        output_file = Path(output_path)
        
        shutil.make_archive(
            str(output_file.with_suffix("")),
            "zip",
            pipeline_dir,
        )
        
        logger.info(f"Exported pipeline to {output_file}")
        return output_file
    
    def delete_pipeline(self, pipeline_name: str, delete_results: bool = False) -> bool:
        pipeline_dir = self.pipelines_dir / pipeline_name
        
        if not pipeline_dir.exists():
            return False
        
        shutil.rmtree(pipeline_dir)
        
        if delete_results:
            results_dir = self.results_dir / pipeline_name
            if results_dir.exists():
                shutil.rmtree(results_dir)
        
        self.runs = [r for r in self.runs if r.pipeline_name != pipeline_name]
        self._save_runs()
        
        logger.info(f"Deleted pipeline: {pipeline_name}")
        return True
