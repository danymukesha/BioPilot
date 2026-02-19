"""Test suite for BioPilot."""

import pytest
import tempfile
import shutil
from pathlib import Path

from biopilot.src.annotation import AnnotationDB, Sample
from biopilot.src.analyzer import DataAnalyzer
from biopilot.src.reproducibility import ReproducibilityLogger


class TestAnnotationDB:
    @pytest.fixture
    def temp_db(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = Path(tmpdir) / "test.db"
            db = AnnotationDB(db_path=str(db_path))
            yield db
    
    def test_add_sample(self, temp_db):
        sample = Sample(
            sample_id="test_001",
            accession="TEST001",
            species="Homo sapiens",
            tissue_type="liver",
            sequencing_type="RNA-seq",
            platform="Illumina",
            condition="treated",
        )
        
        row_id = temp_db.add_sample(sample)
        assert row_id > 0
    
    def test_get_sample(self, temp_db):
        sample = Sample(
            sample_id="test_002",
            accession="TEST002",
            species="Mus musculus",
            tissue_type="brain",
            sequencing_type="RNA-seq",
            platform="Illumina",
            condition="control",
        )
        temp_db.add_sample(sample)
        
        retrieved = temp_db.get_sample("test_002")
        assert retrieved is not None
        assert retrieved.species == "Mus musculus"
        assert retrieved.tissue_type == "brain"
    
    def test_query_samples(self, temp_db):
        samples = [
            Sample(sample_id="q1", accession="A1", species="Homo sapiens", tissue_type="liver", sequencing_type="RNA-seq", platform="Illumina", condition="treated"),
            Sample(sample_id="q2", accession="A2", species="Homo sapiens", tissue_type="brain", sequencing_type="RNA-seq", platform="Illumina", condition="control"),
            Sample(sample_id="q3", accession="A3", species="Mus musculus", tissue_type="liver", sequencing_type="RNA-seq", platform="Illumina", condition="treated"),
        ]
        
        for s in samples:
            temp_db.add_sample(s)
        
        human_samples = temp_db.query_samples(species="Homo sapiens")
        assert len(human_samples) == 2
        
        liver_samples = temp_db.query_samples(tissue_type="liver")
        assert len(liver_samples) == 2
    
    def test_add_tag(self, temp_db):
        sample = Sample(
            sample_id="tag_test",
            accession="T1",
            species="Homo sapiens",
            tissue_type="liver",
            sequencing_type="RNA-seq",
            platform="Illumina",
            condition="treated",
        )
        temp_db.add_sample(sample)
        
        temp_db.add_tag("tag_test", "batch", "batch_001")
        tags = temp_db.get_tags("tag_test")
        
        assert "batch" in tags
        assert tags["batch"] == "batch_001"
    
    def test_statistics(self, temp_db):
        samples = [
            Sample(sample_id="s1", accession="A1", species="Homo sapiens", tissue_type="liver", sequencing_type="RNA-seq", platform="Illumina", condition="treated"),
            Sample(sample_id="s2", accession="A2", species="Homo sapiens", tissue_type="brain", sequencing_type="RNA-seq", platform="Illumina", condition="control"),
            Sample(sample_id="s3", accession="A3", species="Mus musculus", tissue_type="liver", sequencing_type="ChIP-seq", platform="Illumina", condition="treated"),
        ]
        
        for s in samples:
            temp_db.add_sample(s)
        
        stats = temp_db.get_statistics()
        assert stats["total_samples"] == 3
        assert stats["species_count"] == 2


class TestDataAnalyzer:
    @pytest.fixture
    def temp_dirs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            results_dir = Path(tmpdir) / "results"
            logs_dir = Path(tmpdir) / "logs"
            analyzer = DataAnalyzer(
                results_dir=str(results_dir),
                logs_dir=str(logs_dir),
            )
            yield analyzer
    
    def test_normalize_cpm(self, temp_dirs):
        import pandas as pd
        import numpy as np
        
        data = pd.DataFrame(
            np.random.randint(100, 1000, size=(10, 4)),
            columns=["sample1", "sample2", "sample3", "sample4"],
            index=[f"gene_{i}" for i in range(10)],
        )
        
        normalized, params = temp_dirs.normalize(data, method="cpm", log_transform=False)
        
        assert normalized.shape == data.shape
        assert params["method"] == "cpm"
        assert "log2_transformed" not in params
    
    def test_normalize_zscore(self, temp_dirs):
        import pandas as pd
        import numpy as np
        
        data = pd.DataFrame(
            np.random.randn(10, 4),
            columns=["sample1", "sample2", "sample3", "sample4"],
            index=[f"gene_{i}" for i in range(10)],
        )
        
        normalized, params = temp_dirs.normalize(data, method="zscore")
        
        assert normalized.shape == data.shape
        assert params["method"] == "zscore"


class TestReproducibilityLogger:
    @pytest.fixture
    def temp_logger(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            logs_dir = Path(tmpdir) / "logs"
            env_dir = Path(tmpdir) / "env"
            logger = ReproducibilityLogger(
                logs_dir=str(logs_dir),
                env_dir=str(env_dir),
            )
            yield logger
    
    def test_log_command(self, temp_logger):
        log = temp_logger.log_command(
            command="biopilot search test",
            exit_code=0,
            duration=1.5,
        )
        
        assert log.command == "biopilot search test"
        assert log.exit_code == 0
        assert log.duration_seconds == 1.5
    
    def test_capture_environment(self, temp_logger):
        snapshot = temp_logger.capture_environment(
            save_requirements=False,
            save_conda_env=False,
        )
        
        assert snapshot.python_version is not None
        assert snapshot.platform_info is not None
        assert snapshot.hash is not None
    
    def test_get_commands(self, temp_logger):
        temp_logger.log_command("cmd1", exit_code=0)
        temp_logger.log_command("cmd2", exit_code=1)
        
        commands = temp_logger.commands
        assert len(commands) == 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
