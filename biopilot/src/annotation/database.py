"""SQLite-based annotation database for sample metadata."""

import json
import sqlite3
import logging
from pathlib import Path
from typing import Optional, List, Dict, Any
from datetime import datetime
from dataclasses import dataclass

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class Sample:
    sample_id: str
    accession: str
    species: str
    tissue_type: str
    sequencing_type: str
    platform: str
    condition: str
    replicate: Optional[int] = None
    metadata: Optional[Dict[str, Any]] = None
    created_at: Optional[str] = None
    updated_at: Optional[str] = None


class AnnotationDB:
    """SQLite database for storing and querying sample annotations."""
    
    SCHEMA = """
    CREATE TABLE IF NOT EXISTS samples (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        sample_id TEXT UNIQUE NOT NULL,
        accession TEXT,
        species TEXT,
        tissue_type TEXT,
        sequencing_type TEXT,
        platform TEXT,
        condition TEXT,
        replicate INTEGER,
        metadata TEXT,
        created_at TEXT,
        updated_at TEXT
    );
    
    CREATE INDEX IF NOT EXISTS idx_species ON samples(species);
    CREATE INDEX IF NOT EXISTS idx_tissue ON samples(tissue_type);
    CREATE INDEX IF NOT EXISTS idx_seq_type ON samples(sequencing_type);
    CREATE INDEX IF NOT EXISTS idx_accession ON samples(accession);
    CREATE INDEX IF NOT EXISTS idx_condition ON samples(condition);
    
    CREATE TABLE IF NOT EXISTS tags (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        sample_id TEXT NOT NULL,
        tag_name TEXT NOT NULL,
        tag_value TEXT,
        created_at TEXT,
        FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
        UNIQUE(sample_id, tag_name)
    );
    
    CREATE INDEX IF NOT EXISTS idx_tag_name ON tags(tag_name);
    
    CREATE TABLE IF NOT EXISTS projects (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        project_id TEXT UNIQUE NOT NULL,
        name TEXT,
        description TEXT,
        samples TEXT,
        created_at TEXT,
        updated_at TEXT
    );
    
    CREATE TABLE IF NOT EXISTS analysis_log (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        sample_id TEXT,
        analysis_type TEXT,
        parameters TEXT,
        results TEXT,
        timestamp TEXT
    );
    """
    
    def __init__(self, db_path: str = "data/metadata/annotations.db"):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = None
        self._connect()
        self._init_db()
    
    def _connect(self):
        #self.conn = sqlite3.connect(str(self.db_path))
        self.conn = sqlite3.connect(str(self.db_path), check_same_thread=False)
        self.conn.row_factory = sqlite3.Row
    
    def _init_db(self):
        self.conn.executescript(self.SCHEMA)
        self.conn.commit()
        logger.info(f"Initialized database at {self.db_path}")
    
    def add_sample(self, sample: Sample) -> int:
        now = datetime.now().isoformat()
        
        cursor = self.conn.execute(
            """
            INSERT OR REPLACE INTO samples 
            (sample_id, accession, species, tissue_type, sequencing_type, 
             platform, condition, replicate, metadata, created_at, updated_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                sample.sample_id,
                sample.accession,
                sample.species,
                sample.tissue_type,
                sample.sequencing_type,
                sample.platform,
                sample.condition,
                sample.replicate,
                json.dumps(sample.metadata) if sample.metadata else None,
                sample.created_at or now,
                now,
            ),
        )
        self.conn.commit()
        return cursor.lastrowid
    
    def add_samples_batch(self, samples: List[Sample]) -> int:
        now = datetime.now().isoformat()
        count = 0
        
        for sample in samples:
            self.conn.execute(
                """
                INSERT OR REPLACE INTO samples 
                (sample_id, accession, species, tissue_type, sequencing_type, 
                 platform, condition, replicate, metadata, created_at, updated_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    sample.sample_id,
                    sample.accession,
                    sample.species,
                    sample.tissue_type,
                    sample.sequencing_type,
                    sample.platform,
                    sample.condition,
                    sample.replicate,
                    json.dumps(sample.metadata) if sample.metadata else None,
                    sample.created_at or now,
                    now,
                ),
            )
            count += 1
        
        self.conn.commit()
        logger.info(f"Added {count} samples to database")
        return count
    
    def get_sample(self, sample_id: str) -> Optional[Sample]:
        row = self.conn.execute(
            "SELECT * FROM samples WHERE sample_id = ?", (sample_id,)
        ).fetchone()
        
        if not row:
            return None
        
        return self._row_to_sample(row)
    
    def _row_to_sample(self, row: sqlite3.Row) -> Sample:
        return Sample(
            sample_id=row["sample_id"],
            accession=row["accession"],
            species=row["species"],
            tissue_type=row["tissue_type"],
            sequencing_type=row["sequencing_type"],
            platform=row["platform"],
            condition=row["condition"],
            replicate=row["replicate"],
            metadata=json.loads(row["metadata"]) if row["metadata"] else None,
            created_at=row["created_at"],
            updated_at=row["updated_at"],
        )
    
    def query_samples(self, **filters) -> List[Sample]:
        conditions = []
        params = []
        
        for key, value in filters.items():
            if value is not None and hasattr(Sample, key):
                conditions.append(f"{key} = ?")
                params.append(value)
        
        where_clause = " AND ".join(conditions) if conditions else "1=1"
        query = f"SELECT * FROM samples WHERE {where_clause}"
        
        rows = self.conn.execute(query, params).fetchall()
        return [self._row_to_sample(row) for row in rows]
    
    def search_samples(self, text: str) -> List[Sample]:
        query = """
            SELECT * FROM samples 
            WHERE sample_id LIKE ? 
               OR accession LIKE ?
               OR species LIKE ?
               OR tissue_type LIKE ?
               OR condition LIKE ?
        """
        pattern = f"%{text}%"
        rows = self.conn.execute(query, [pattern] * 5).fetchall()
        return [self._row_to_sample(row) for row in rows]
    
    def add_tag(self, sample_id: str, tag_name: str, tag_value: str = "") -> int:
        now = datetime.now().isoformat()
        
        cursor = self.conn.execute(
            """
            INSERT OR REPLACE INTO tags (sample_id, tag_name, tag_value, created_at)
            VALUES (?, ?, ?, ?)
            """,
            (sample_id, tag_name, tag_value, now),
        )
        self.conn.commit()
        return cursor.lastrowid
    
    def get_tags(self, sample_id: str) -> Dict[str, str]:
        rows = self.conn.execute(
            "SELECT tag_name, tag_value FROM tags WHERE sample_id = ?",
            (sample_id,),
        ).fetchall()
        
        return {row["tag_name"]: row["tag_value"] for row in rows}
    
    def get_samples_by_tag(self, tag_name: str, tag_value: Optional[str] = None) -> List[Sample]:
        if tag_value:
            query = """
                SELECT s.* FROM samples s
                JOIN tags t ON s.sample_id = t.sample_id
                WHERE t.tag_name = ? AND t.tag_value = ?
            """
            rows = self.conn.execute(query, (tag_name, tag_value)).fetchall()
        else:
            query = """
                SELECT s.* FROM samples s
                JOIN tags t ON s.sample_id = t.sample_id
                WHERE t.tag_name = ?
            """
            rows = self.conn.execute(query, (tag_name,)).fetchall()
        
        return [self._row_to_sample(row) for row in rows]
    
    def list_species(self) -> List[str]:
        rows = self.conn.execute(
            "SELECT DISTINCT species FROM samples WHERE species IS NOT NULL"
        ).fetchall()
        return [row["species"] for row in rows]
    
    def list_tissue_types(self) -> List[str]:
        rows = self.conn.execute(
            "SELECT DISTINCT tissue_type FROM samples WHERE tissue_type IS NOT NULL"
        ).fetchall()
        return [row["tissue_type"] for row in rows]
    
    def list_sequencing_types(self) -> List[str]:
        rows = self.conn.execute(
            "SELECT DISTINCT sequencing_type FROM samples WHERE sequencing_type IS NOT NULL"
        ).fetchall()
        return [row["sequencing_type"] for row in rows]
    
    def get_statistics(self) -> Dict[str, Any]:
        stats = {}
        
        stats["total_samples"] = self.conn.execute(
            "SELECT COUNT(*) FROM samples"
        ).fetchone()[0]
        
        stats["species_count"] = self.conn.execute(
            "SELECT COUNT(DISTINCT species) FROM samples"
        ).fetchone()[0]
        
        stats["by_species"] = dict(
            self.conn.execute(
                "SELECT species, COUNT(*) FROM samples GROUP BY species"
            ).fetchall()
        )
        
        stats["by_tissue"] = dict(
            self.conn.execute(
                "SELECT tissue_type, COUNT(*) FROM samples GROUP BY tissue_type"
            ).fetchall()
        )
        
        stats["by_sequencing_type"] = dict(
            self.conn.execute(
                "SELECT sequencing_type, COUNT(*) FROM samples GROUP BY sequencing_type"
            ).fetchall()
        )
        
        return stats
    
    def create_project(self, project_id: str, name: str, description: str = "") -> int:
        now = datetime.now().isoformat()
        
        cursor = self.conn.execute(
            """
            INSERT OR REPLACE INTO projects (project_id, name, description, samples, created_at, updated_at)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (project_id, name, description, "[]", now, now),
        )
        self.conn.commit()
        return cursor.lastrowid
    
    def add_samples_to_project(self, project_id: str, sample_ids: List[str]) -> None:
        row = self.conn.execute(
            "SELECT samples FROM projects WHERE project_id = ?", (project_id,)
        ).fetchone()
        
        if not row:
            raise ValueError(f"Project {project_id} not found")
        
        current_samples = json.loads(row["samples"])
        current_samples.extend(sample_ids)
        unique_samples = list(set(current_samples))
        
        now = datetime.now().isoformat()
        self.conn.execute(
            "UPDATE projects SET samples = ?, updated_at = ? WHERE project_id = ?",
            (json.dumps(unique_samples), now, project_id),
        )
        self.conn.commit()
    
    def get_project_samples(self, project_id: str) -> List[Sample]:
        row = self.conn.execute(
            "SELECT samples FROM projects WHERE project_id = ?", (project_id,)
        ).fetchone()
        
        if not row:
            return []
        
        sample_ids = json.loads(row["samples"])
        samples = []
        
        for sample_id in sample_ids:
            sample = self.get_sample(sample_id)
            if sample:
                samples.append(sample)
        
        return samples
    
    def delete_sample(self, sample_id: str) -> bool:
        self.conn.execute("DELETE FROM tags WHERE sample_id = ?", (sample_id,))
        cursor = self.conn.execute(
            "DELETE FROM samples WHERE sample_id = ?", (sample_id,)
        )
        self.conn.commit()
        return cursor.rowcount > 0
    
    def export_to_csv(self, output_path: str) -> Path:
        import csv
        
        output_file = Path(output_path)
        rows = self.conn.execute("SELECT * FROM samples").fetchall()
        
        if not rows:
            logger.warning("No samples to export")
            return output_file
        
        with open(output_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            for row in rows:
                writer.writerow(dict(row))
        
        logger.info(f"Exported {len(rows)} samples to {output_file}")
        return output_file
    
    def close(self):
        if self.conn:
            self.conn.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
