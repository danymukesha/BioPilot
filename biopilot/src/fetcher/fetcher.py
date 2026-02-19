"""Dataset fetcher for GEO, NCBI, and ENA databases."""

import os
import json
import time
import hashlib
import logging
import subprocess
from pathlib import Path
from typing import Optional, List, Dict, Any
from datetime import datetime
from dataclasses import dataclass, asdict
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
import urllib.parse

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class DatasetInfo:
    accession: str
    title: str
    platform: str
    organism: str
    sample_count: int
    study_type: str
    source_db: str
    download_url: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None
    download_date: Optional[str] = None
    local_path: Optional[str] = None


class DatasetFetcher:
    """Fetch datasets from GEO, NCBI SRA, and ENA databases."""
    
    BASE_URLS = {
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi",
        "ncbi_esearch": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        "ncbi_esummary": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        "ncbi_efetch": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        "ena": "https://www.ebi.ac.uk/ena/portal/api",
    }
    
    def __init__(self, data_dir: str = "data", cache_dir: Optional[str] = None):
        self.data_dir = Path(data_dir)
        self.raw_dir = self.data_dir / "raw"
        self.metadata_dir = self.data_dir / "metadata"
        self.cache_dir = Path(cache_dir) if cache_dir else self.data_dir / ".cache"
        
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        self.metadata_dir.mkdir(parents=True, exist_ok=True)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        self.request_delay = 0.34
        self.last_request_time = 0
    
    def _rate_limit(self):
        elapsed = time.time() - self.last_request_time
        if elapsed < self.request_delay:
            time.sleep(self.request_delay - elapsed)
        self.last_request_time = time.time()
    
    def _fetch_url(self, url: str, params: Optional[Dict] = None) -> str:
        self._rate_limit()
        
        full_url = url
        if params:
            #param_str = "&".join(f"{k}={v}" for k, v in params.items())
            #full_url = f"{url}?{param_str}"
            full_url = f"{url}?{urllib.parse.urlencode(params)}"
        cache_key = hashlib.md5(full_url.encode()).hexdigest()
        cache_file = self.cache_dir / f"{cache_key}.txt"
        
        if cache_file.exists():
            return cache_file.read_text(encoding="utf-8")
        
        try:
            with urllib.request.urlopen(full_url, timeout=60) as response:
                content = response.read().decode("utf-8")
            cache_file.write_text(content, encoding="utf-8")
            return content
        except urllib.error.URLError as e:
            logger.error(f"Failed to fetch URL {full_url}: {e}")
            raise
    
    def search_geo(self, query: str, max_results: int = 20) -> List[DatasetInfo]:
        """Search GEO database for datasets."""
        logger.info(f"Searching GEO for: {query}")
        
        params = {
            "db": "gds",
            "term": query,
            "retmax": max_results,
            "retmode": "json",
        }
        
        try:
            content = self._fetch_url(self.BASE_URLS["ncbi_esearch"], params)
            data = json.loads(content)
            
            if "esearchresult" not in data:
                logger.warning("No results found in GEO search")
                return []
            
            ids = data["esearchresult"].get("idlist", [])
            return self._fetch_geo_summaries(ids)
        except Exception as e:
            logger.error(f"GEO search failed: {e}")
            return []
    
    def _fetch_geo_summaries(self, ids: List[str]) -> List[DatasetInfo]:
        if not ids:
            return []
        
        params = {
            "db": "gds",
            "id": ",".join(ids),
            "retmode": "json",
        }
        
        content = self._fetch_url(self.BASE_URLS["ncbi_esummary"], params)
        data = json.loads(content)
        
        datasets = []
        result = data.get("result", {})
        
        for uid in ids:
            if uid in result and uid != "uids":
                item = result[uid]
                datasets.append(DatasetInfo(
                    accession=item.get("accession", "N/A"),
                    title=item.get("title", "N/A"),
                    platform=item.get("platfm", {}).get("name", "Unknown"),
                    organism=item.get("taxon", "Unknown"),
                    sample_count=item.get("n_samples", 0),
                    study_type=item.get("gdstype", "Unknown"),
                    source_db="GEO",
                    metadata=item,
                ))
        
        return datasets
    
    def get_geo_dataset(self, accession: str) -> Optional[DatasetInfo]:
        """Get detailed information about a GEO dataset."""
        logger.info(f"Fetching GEO dataset: {accession}")
        
        params = {
            "db": "gds",
            "term": accession,
            "retmax": 1,
        }
        
        content = self._fetch_url(self.BASE_URLS["ncbi_esearch"], params)
        data = json.loads(content)
        
        ids = data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            logger.warning(f"Dataset {accession} not found")
            return None
        
        results = self._fetch_geo_summaries(ids)
        return results[0] if results else None
    
    def search_sra(self, query: str, max_results: int = 20) -> List[DatasetInfo]:
        """Search NCBI SRA for sequencing data."""
        logger.info(f"Searching SRA for: {query}")
        
        params = {
            "db": "sra",
            "term": query,
            "retmax": max_results,
        }
        
        content = self._fetch_url(self.BASE_URLS["ncbi_esearch"], params)
        data = json.loads(content)
        
        ids = data.get("esearchresult", {}).get("idlist", [])
        return self._fetch_sra_summaries(ids)
    
    def _fetch_sra_summaries(self, ids: List[str]) -> List[DatasetInfo]:
        if not ids:
            return []
        
        params = {
            "db": "sra",
            "id": ",".join(ids[:10]),
            "rettype": "runinfo",
            "retmode": "text",
        }
        
        try:
            content = self._fetch_url(self.BASE_URLS["ncbi_efetch"], params)
            return self._parse_sra_runinfo(content)
        except Exception as e:
            logger.error(f"Failed to fetch SRA summaries: {e}")
            return []
    
    def _parse_sra_runinfo(self, content: str) -> List[DatasetInfo]:
        datasets = []
        lines = content.strip().split("\n")
        
        if len(lines) < 2:
            return datasets
        
        headers = lines[0].split(",")
        
        for line in lines[1:]:
            values = line.split(",")
            row = dict(zip(headers, values))
            
            datasets.append(DatasetInfo(
                accession=row.get("Run", "N/A"),
                title=row.get("SampleName", "N/A"),
                platform=row.get("Platform", "Unknown"),
                organism=row.get("ScientificName", "Unknown"),
                sample_count=1,
                study_type=row.get("LibraryStrategy", "Unknown"),
                source_db="SRA",
                download_url=row.get("download_path"),
                metadata=row,
            ))
        
        return datasets
    
    def search_ena(self, query: str, max_results: int = 20) -> List[DatasetInfo]:
        """Search ENA database for datasets."""
        logger.info(f"Searching ENA for: {query}")
        
        url = f"{self.BASE_URLS['ena']}/search"
        params = {
            "dataPortal": "ena",
            "query": query,
            "result": "read_run",
            "fields": "run_accession,study_title,scientific_name,instrument_model,library_strategy,fastq_ftp",
            "format": "json",
            "limit": max_results,
        }
        
        try:
            content = self._fetch_url(url, params)
            data = json.loads(content)
            return self._parse_ena_results(data)
        except Exception as e:
            logger.error(f"ENA search failed: {e}")
            return []
    
    def _parse_ena_results(self, data: List[Dict]) -> List[DatasetInfo]:
        datasets = []
        
        for item in data:
            datasets.append(DatasetInfo(
                accession=item.get("run_accession", "N/A"),
                title=item.get("study_title", "N/A"),
                platform=item.get("instrument_model", "Unknown"),
                organism=item.get("scientific_name", "Unknown"),
                sample_count=1,
                study_type=item.get("library_strategy", "Unknown"),
                source_db="ENA",
                download_url=item.get("fastq_ftp"),
                metadata=item,
            ))
        
        return datasets
    
    def download_sra(self, accession: str, output_dir: Optional[str] = None) -> Optional[Path]:
        """Download SRA data using prefetch/fasterq-dump."""
        output_path = Path(output_dir) if output_dir else self.raw_dir
        output_path.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Downloading SRA: {accession}")
        
        try:
            result = subprocess.run(
                ["which", "prefetch"],
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                logger.warning("SRA Toolkit not found. Please install it for SRA downloads.")
                return None
            
            subprocess.run(
                ["prefetch", accession, "-O", str(output_path)],
                check=True,
                capture_output=True,
            )
            
            sra_file = output_path / accession / f"{accession}.sra"
            
            subprocess.run(
                ["fasterq-dump", str(sra_file), "-O", str(output_path)],
                check=True,
                capture_output=True,
            )
            
            logger.info(f"Successfully downloaded {accession}")
            return output_path
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to download {accession}: {e}")
            return None
    
    def download_ena_fastq(self, url: str, output_dir: Optional[str] = None) -> Optional[Path]:
        """Download FASTQ from ENA FTP server."""
        if not url:
            logger.error("No download URL provided")
            return None
        
        output_path = Path(output_dir) if output_dir else self.raw_dir
        output_path.mkdir(parents=True, exist_ok=True)
        
        filename = url.split("/")[-1]
        output_file = output_path / filename
        
        logger.info(f"Downloading from ENA: {filename}")
        
        try:
            self._rate_limit()
            urllib.request.urlretrieve(url, output_file)
            logger.info(f"Successfully downloaded {filename}")
            return output_file
        except Exception as e:
            logger.error(f"Failed to download {url}: {e}")
            return None
    
    def download_geo_supplementary(self, accession: str, output_dir: Optional[str] = None) -> Optional[Path]:
        """Download GEO supplementary files."""
        output_path = Path(output_dir) if output_dir else self.raw_dir
        output_path.mkdir(parents=True, exist_ok=True)
        
        geo_acc = accession.replace("GDS", "GSE")
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_acc[:-3]}nnn/{geo_acc}/suppl/"
        
        logger.info(f"Attempting to download GEO supplementary files for {accession}")
        logger.info(f"URL: {url}")
        
        return output_path
    
    def save_metadata(self, dataset: DatasetInfo, output_dir: Optional[str] = None) -> Path:
        """Save dataset metadata to JSON file."""
        output_path = Path(output_dir) if output_dir else self.metadata_dir
        output_path.mkdir(parents=True, exist_ok=True)
        
        dataset.download_date = datetime.now().isoformat()
        
        metadata_file = output_path / f"{dataset.accession}_metadata.json"
        
        with open(metadata_file, "w", encoding="utf-8") as f:
            json.dump(asdict(dataset), f, indent=2, default=str)
        
        logger.info(f"Saved metadata to {metadata_file}")
        return metadata_file
    
    def load_metadata(self, accession: str) -> Optional[DatasetInfo]:
        """Load previously saved metadata."""
        metadata_file = self.metadata_dir / f"{accession}_metadata.json"
        
        if not metadata_file.exists():
            return None
        
        with open(metadata_file, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        return DatasetInfo(**data)
    
    def list_downloaded(self) -> List[Dict[str, Any]]:
        """List all downloaded datasets."""
        datasets = []
        
        for metadata_file in self.metadata_dir.glob("*_metadata.json"):
            with open(metadata_file, "r", encoding="utf-8") as f:
                datasets.append(json.load(f))
        
        return datasets
    
    def search_all(self, query: str, max_results: int = 10) -> Dict[str, List[DatasetInfo]]:
        """Search all databases simultaneously."""
        return {
            "geo": self.search_geo(query, max_results),
            "sra": self.search_sra(query, max_results),
            "ena": self.search_ena(query, max_results),
        }
