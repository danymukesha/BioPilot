"""Reproducibility logger for tracking commands, versions, and environment."""

import json
import logging
import platform
import subprocess
import sys
from pathlib import Path
from typing import Optional, List, Dict, Any
from datetime import datetime
from dataclasses import dataclass, asdict
import hashlib

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class CommandLog:
    command: str
    timestamp: str
    working_dir: str
    exit_code: Optional[int] = None
    duration_seconds: Optional[float] = None
    output: Optional[str] = None
    error: Optional[str] = None


@dataclass
class EnvironmentSnapshot:
    timestamp: str
    python_version: str
    platform_info: Dict[str, str]
    packages: Dict[str, str]
    environment_file: Optional[str] = None
    hash: Optional[str] = None


class ReproducibilityLogger:
    """Log commands, library versions, and environment snapshots."""
    
    def __init__(self, logs_dir: str = "logs", env_dir: str = "env"):
        self.logs_dir = Path(logs_dir)
        self.env_dir = Path(env_dir)
        
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        self.env_dir.mkdir(parents=True, exist_ok=True)
        
        self.commands_file = self.logs_dir / "commands.json"
        self.snapshot_file = self.logs_dir / "environment_snapshot.json"
        self.session_file = self.logs_dir / "session.json"
        
        self.commands: List[CommandLog] = self._load_commands()
        self.session_start = datetime.now().isoformat()
        self.session_commands: List[str] = []
    
    def _load_commands(self) -> List[CommandLog]:
        if not self.commands_file.exists():
            return []
        with open(self.commands_file, "r", encoding="utf-8") as f:
            data = json.load(f)
        return [CommandLog(**cmd) for cmd in data]
    
    def _save_commands(self):
        with open(self.commands_file, "w", encoding="utf-8") as f:
            json.dump([asdict(cmd) for cmd in self.commands], f, indent=2, default=str)
    
    def log_command(
        self,
        command: str,
        working_dir: Optional[str] = None,
        exit_code: Optional[int] = None,
        duration: Optional[float] = None,
        output: Optional[str] = None,
        error: Optional[str] = None,
    ) -> CommandLog:
        log_entry = CommandLog(
            command=command,
            timestamp=datetime.now().isoformat(),
            working_dir=working_dir or str(Path.cwd()),
            exit_code=exit_code,
            duration_seconds=duration,
            output=output[:1000] if output else None,
            error=error[:1000] if error else None,
        )
        
        self.commands.append(log_entry)
        self.session_commands.append(command)
        self._save_commands()
        
        logger.info(f"Logged command: {command[:50]}...")
        return log_entry
    
    def get_python_version(self) -> str:
        return f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    
    def get_platform_info(self) -> Dict[str, str]:
        return {
            "system": platform.system(),
            "node": platform.node(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
            "processor": platform.processor(),
        }
    
    def get_installed_packages(self) -> Dict[str, str]:
        packages = {}
        
        try:
            import importlib.metadata
            
            for dist in importlib.metadata.distributions():
                name = dist.metadata.get("Name", "")
                version = dist.version
                if name:
                    packages[name] = version
        except Exception:
            try:
                result = subprocess.run(
                    [sys.executable, "-m", "pip", "list", "--format=json"],
                    capture_output=True,
                    text=True,
                )
                if result.returncode == 0:
                    for pkg in json.loads(result.stdout):
                        packages[pkg["name"]] = pkg["version"]
            except Exception as e:
                logger.warning(f"Could not get package list: {e}")
        
        return packages
    
    def capture_environment(
        self,
        include_packages: bool = True,
        save_requirements: bool = True,
        save_conda_env: bool = True,
    ) -> EnvironmentSnapshot:
        snapshot = EnvironmentSnapshot(
            timestamp=datetime.now().isoformat(),
            python_version=self.get_python_version(),
            platform_info=self.get_platform_info(),
            packages=self.get_installed_packages() if include_packages else {},
        )
        
        content = json.dumps(asdict(snapshot), sort_keys=True)
        snapshot.hash = hashlib.sha256(content.encode()).hexdigest()[:16]
        
        with open(self.snapshot_file, "w", encoding="utf-8") as f:
            json.dump(asdict(snapshot), f, indent=2, default=str)
        
        if save_requirements:
            self._save_requirements(snapshot.packages)
            snapshot.environment_file = str(self.env_dir / "requirements.txt")
        
        if save_conda_env:
            self._save_conda_env(snapshot.packages)
        
        logger.info(f"Captured environment snapshot: {snapshot.hash}")
        return snapshot
    
    def _save_requirements(self, packages: Dict[str, str]):
        requirements_file = self.env_dir / "requirements.txt"
        
        with open(requirements_file, "w", encoding="utf-8") as f:
            for name, version in sorted(packages.items()):
                f.write(f"{name}=={version}\n")
        
        logger.info(f"Saved requirements to {requirements_file}")
    
    def _save_conda_env(self, packages: Dict[str, str]):
        env_file = self.env_dir / "environment.yml"
        
        content = f"""name: biopilot
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python={self.get_python_version()}
  - pip
  - pip:
"""
        for name, version in sorted(packages.items()):
            content += f"    - {name}=={version}\n"
        
        with open(env_file, "w", encoding="utf-8") as f:
            f.write(content)
        
        logger.info(f"Saved conda environment to {env_file}")
    
    def load_last_snapshot(self) -> Optional[EnvironmentSnapshot]:
        if not self.snapshot_file.exists():
            return None
        
        with open(self.snapshot_file, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        return EnvironmentSnapshot(**data)
    
    def compare_environments(
        self,
        snapshot1: EnvironmentSnapshot,
        snapshot2: EnvironmentSnapshot,
    ) -> Dict[str, Any]:
        pkgs1 = snapshot1.packages
        pkgs2 = snapshot2.packages
        
        all_packages = set(pkgs1.keys()) | set(pkgs2.keys())
        
        differences = {
            "added": {},
            "removed": {},
            "changed": {},
            "unchanged": [],
        }
        
        for pkg in all_packages:
            v1 = pkgs1.get(pkg)
            v2 = pkgs2.get(pkg)
            
            if v1 is None and v2 is not None:
                differences["added"][pkg] = v2
            elif v1 is not None and v2 is None:
                differences["removed"][pkg] = v1
            elif v1 != v2:
                differences["changed"][pkg] = {"from": v1, "to": v2}
            else:
                differences["unchanged"].append(pkg)
        
        return differences
    
    def get_session_log(self) -> Dict[str, Any]:
        return {
            "session_start": self.session_start,
            "session_end": datetime.now().isoformat(),
            "commands_executed": len(self.session_commands),
            "commands": self.session_commands,
        }
    
    def save_session_log(self) -> Path:
        session_log = self.get_session_log()
        
        with open(self.session_file, "w", encoding="utf-8") as f:
            json.dump(session_log, f, indent=2, default=str)
        
        logger.info(f"Saved session log to {self.session_file}")
        return self.session_file
    
    def get_commands_in_range(
        self,
        start_time: Optional[str] = None,
        end_time: Optional[str] = None,
    ) -> List[CommandLog]:
        filtered = []
        
        for cmd in self.commands:
            if start_time and cmd.timestamp < start_time:
                continue
            if end_time and cmd.timestamp > end_time:
                continue
            filtered.append(cmd)
        
        return filtered
    
    def export_log(self, output_path: str, format: str = "json") -> Path:
        output_file = Path(output_path)
        
        if format == "json":
            data = {
                "session": self.get_session_log(),
                "commands": [asdict(cmd) for cmd in self.commands],
                "last_snapshot": asdict(self.load_last_snapshot()) if self.load_last_snapshot() else None,
            }
            with open(output_file, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2, default=str)
        
        elif format == "txt":
            with open(output_file, "w", encoding="utf-8") as f:
                f.write("BioPilot Command Log\n")
                f.write("=" * 50 + "\n\n")
                f.write(f"Session started: {self.session_start}\n\n")
                f.write("Commands:\n")
                f.write("-" * 50 + "\n")
                for cmd in self.commands:
                    f.write(f"[{cmd.timestamp}] {cmd.command}\n")
                    if cmd.exit_code is not None:
                        f.write(f"  Exit code: {cmd.exit_code}\n")
                    f.write("\n")
        
        logger.info(f"Exported log to {output_file}")
        return output_file
    
    def generate_report(self) -> str:
        snapshot = self.load_last_snapshot()
        
        report = []
        report.append("# BioPilot Reproducibility Report")
        report.append(f"\nGenerated: {datetime.now().isoformat()}")
        report.append(f"\n## Session Summary")
        report.append(f"- Session started: {self.session_start}")
        report.append(f"- Commands executed: {len(self.session_commands)}")
        
        if snapshot:
            report.append(f"\n## Environment Snapshot")
            report.append(f"- Timestamp: {snapshot.timestamp}")
            report.append(f"- Python: {snapshot.python_version}")
            report.append(f"- Platform: {snapshot.platform_info.get('system', 'Unknown')}")
            report.append(f"- Total packages: {len(snapshot.packages)}")
        
        report.append(f"\n## Commands")
        for i, cmd in enumerate(self.commands[-20:], 1):
            status = "✓" if cmd.exit_code == 0 else "✗" if cmd.exit_code else "?"
            report.append(f"{i}. [{status}] {cmd.command}")
        
        return "\n".join(report)
