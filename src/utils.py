from __future__ import annotations

import logging
import re
import shutil
import subprocess
from pathlib import Path
from typing import Iterable


def get_logger(verbose: bool = True) -> logging.Logger:
    logger = logging.getLogger("sprb_scan")
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        logger.addHandler(handler)
    logger.setLevel(logging.INFO if verbose else logging.WARNING)
    logger.propagate = False
    return logger


def run_command(cmd: list[str], logger: logging.Logger | None = None) -> None:
    if logger:
        logger.info("CMD %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed:\n{' '.join(cmd)}\n\nSTDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"
        )


def remove_prefix_artifacts(prefix: Path) -> None:
    parent = prefix.parent
    stem = prefix.name
    for path in parent.glob(stem + "*"):
        if path.is_file():
            path.unlink()


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    normalized = sorted((min(a, b), max(a, b)) for a, b in intervals)
    if not normalized:
        return []

    merged = [normalized[0]]
    for start, end in normalized[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def extract_family(query_id: str) -> str:
    if "__" in query_id:
        return query_id.split("__", 1)[0]
    return query_id


def extract_member(query_id: str) -> str:
    if "__" in query_id:
        return query_id.split("__", 1)[1]
    return query_id


def short_module_name(module_name: str) -> str:
    return re.sub(r"^domain_", "", str(module_name))


def safe_filename(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def prepare_output_path(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def remove_tree(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)

