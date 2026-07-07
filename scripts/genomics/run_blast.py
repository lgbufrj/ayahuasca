"""
run_blastp.py
─────────────────────────────────────────────────────────────────────────────
Runs BLASTP searches for every (protein × reference-species × target-organism)
combination, skipping jobs whose output already exists.

Features
  • Structured logging with timestamps and colour-coded levels
  • Parallel execution via concurrent.futures.ThreadPoolExecutor
  • Per-job retry logic with configurable attempts
  • Dry-run mode (--dry-run) to preview work without touching the filesystem
  • Summary table printed on completion
  • Clean separation of concerns: discovery → planning → execution → reporting
─────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import Iterator

from data import GENOME_PATH, PROTEINS_PATH, organisms, proteins
from logger import build_logger, COLORS

# ──────────────────────────────────────────────────────────────────────────────
# Logging
# ──────────────────────────────────────────────────────────────────────────────

log = build_logger("blastp")

RESET  = COLORS["RESET"]
BOLD   = COLORS["BOLD"]
GREEN  = COLORS["GREEN"]
YELLOW = COLORS["YELLOW"]
RED    = COLORS["RED"]
CYAN   = COLORS["CYAN"]
GREY   = COLORS["GREY"]

# ──────────────────────────────────────────────────────────────────────────────
# Data model
# ──────────────────────────────────────────────────────────────────────────────

class JobStatus(Enum):
    PENDING  = auto()
    SKIPPED  = auto()
    SUCCESS  = auto()
    FAILED   = auto()


@dataclass
class BlastJob:
    """All paths and metadata needed for a single BLASTP call."""

    prot_name:   str
    ref_species: str
    fasta_name:  str
    organism:    str
    genome_type: str
    query_path:  Path
    db_path:     Path
    out_path:    Path

    # Populated after execution
    status:   JobStatus = field(default=JobStatus.PENDING, compare=False)
    duration: float     = field(default=0.0,               compare=False)
    error:    str       = field(default="",                compare=False)

    @property
    def label(self) -> str:
        return f"{self.ref_species}/{self.prot_name} → {self.organism}"


# ──────────────────────────────────────────────────────────────────────────────
# Job discovery
# ──────────────────────────────────────────────────────────────────────────────

def _iter_ref_proteins() -> Iterator[tuple[str, str, str]]:
    """Yield (prot_name, species, fasta_name) for every reference protein."""
    for prot_name, info in proteins.items():
        for species, org_info in info["organisms"].items():
            uniprot    = org_info["uniprot_id"]
            fasta_name = f"{prot_name}_{uniprot}"
            yield prot_name, species, fasta_name


def _target_organisms() -> list[str]:
    """Return organism keys that are NOT reference genomes."""
    return [
        org
        for org, meta in organisms.items()
        if not meta["ref"]
    ]


def _db_path(organism: str, genome_type: str) -> Path:
    meta = organisms[organism]
    if genome_type == "phased":
        if not meta.get("genome_files", {}).get("prot", {}).get("phased"):
            raise ValueError(
                f"Organism '{organism}' has no phased genome. "
                f"Use --genome-type non_phased."
            )
        return Path(GENOME_PATH) / organism / "phased"     / "blast" / "prot_phased_db"
    return     Path(GENOME_PATH) / organism / "non_phased" / "blast" / "prot_non_phased_db"


def discover_jobs(*, genome_type: str = "phased") -> list[BlastJob]:
    """Build the full list of BLASTP jobs from the data dictionaries."""
    jobs: list[BlastJob] = []

    for prot_name, species, fasta_name in _iter_ref_proteins():
        query_path = (
            Path(PROTEINS_PATH)
            / prot_name / "reference" / species / "structure"
            / f"{fasta_name}.fasta"
        )

        for organism in _target_organisms():
            out_path = (
                Path(PROTEINS_PATH)
                / prot_name / "analysis" / "genomic" / "blast"
                / organism / species
                / f"{fasta_name}.xml"
            )

            jobs.append(
                BlastJob(
                    prot_name   = prot_name,
                    ref_species = species,
                    fasta_name  = fasta_name,
                    organism    = organism,
                    genome_type = genome_type,
                    query_path  = query_path,
                    db_path     = _db_path(organism, genome_type),
                    out_path    = out_path,
                )
            )

    return jobs


# ──────────────────────────────────────────────────────────────────────────────
# Execution
# ──────────────────────────────────────────────────────────────────────────────

BLASTP_DEFAULTS: dict[str, str] = {
    "-outfmt": "5",   # XML output
}


def run_job(
    job:       BlastJob,
    *,
    dry_run:   bool = False,
    overwrite: bool = False,
    max_retry: int  = 2,
) -> BlastJob:
    """
    Execute (or skip) a single BLASTP job in-place, returning the same object
    with its status, duration, and error fields populated.
    """
    # ── Already done? ────────────────────────────────────────────────────────
    if job.out_path.exists() and not overwrite:
        log.debug("skip  %s (output exists)", job.label)
        job.status = JobStatus.SKIPPED
        return job

    # ── Validate inputs ──────────────────────────────────────────────────────
    if not job.query_path.exists():
        job.status = JobStatus.FAILED
        job.error  = f"Query file not found: {job.query_path}"
        log.error("%s — %s", job.label, job.error)
        return job

    # ── Prepare output directory ─────────────────────────────────────────────
    if not dry_run:
        job.out_path.parent.mkdir(parents=True, exist_ok=True)

    # ── Build command ────────────────────────────────────────────────────────
    cmd = [
        "blastp",
        "-query", str(job.query_path),
        "-db",    str(job.db_path),
        "-out",   str(job.out_path),
        *[token for pair in BLASTP_DEFAULTS.items() for token in pair],
    ]

    if dry_run:
        log.info("dry-run  %s\n         %s", job.label, " ".join(cmd))
        job.status = JobStatus.SUCCESS
        return job

    # ── Run with retry ───────────────────────────────────────────────────────
    log.info("start  %s", job.label)
    t0 = time.perf_counter()

    for attempt in range(1, max_retry + 1):
        try:
            subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
            )
            job.status   = JobStatus.SUCCESS
            job.duration = time.perf_counter() - t0
            log.info(
                "done   %s  (%.1fs)",
                job.label,
                job.duration,
            )
            return job

        except subprocess.CalledProcessError as exc:
            if attempt < max_retry:
                log.warning(
                    "retry %d/%d  %s — %s",
                    attempt, max_retry, job.label, exc.stderr.strip(),
                )
                time.sleep(2 ** attempt)  # exponential back-off
            else:
                job.status = JobStatus.FAILED
                job.error  = exc.stderr.strip() or str(exc)
                job.duration = time.perf_counter() - t0
                log.error("FAIL   %s — %s", job.label, job.error)

    return job


def run_all(
    jobs:      list[BlastJob],
    *,
    workers:   int  = 4,
    dry_run:   bool = False,
    overwrite: bool = False,
    max_retry: int  = 2,
) -> list[BlastJob]:
    """Run all jobs concurrently and return them with updated statuses."""
    if not jobs:
        log.warning("No jobs to run.")
        return jobs

    log.info(
        "%s%s total jobs%s — workers=%d  dry_run=%s",
        COLORS['BOLD'], len(jobs), COLORS['RESET'], workers, dry_run,
    )

    completed: list[BlastJob] = []

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(run_job, job, dry_run=dry_run, overwrite=overwrite, max_retry=max_retry): job
            for job in jobs
        }

        for future in as_completed(futures):
            try:
                completed.append(future.result())
            except Exception as exc:          # pragma: no cover
                original = futures[future]
                original.status = JobStatus.FAILED
                original.error  = str(exc)
                log.error("Unexpected error for %s: %s", original.label, exc)
                completed.append(original)

    return completed


# ──────────────────────────────────────────────────────────────────────────────
# Reporting
# ──────────────────────────────────────────────────────────────────────────────

def print_summary(jobs: list[BlastJob]) -> None:
    """Print a coloured summary table to stdout."""
    counts = {s: 0 for s in JobStatus}
    for job in jobs:
        counts[job.status] += 1

    failed = [j for j in jobs if j.status == JobStatus.FAILED]
    total_time = sum(j.duration for j in jobs)

    width = 60
    sep   = "─" * width

    print(f"\n{BOLD}{CYAN}{sep}{RESET}")
    print(f"{BOLD}  BLASTP Summary{RESET}")
    print(f"{CYAN}{sep}{RESET}")
    print(f"  {'Total jobs':<20} {len(jobs)}")
    print(f"  {GREEN}{'Completed':<20}{RESET} {counts[JobStatus.SUCCESS]}")
    print(f"  {GREY}{'Skipped':<20}{RESET} {counts[JobStatus.SKIPPED]}")
    print(f"  {RED}{'Failed':<20}{RESET} {counts[JobStatus.FAILED]}")
    print(f"  {'Wall time':<20} {total_time:.1f}s")

    if failed:
        print(f"\n{RED}{BOLD}  Failed jobs:{RESET}")
        for j in failed:
            print(f"    {RED}✗{RESET}  {j.label}")
            if j.error:
                print(f"       {GREY}{j.error[:120]}{RESET}")

    print(f"{CYAN}{sep}{RESET}\n")


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run BLASTP for all protein / organism combinations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--workers", "-w",
        type=int, default=4,
        help="Number of parallel BLASTP processes",
    )
    p.add_argument(
        "--max-retry", "-r",
        type=int, default=2,
        help="Max retries per failed job (exponential back-off)",
    )
    p.add_argument(
        "--dry-run", "-n",
        action="store_true",
        help="Preview commands without running them or writing files",
    )
    p.add_argument(
        "--overwrite", "-o",
        action="store_true",
        help="Re-run BLAST even if output file already exists",
    )
    p.add_argument(
        "--genome-type", "-g",
        choices=["phased", "non_phased"],
        default="phased",
        help="Which genome assembly to query against",
    )
    p.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable DEBUG-level logging",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if args.verbose:
        log.setLevel(logging.DEBUG)

    log.info("Discovering jobs (genome: %s)…", args.genome_type)
    jobs = discover_jobs(genome_type=args.genome_type)
    log.info("Found %d job(s)", len(jobs))

    results = run_all(
        jobs,
        workers   = args.workers,
        dry_run   = args.dry_run,
        overwrite = args.overwrite,
        max_retry = args.max_retry,
    )

    print_summary(results)

    failed = sum(1 for j in results if j.status == JobStatus.FAILED)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())