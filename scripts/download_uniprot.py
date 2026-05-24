"""
download_structures.py
─────────────────────────────────────────────────────────────────────────────
Downloads per-organism reference structures for every protein in data.proteins:

  1. FASTA         — from UniProt / UniParc
  2. PDB           — from UniProt (AlphaFold or experimental, whatever UniProt
                     serves for that accession)
  3. X-ray PDB     — best-resolution X-ray crystal structure from RCSB,
                     queried by UniProt ID  (skipped if none available)

Output layout
    <PROTEINS_PATH>/<ptn>/<reference>/<species>/structure/
        <ptn>_<uniprot_id>.fasta
        <ptn>_<uniprot_id>.pdb
        <ptn>_<uniprot_id>_xray.pdb      ← best X-ray structure, if any

Usage
    python download_structures.py                        # all proteins
    python download_structures.py --protein asmt         # single protein
    python download_structures.py --force                # re-download all
    python download_structures.py --skip-xray            # FASTA + PDB only
─────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Optional

import requests

import data
from logger import build_logger

log = build_logger("download_structures")

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────

UNIPROT_URL  = data.databases["uniprot"]["url"]        # e.g. https://rest.uniprot.org
RCSB_SEARCH  = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_FETCH   = "https://files.rcsb.org/download/{pdb_id}.pdb"

RETRY_DELAY  = 2
MAX_RETRIES  = 3


# ──────────────────────────────────────────────────────────────────────────────
# Shared HTTP helper
# ──────────────────────────────────────────────────────────────────────────────

def _get(url: str, *, label: str, **kwargs) -> requests.Response:
    """GET with retry and exponential back-off. Raises on final failure."""
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.get(url, timeout=30, **kwargs)
            resp.raise_for_status()
            return resp
        except requests.HTTPError as exc:
            status = exc.response.status_code if exc.response is not None else "?"
            if status == 404 or attempt == MAX_RETRIES:
                raise
            log.warning(
                "attempt %d/%d  %s: HTTP %s — retrying in %ds",
                attempt, MAX_RETRIES, label, status, RETRY_DELAY * attempt,
            )
            time.sleep(RETRY_DELAY * attempt)
        except requests.RequestException as exc:
            if attempt == MAX_RETRIES:
                raise
            log.warning(
                "attempt %d/%d  %s: %s — retrying in %ds",
                attempt, MAX_RETRIES, label, exc, RETRY_DELAY * attempt,
            )
            time.sleep(RETRY_DELAY * attempt)

    raise RuntimeError("unreachable")  # pragma: no cover


# ──────────────────────────────────────────────────────────────────────────────
# Path helpers
# ──────────────────────────────────────────────────────────────────────────────

def _struct_dir(ptn_name: str, species: str) -> Path:
    return (
        Path(data.PROTEINS_PATH)
        / ptn_name / "reference" / species / "structure"
    )


def _fasta_path(ptn_name: str, species: str, uniprot_id: str) -> Path:
    return _struct_dir(ptn_name, species) / f"{ptn_name}_{uniprot_id}.fasta"


def _pdb_path(ptn_name: str, species: str, uniprot_id: str) -> Path:
    return _struct_dir(ptn_name, species) / f"{ptn_name}_{uniprot_id}.pdb"


def _xray_path(ptn_name: str, species: str, uniprot_id: str) -> Path:
    return _struct_dir(ptn_name, species) / f"{ptn_name}_{uniprot_id}_xray.pdb"


# ──────────────────────────────────────────────────────────────────────────────
# 1. FASTA
# ──────────────────────────────────────────────────────────────────────────────

def download_fasta(
    ptn_name: str,
    species: str,
    uniprot_id: str,
    *,
    force: bool = False,
) -> bool:
    dest = _fasta_path(ptn_name, species, uniprot_id)

    if dest.exists() and not force:
        log.info("skip  FASTA  %s/%s (%s)", ptn_name, species, uniprot_id)
        return True

    db  = "uniparc" if uniprot_id.startswith("UPI") else "uniprot"
    url = f"{UNIPROT_URL}/{db}/{uniprot_id}.fasta"

    try:
        resp = _get(url, label=f"FASTA {uniprot_id}")
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(resp.text)
        log.info("saved FASTA  %s/%s → %s", ptn_name, species, dest.name)
        return True
    except Exception as exc:
        log.error("FASTA  %s/%s (%s): %s", ptn_name, species, uniprot_id, exc)
        return False


# ──────────────────────────────────────────────────────────────────────────────
# 2. UniProt PDB  (AlphaFold or experimental)
# ──────────────────────────────────────────────────────────────────────────────

def download_uniprot_pdb(
    ptn_name: str,
    species: str,
    uniprot_id: str,
    *,
    force: bool = False,
) -> bool:
    dest = _pdb_path(ptn_name, species, uniprot_id)

    if dest.exists() and not force:
        log.info("skip  PDB    %s/%s (%s)", ptn_name, species, uniprot_id)
        return True

    # UniParc accessions don't have PDB structures
    if uniprot_id.startswith("UPI"):
        log.warning("PDB    %s/%s: UniParc IDs have no PDB — skipping", ptn_name, species)
        return True

    url = f"{UNIPROT_URL}/uniprot/{uniprot_id}.pdb"

    try:
        resp = _get(url, label=f"PDB {uniprot_id}")
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(resp.content)
        log.info("saved PDB    %s/%s → %s", ptn_name, species, dest.name)
        return True
    except requests.HTTPError as exc:
        status = exc.response.status_code if exc.response is not None else "?"
        if status == 404:
            log.warning("PDB    %s/%s (%s): not available on UniProt", ptn_name, species, uniprot_id)
            return True   # not a hard failure — structure simply doesn't exist
        log.error("PDB    %s/%s (%s): HTTP %s", ptn_name, species, uniprot_id, status)
        return False
    except Exception as exc:
        log.error("PDB    %s/%s (%s): %s", ptn_name, species, uniprot_id, exc)
        return False


# ──────────────────────────────────────────────────────────────────────────────
# 3. Best X-ray structure from RCSB
# ──────────────────────────────────────────────────────────────────────────────

def _find_best_xray(uniprot_id: str) -> Optional[str]:
    """
    Query RCSB for X-ray crystal structures mapped to *uniprot_id*.
    Returns the PDB ID with the best (lowest) resolution, or None.
    """
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator":  "exact_match",
                        "value":     uniprot_id,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator":  "exact_match",
                        "value":     "X-RAY DIFFRACTION",
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "sort": [
                {
                    "sort_by":    "rcsb_entry_info.resolution_combined",
                    "direction":  "asc",
                }
            ],
            "paginate": {"start": 0, "rows": 1},
        },
    }

    try:
        resp = requests.post(RCSB_SEARCH, json=query, timeout=30)
        resp.raise_for_status()
        results = resp.json().get("result_set", [])
        if results:
            return results[0]["identifier"]
        return None
    except Exception as exc:
        log.warning("RCSB search for %s: %s", uniprot_id, exc)
        return None


def download_xray_pdb(
    ptn_name: str,
    species: str,
    uniprot_id: str,
    *,
    force: bool = False,
) -> bool:
    dest = _xray_path(ptn_name, species, uniprot_id)

    if dest.exists() and not force:
        log.info("skip  X-ray  %s/%s (%s)", ptn_name, species, uniprot_id)
        return True

    if uniprot_id.startswith("UPI"):
        log.warning("X-ray  %s/%s: UniParc IDs have no RCSB mapping — skipping", ptn_name, species)
        return True

    log.info("search X-ray %s/%s (%s)…", ptn_name, species, uniprot_id)
    pdb_id = _find_best_xray(uniprot_id)

    if pdb_id is None:
        log.info("X-ray  %s/%s (%s): no X-ray structure found on RCSB", ptn_name, species, uniprot_id)
        return True   # not a failure — structure simply doesn't exist

    url = RCSB_FETCH.format(pdb_id=pdb_id)
    try:
        resp = _get(url, label=f"X-ray {pdb_id}")
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(resp.content)
        log.info(
            "saved X-ray  %s/%s → %s  (RCSB: %s)",
            ptn_name, species, dest.name, pdb_id,
        )
        return True
    except Exception as exc:
        log.error("X-ray  %s/%s (%s / %s): %s", ptn_name, species, uniprot_id, pdb_id, exc)
        return False


# ──────────────────────────────────────────────────────────────────────────────
# Orchestration
# ──────────────────────────────────────────────────────────────────────────────

def download_protein(
    ptn_name: str,
    *,
    force:     bool = False,
    skip_xray: bool = False,
) -> bool:
    ptn_data = data.proteins[ptn_name]
    ok = True

    for species, org_data in ptn_data.get("organisms", {}).items():
        if not org_data or "uniprot_id" not in org_data:
            log.warning("%s/%s: no uniprot_id — skipping", ptn_name, species)
            continue

        uniprot_id = org_data["uniprot_id"]
        log.info("── %s / %s  (%s)", ptn_name, species, uniprot_id)

        ok &= download_fasta(ptn_name, species, uniprot_id, force=force)
        ok &= download_uniprot_pdb(ptn_name, species, uniprot_id, force=force)

        if not skip_xray:
            ok &= download_xray_pdb(ptn_name, species, uniprot_id, force=force)

    return ok


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Download FASTA, PDB, and X-ray structures for all proteins.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--protein", "-p",
        help="Single protein key to process (default: all)",
    )
    p.add_argument(
        "--force", "-f",
        action="store_true",
        help="Re-download files that already exist",
    )
    p.add_argument(
        "--skip-xray",
        action="store_true",
        help="Skip the RCSB X-ray structure search",
    )
    p.add_argument("--verbose", "-v", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    import logging
    args = parse_args(argv)

    if args.verbose:
        log.setLevel(logging.DEBUG)

    targets = [args.protein] if args.protein else list(data.proteins.keys())

    failed = []
    for ptn in targets:
        if ptn not in data.proteins:
            log.error("Unknown protein key: %s", ptn)
            return 1
        if not download_protein(ptn, force=args.force, skip_xray=args.skip_xray):
            failed.append(ptn)

    if failed:
        log.error("Completed with errors: %s", ", ".join(failed))
        return 1

    log.info("All done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())