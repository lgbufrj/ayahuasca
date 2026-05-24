"""
download_sdf.py
─────────────────────────────────────────────────────────────────────────────
Downloads PubChem SDF files for every compound in data.compounds.
Tries 3D first; falls back to 2D if unavailable.
Skips compounds whose SDF already exists on disk.

Usage
    python download_sdf.py                    # all compounds
    python download_sdf.py --compound sam     # single compound
    python download_sdf.py --force            # re-download existing files
─────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import requests

import data
from logger import build_logger

log = build_logger("download_sdf")

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────

PUBCHEM_URL  = data.databases["pubchem"]["url"]
RETRY_DELAY  = 2      # seconds between retries
MAX_RETRIES  = 3
RECORD_TYPES = ("3d", "2d")   # preference order


# ──────────────────────────────────────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────────────────────────────────────

def _sdf_path(compound_key: str, pubchem_id: str | int) -> Path:
    return (
        Path(data.COMPOUNDS_PATH)
        / compound_key / "structure"
        / f"{compound_key}_{pubchem_id}.sdf"
    )


# ──────────────────────────────────────────────────────────────────────────────
# Download
# ──────────────────────────────────────────────────────────────────────────────

def _fetch_sdf(cid: str | int, record_type: str) -> bytes:
    """
    Fetch SDF bytes from PubChem. Retries on transient errors.
    Raises requests.HTTPError on a final failure.
    """
    url = (
        f"{PUBCHEM_URL}/rest/pug/compound/cid/{cid}/SDF"
        f"?record_type={record_type}"
    )

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            return response.content
        except requests.HTTPError as exc:
            # 404 means this record type genuinely doesn't exist — don't retry
            if exc.response is not None and exc.response.status_code == 404:
                raise
            if attempt < MAX_RETRIES:
                log.warning(
                    "attempt %d/%d failed for CID %s (%s): %s — retrying in %ds",
                    attempt, MAX_RETRIES, cid, record_type, exc, RETRY_DELAY,
                )
                time.sleep(RETRY_DELAY * attempt)
            else:
                raise
        except requests.RequestException as exc:
            if attempt < MAX_RETRIES:
                log.warning(
                    "attempt %d/%d failed for CID %s (%s): %s — retrying in %ds",
                    attempt, MAX_RETRIES, cid, record_type, exc, RETRY_DELAY,
                )
                time.sleep(RETRY_DELAY * attempt)
            else:
                raise

    raise RuntimeError("unreachable")   # pragma: no cover


def download_compound(
    compound_key: str,
    *,
    force: bool = False,
) -> bool:
    """
    Download the SDF for *compound_key*. Returns True on success.
    Tries 3D first, falls back to 2D.
    """
    cdata      = data.compounds[compound_key]
    pubchem_id = cdata["pubchem_id"]
    dest       = _sdf_path(compound_key, pubchem_id)

    if dest.exists() and not force:
        log.info("skip  %s (CID %s) — file exists", compound_key, pubchem_id)
        return True

    dest.parent.mkdir(parents=True, exist_ok=True)

    for record_type in RECORD_TYPES:
        try:
            log.info("fetch %s (CID %s) [%s]…", compound_key, pubchem_id, record_type)
            content = _fetch_sdf(pubchem_id, record_type)
            dest.write_bytes(content)
            log.info("saved %s → %s", compound_key, dest)
            return True
        except requests.HTTPError as exc:
            status = exc.response.status_code if exc.response is not None else "?"
            if status == 404:
                log.warning(
                    "%s (CID %s): no %s structure on PubChem, trying next…",
                    compound_key, pubchem_id, record_type,
                )
            else:
                log.error(
                    "%s (CID %s) [%s]: HTTP %s — %s",
                    compound_key, pubchem_id, record_type, status, exc,
                )
                return False
        except Exception as exc:
            log.error(
                "%s (CID %s) [%s]: %s",
                compound_key, pubchem_id, record_type, exc,
            )
            return False

    log.error("%s (CID %s): no structure available (tried %s)", compound_key, pubchem_id, ", ".join(RECORD_TYPES))
    return False


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Download PubChem SDF files for all compounds.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--compound", "-c",
        help="Single compound key to download (default: all)",
    )
    p.add_argument(
        "--force", "-f",
        action="store_true",
        help="Re-download even if the SDF already exists",
    )
    p.add_argument("--verbose", "-v", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    import logging
    args = parse_args(argv)

    if args.verbose:
        log.setLevel(logging.DEBUG)

    targets = [args.compound] if args.compound else list(data.compounds.keys())

    failed = []
    for key in targets:
        if key not in data.compounds:
            log.error("Unknown compound key: %s", key)
            return 1
        if not download_compound(key, force=args.force):
            failed.append(key)

    if failed:
        log.error("Failed downloads: %s", ", ".join(failed))
        return 1

    log.info("All done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())