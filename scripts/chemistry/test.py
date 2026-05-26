"""Minimal pipeline runner wrapper.

Usage:
    python -m chemistry.test [PROTEIN] [OUTDIR]

PROTEIN may also be set via the environment variable CHEM_PROTEIN.
"""
import os
import sys
from chemistry.runner import run_pipeline
from data import PROTEINS_PATH, proteins, compounds


def _main(argv):
        
    protein      = "sard4"
    ref_organism = "tabaco"
    ooi          = "tucunaca"
    reaction     = "r1"
    substrate    = "tetrahydroharmol"
    cofactors    = proteins[protein]["cofactors"]

    boltz_path         = f"{PROTEINS_PATH}/{protein}/analysis/structural/boltz"
    structure_filename = "sample_0_predicted_structure.cif"

    if ooi:
        structure_path = (
            f"{boltz_path}/{ooi}/{ref_organism}"
            f"/{protein}_{ooi}_{ref_organism}_{reaction}/{structure_filename}"
        )
        out_dir = f"{boltz_path}/{ooi}/{ref_organism}/nac"
    else:
        structure_path = (
            f"{boltz_path}/{ref_organism}"
            f"/{protein}_{ref_organism}_{reaction}/{structure_filename}"
        )
        out_dir = f"{boltz_path}/{ref_organism}/nac"
        
    paths = run_pipeline(protein, outdir=out_dir)
    print('Wrote', len(paths), 'files')


if __name__ == '__main__':
        _main(sys.argv)
