# CLAUDE.md

Project overview and setup: see [README.md](README.md).

## Where results live: `pathways/`

`pathways/` holds every result the scripts produce. It is **git-ignored** (regenerable, too large for git) but present on disk — read it directly; don't assume a result is missing because git doesn't track it.

- `scripts/data.py` is the source of truth for paths. `CURRENT_PROJECT` selects the pathway (`harmine_biosynthesis`, `yuremamine_biosynthesis`, `benchmark`, …) and `PROJECT_PATH`, `PAPER_PATH`, `COMPOUNDS_PATH`, `PROTEINS_PATH`, `GENOME_PATH`, `TRANSCRIPTOME_PATH` derive from it. To find where a script writes or reads, check the path it builds from these constants.
- Compounds and proteins are keyed by the IDs in `inputs/<pathway>.yaml`.

Layout of `pathways/<pathway>/`:

```
compounds/<cpd>/
  structure/<cpd>_<pubchem_id>.{sdf,pdb,pdbqt,png}
  quantum/<cpd>_<pubchem_id>{,_optfreq,_sp}.inp, .out   # ORCA inputs/outputs (_optfreq = geometry actually computed)
  thermo/<cpd>_<pubchem_id>.json                         # parsed ORCA Gibbs data
proteins/<prot>/
  reference/<organism>/                                  # reference sequences/structures
  analysis/genomic/blast/<target_org>/<ref_org>/         # BLAST hits (*_blast_filtered.csv)
  analysis/phylo/, analysis/structural/{alignment,alphafold,boltz,vina}/
genome/<organism>/, transcriptome/<organism>/
paper/tables/   # CSV tables for the paper (generate_tables.py, <pathway>_gibbs.csv, …)
paper/figures/  # PNG figures
```

Bulky binary/raw formats (`.npz`, `.pkl`, `.cif`, `.a3m`, BLAST DBs, archives) are excluded from search via `.ignore`; open them by path when needed.

## Running scripts

Run from `scripts/` with its venv (`scripts/.venv`, Python 3.10). Subpackage scripts import `data` from `scripts/`, so run them as modules, e.g. `python -m thermo.reaction_gibbs`.
