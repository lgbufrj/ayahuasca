from data import PROTEINS_PATH, SCRIPTS_PATH
import re
import shutil
from pathlib import Path

results_dir = Path(f"{SCRIPTS_PATH}/structural/boltz/ayahuasca_paper/ayahuasca_boltz_results")


def parse_folder_name(name: str) -> tuple[dict, bool] | tuple[None, None]:
    """
    Try to match the folder name against both templates.
    Returns (params_dict, is_ooi) or (None, None) if no match.

    Anchors on rxn_id (r<digits> at the end) and pubchem_id (the segment
    immediately before rxn_id). Everything before those two segments is
    either [ref_organism, protein] (REF) or [ooi, ref_organism, protein] (OOI).
    """
    suffix = name[len("boltz_results_"):]
    parts = suffix.split("_")

    # rxn_id is always the last segment and is one of: rr, r1, r2, r3
    if not re.match(r"^r[r123]$", parts[-1]):
        print(f"  [SKIP] Last segment doesn't look like rxn_id (rr/r1/r2/r3): {name}")
        return None, None

    rxn_id    = parts[-1]
    pubchem_id = parts[-2]
    leading   = parts[:-2]  # everything before pubchem_id

    n_leading = len(leading)
    if n_leading == 2:       # ref_organism + protein  -> REF template
        params = dict(ref_organism=leading[0], protein=leading[1],
                      pubchem_id=pubchem_id, rxn_id=rxn_id)
        return params, False
    elif n_leading == 3:     # ooi + ref_organism + protein -> OOI template
        params = dict(ooi=leading[0], ref_organism=leading[1], protein=leading[2],
                      pubchem_id=pubchem_id, rxn_id=rxn_id)
        return params, True
    else:
        print(f"  [SKIP] Unexpected number of leading segments ({n_leading}) in: {name}")
        return None, None


def build_output_dir(params: dict, is_ooi: bool) -> Path:
    if is_ooi:
        return Path(
            f"{PROTEINS_PATH}/{params['protein']}/analysis/structural/boltz"
            f"/{params['ooi']}/{params['ref_organism']}/{params['rxn_id']}"
        )
    else:
        return Path(
            f"{PROTEINS_PATH}/{params['protein']}/analysis/structural/boltz"
            f"/{params['ref_organism']}/{params['rxn_id']}"
        )


def copy_folder_contents(src: Path, dst: Path) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    for item in src.iterdir():
        dest_item = dst / item.name
        if item.is_dir():
            shutil.copytree(item, dest_item, dirs_exist_ok=True)
        else:
            shutil.copy2(item, dest_item)


def main():
    if not results_dir.exists():
        print(f"Results directory not found: {results_dir}")
        return

    folders = [f for f in results_dir.iterdir() if f.is_dir()]
    if not folders:
        print("No subdirectories found.")
        return

    print(f"Found {len(folders)} folder(s) in {results_dir}\n")

    for folder in sorted(folders):
        name = folder.name
        print(f"Processing: {name}")

        params, is_ooi = parse_folder_name(name)
        if params is None:
            print(f"  [SKIP] Could not parse folder name.\n")
            continue

        print(f"  Params: {params}  |  ooi={'yes' if is_ooi else 'no'}")

        output_dir = build_output_dir(params, is_ooi)
        print(f"  -> Output: {output_dir}")

        copy_folder_contents(folder, output_dir)
        print(f"  [OK] Copied.\n")


if __name__ == "__main__":
    main()