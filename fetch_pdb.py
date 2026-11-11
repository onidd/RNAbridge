import sys
from rnabridge import PDBDownloader
from config import CIF_DIR, JSON_DIR, RESULTS_DIR, OUTPUT_DIR, DATASET

MAX_SAFE_DELETE = 100       
MIN_EXPECTED_ACTIVE_IDS = 5000  


def sync_pdb_data():
    """
    Synchronizes local CIF files with the active RNA structures from RCSB PDB.
    1. Fetches all active RNA IDs.
    2. Safely removes local files for IDs that are no longer active.
    3. Downloads missing CIF files.
    """
    print("--- STARTING PDB SYNCHRONIZATION ---")
    if DATASET == "rnasolo":
        print(
            "[INFO] RNAsolo dataset selected. "
            "Skipping RCSB PDB synchronization; no files were changed."
        )
        return
    active_ids = PDBDownloader.get_active_rna_ids()
    if not active_ids or len(active_ids) < MIN_EXPECTED_ACTIVE_IDS:
        print(
            f"ERROR: Could not fetch complete active RNA IDs list (got {len(active_ids) if active_ids else 0}). "
            f"Aborting synchronization!"
        )
        sys.exit(1)

    print("Checking for obsolete local files...")
    obsolete_files = []
    for cif_file in CIF_DIR.glob("*.cif"):
        pdb_id = cif_file.stem.upper()
        if len(pdb_id) != 4 or pdb_id not in active_ids:
            obsolete_files.append(cif_file)

    if len(obsolete_files) > MAX_SAFE_DELETE:
        print(
            f"WARNING [SAFETY GUARD]: {len(obsolete_files)} files marked as obsolete "
            f"(limit is {MAX_SAFE_DELETE}). Skipping deletion to prevent accidental data loss!"
        )
    else:
        for cif_file in obsolete_files:
            pdb_id = cif_file.stem.upper()
            print(f"Removing data for obsolete ID: {pdb_id}")
            cif_file.unlink()
            (JSON_DIR / f"{cif_file.stem}.json").unlink(missing_ok=True)
            (RESULTS_DIR / f"{cif_file.stem}.json").unlink(missing_ok=True)

            if OUTPUT_DIR.exists():
                for sub in OUTPUT_DIR.iterdir():
                    if sub.is_dir():
                        for f in sub.glob(f"*{cif_file.stem}*"):
                            f.unlink()

    print(f"Synchronizing {len(active_ids)} active structures...")
    downloaded_count = 0
    for pdb_id in sorted(active_ids):
        target = CIF_DIR / f"{pdb_id}.cif"
        if not target.exists():
            print(f"Downloading missing structure: {pdb_id}")
            if PDBDownloader.download_cif(pdb_id, str(target)):
                downloaded_count += 1

    print(f"--- SYNCHRONIZATION FINISHED. Downloaded {downloaded_count} new files. ---")


if __name__ == "__main__":
    sync_pdb_data()
