from rnabridge import PDBDownloader
from config import CIF_DIR, JSON_DIR, RESULTS_DIR, OUTPUT_DIR


def sync_pdb_data():
    """
    Synchronizes local CIF files with the active RNA structures from RCSB PDB.
    1. Fetches all active RNA IDs.
    2. Removes local files for IDs that are no longer active or correctly named.
    3. Downloads missing CIF files.
    """
    print("--- STARTING PDB SYNCHRONIZATION ---")
    active_ids = PDBDownloader.get_active_rna_ids()
    if not active_ids:
        print("ERROR: Could not fetch active RNA IDs.")
        return

    # 1. Cleanup obsolete or malformed files
    print("Cleaning up obsolete local files...")
    for cif_file in CIF_DIR.glob("*.cif"):
        pdb_id = cif_file.stem.upper()
        if len(pdb_id) != 4 or pdb_id not in active_ids:
            print(f"Removing data for obsolete ID: {pdb_id}")
            cif_file.unlink()
            (JSON_DIR / f"{cif_file.stem}.json").unlink(missing_ok=True)
            (RESULTS_DIR / f"{cif_file.stem}.json").unlink(missing_ok=True)

            if OUTPUT_DIR.exists():
                for sub in OUTPUT_DIR.iterdir():
                    if sub.is_dir():
                        for f in sub.glob(f"*{cif_file.stem}*"):
                            f.unlink()

    # 2. Download missing CIF files
    print(f"Synchronizing {len(active_ids)} active structures...")
    downloaded_count = 0
    for pdb_id in active_ids:
        target = CIF_DIR / f"{pdb_id}.cif"
        if not target.exists():
            print(f"Downloading missing structure: {pdb_id}")
            if PDBDownloader.download_cif(pdb_id, str(target)):
                downloaded_count += 1

    print(f"--- SYNCHRONIZATION FINISHED. Downloaded {downloaded_count} new files. ---")


if __name__ == "__main__":
    sync_pdb_data()
