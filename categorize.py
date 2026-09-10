import os
import json
import glob
import subprocess
from concurrent.futures import ThreadPoolExecutor
from rnabridge import Visualizer
from config import INPUT_DIR, OUTPUT_DIR, CIF_DIR, JSON_DIR, VARNA_CLI_PATH


def run_varna_task(v_json_path):
    """
    Executes the VARNA CLI tool to generate SVG visualizations from a VARNA JSON file.
    """
    cli_path = VARNA_CLI_PATH
    subprocess.run([cli_path, "varna-tz", str(v_json_path)])


def run_categorization():
    """
    Main categorization and visualization script:
    1. Reads analyzed helices and junctions from INPUT_DIR.
    2. Groups helices by the number of segments (stems).
    3. Generates PyMOL selections, .pml scripts, and .cif segments.
    4. Exports VARNA JSON data and triggers SVG generation.
    """
    print("--- STARTING VISUALIZATION GENERATION ---")
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR, exist_ok=True)

    files = glob.glob(os.path.join(str(INPUT_DIR), "*.json"))
    mapping_cache, total_saved, varna_tasks = {}, 0, []

    for index, filepath in enumerate(files, 1):
        base_name = os.path.splitext(os.path.basename(filepath))[0]
        if "_test" in base_name:
            continue

        cif_path = os.path.join(str(CIF_DIR), f"{base_name}.cif")
        raw_json_path = os.path.join(str(JSON_DIR), f"{base_name}.json")

        try:
            with open(filepath, "r", encoding="UTF-8") as f:
                data = json.load(f)

            if base_name not in mapping_cache and os.path.exists(raw_json_path):
                with open(raw_json_path, "r", encoding="UTF-8") as f_raw:
                    mapping_cache[base_name] = json.load(f_raw).get("bpseq_index")
                if len(mapping_cache) > 200:
                    mapping_cache.pop(next(iter(mapping_cache)))

            # 1. Process Helices
            for helix in data.get("helices", []):
                # Count segments robustly
                unique_stems = set()
                strands = helix.get("strands", {})
                for key in ["upstream", "downstream"]:
                    u_f = strands.get(key, {}).get("strand5p", {}).get("first")
                    if u_f and u_f.get("serial") is not None:
                        unique_stems.add(u_f.get("serial"))
                for comp in helix.get("components", []):
                    i_f = comp.get("internal_stem", {}).get("strand5p", {}).get("first")
                    if i_f and i_f.get("serial") is not None:
                        unique_stems.add(i_f.get("serial"))

                n = len(unique_stems)
                if n < 1:
                    continue

                h_id, save_dir = (
                    helix.get("h_id", "unknown"),
                    os.path.join(OUTPUT_DIR, f"{n}-segment-helis"),
                )
                os.makedirs(save_dir, exist_ok=True)
                out_base = f"{base_name}_h{h_id}_{n}-segment"

                svg_path = os.path.join(
                    save_dir, f"varna-tz-{out_base}_varna-clean.svg"
                )
                if not os.path.exists(svg_path):
                    full_sel = Visualizer.get_continuous_selection(
                        helix, is_junction=False
                    )
                    if full_sel != "none":
                        try:
                            from pymol import cmd

                            cmd.reinitialize()
                            cmd.load(cif_path, "target_obj")
                            cmd.select("target_sel", full_sel)
                            cmd.save(
                                os.path.join(save_dir, f"{out_base}.cif"), "target_sel"
                            )
                            Visualizer.generate_pml_script(
                                os.path.join(save_dir, f"{out_base}.pml"),
                                cif_path,
                                helix,
                                full_sel,
                                is_junction=False,
                            )
                            v_json = os.path.join(save_dir, f"{out_base}_varna.json")
                            Visualizer.export_varna_json(
                                helix,
                                v_json,
                                mapping=mapping_cache.get(base_name),
                                is_junction=False,
                            )
                            varna_tasks.append(v_json)
                            total_saved += 1
                            with open(
                                os.path.join(save_dir, f"{out_base}.json"), "w"
                            ) as out_f:
                                json.dump({"helices": [helix]}, out_f, indent=4)
                        except Exception as e:
                            print(f"   -> ERROR processing helix {out_base}: {e}")

            # 2. Process Junctions
            for junction in data.get("junctions", []):
                n, j_id = (
                    len(junction.get("location", [])),
                    junction.get("j_id", "unknown"),
                )
                save_dir = os.path.join(OUTPUT_DIR, f"{n}-way-junctions")
                os.makedirs(save_dir, exist_ok=True)
                out_base = f"{base_name}_j{j_id}_{n}way"

                svg_path = os.path.join(
                    save_dir, f"varna-tz-{out_base}_varna-clean.svg"
                )
                if not os.path.exists(svg_path):
                    full_sel = Visualizer.get_continuous_selection(
                        junction, is_junction=True
                    )
                    if full_sel != "none":
                        try:
                            from pymol import cmd

                            cmd.reinitialize()
                            cmd.load(cif_path, "target_obj")
                            cmd.select("target_sel", full_sel)
                            cmd.save(
                                os.path.join(save_dir, f"{out_base}.cif"), "target_sel"
                            )
                            Visualizer.generate_pml_script(
                                os.path.join(save_dir, f"{out_base}.pml"),
                                cif_path,
                                junction,
                                full_sel,
                                is_junction=True,
                            )
                            v_json = os.path.join(save_dir, f"{out_base}_varna.json")
                            Visualizer.export_varna_json(
                                junction,
                                v_json,
                                mapping=mapping_cache.get(base_name),
                                is_junction=True,
                            )
                            varna_tasks.append(v_json)
                            total_saved += 1
                            with open(
                                os.path.join(save_dir, f"{out_base}.json"), "w"
                            ) as out_f:
                                json.dump({"junctions": [junction]}, out_f, indent=4)
                        except Exception as e:
                            print(f"   -> ERROR processing junction {out_base}: {e}")

        except Exception as e:
            print(f"ERROR in file {filepath}: {e}")

        if index % 200 == 0:
            print(f"Progress {index}/{len(files)}...")

    if varna_tasks:
        print(f"Triggering VARNA image generation for {len(varna_tasks)} items...")
        max_workers = os.cpu_count() or 1
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            executor.map(run_varna_task, varna_tasks)

    print(f"\n--- FINISHED: Saved {total_saved} new elements ---")


if __name__ == "__main__":
    run_categorization()
