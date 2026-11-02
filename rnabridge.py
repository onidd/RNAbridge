import time
import itertools
from typing import Dict, List, Set, Tuple, Optional, Any
from collections import Counter
import os
import math
import tempfile
import subprocess
import pymol
import json
from pymol import cmd

import urllib.request
import requests


class PDBDownloader:
    """Handles synchronization and downloading of RNA structures from RCSB PDB."""

    @staticmethod
    def get_active_rna_ids() -> Set[str]:
        """Fetches the list of all active RNA entry IDs from RCSB PDB."""
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
        all_ids, start, rows = set(), 0, 1000
        expected_total = None

        while True:
            query = {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_entity_polymer_type",
                        "operator": "exact_match",
                        "value": "RNA",
                    },
                },
                "return_type": "entry",
                "request_options": {"paginate": {"start": start, "rows": rows}},
            }

            success = False
            for attempt in range(1, 4):
                try:
                    r = requests.post(url, json=query, timeout=30)
                    r.raise_for_status()
                    data = r.json()
                    if expected_total is None:
                        expected_total = data.get("total_count", 0)
                    results = data.get("result_set", [])
                    for res in results:
                        all_ids.add(res["identifier"].upper())
                    success = True
                    break
                except Exception as e:
                    print(f"[WARNING] RCSB API error at start={start} (attempt {attempt}/3): {e}")
                    time.sleep(5)

            if not success:
                print(f"[CRITICAL] Failed to fetch full list from RCSB PDB! Stopped at {len(all_ids)}/{expected_total}.")
                return set()

            if not results or (expected_total and len(all_ids) >= expected_total):
                break
            start += rows

        if not expected_total or len(all_ids) < expected_total:
            print(f"[CRITICAL] Incomplete PDB ID list fetched ({len(all_ids)}/{expected_total}).")
            return set()

        return all_ids

    @staticmethod
    def download_cif(pdb_id: str, target_path: str) -> bool:
        """Downloads a CIF file for a given PDB ID."""
        for attempt in range(1, 4):
            try:
                r = requests.get(
                    f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", timeout=20
                )
                if r.status_code == 200 and len(r.content) > 100:
                    with open(target_path, "wb") as f:
                        f.write(r.content)
                    return True
            except Exception:
                time.sleep(2)
        return False


class MetadataProvider:
    """Fetches biological and experimental metadata from RCSB PDB API."""

    @staticmethod
    def fetch_metadata(pdb_id: str) -> Tuple[str, str, Optional[float], str]:
        """Returns (title, organism, resolution, method) for a given PDB ID."""
        clean_id = pdb_id.upper()[-4:]
        try:
            url = f"https://data.rcsb.org/rest/v1/core/entry/{clean_id}"
            with urllib.request.urlopen(url, timeout=5) as resp:
                data = json.loads(resp.read().decode())
                title = data.get("struct", {}).get("title", "Unknown Molecule")
                exptl = data.get("exptl", [])
                method = exptl[0].get("method", "Unknown") if exptl else "Unknown"
                res_val = data.get("rcsb_entry_info", {}).get("resolution_combined", [])
                res = res_val[0] if res_val else None

            url_ent = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{clean_id}/1"
            with urllib.request.urlopen(url_ent, timeout=5) as resp:
                ent_data = json.loads(resp.read().decode())
                src = (
                    ent_data.get("rcsb_entity_source_organism", [])
                    or ent_data.get("entity_src_gen", [])
                    or ent_data.get("entity_src_nat", [])
                )
                org = (
                    src[0].get("scientific_name") or "Unknown Organism"
                    if src
                    else "Unknown Organism"
                )
            return title[:200], org[:100], res, method
        except Exception:
            return "Unknown Molecule", "Unknown Organism", None, "Unknown Method"


class Utils:
    """General utility functions for RNA data processing."""

    @staticmethod
    def to_int(val: Any) -> Optional[int]:
        """Converts various input types (dict with serial, string, int) to an integer."""
        if isinstance(val, dict):
            return val.get("serial")
        try:
            return int(val)
        except Exception:
            return None

    @staticmethod
    def get_full_sequence_from_data(data_obj: Dict) -> str:
        """Reconstructs the full RNA sequence from a complex nested data object (helices/junctions)."""
        chains = {}

        def process_strand(strand):
            seq, f = strand.get("sequence", ""), strand.get("first", {})
            chain, start = f.get("chain"), f.get("serial")
            if seq and chain and start is not None:
                if chain not in chains:
                    chains[chain] = {}
                for i, char in enumerate(seq):
                    chains[chain][start + i] = char

        # Traverse all possible location fields
        for loc in data_obj.get("location", []):
            process_strand(loc)
        for stem in data_obj.get("context", {}).values():
            if isinstance(stem, dict):
                process_strand(stem.get("strand5p", {}))
                process_strand(stem.get("strand3p", {}))
        for ext in data_obj.get("extended_helices", []):
            strands = ext.get("strands", {})
            for s in ["upstream", "downstream"]:
                process_strand(strands.get(s, {}).get("strand5p", {}))
                process_strand(strands.get(s, {}).get("strand3p", {}))
            for comp in ext.get("components", []):
                for loc in comp.get("location", []):
                    process_strand(loc)
                process_strand(comp.get("internal_stem", {}).get("strand5p", {}))
                process_strand(comp.get("internal_stem", {}).get("strand3p", {}))
        for s in ["upstream", "downstream"]:
            process_strand(data_obj.get("strands", {}).get(s, {}).get("strand5p", {}))
            process_strand(data_obj.get("strands", {}).get(s, {}).get("strand3p", {}))
        for c in data_obj.get("components", []):
            for loc in c.get("location", []):
                process_strand(loc)
            process_strand(c.get("internal_stem", {}).get("strand5p", {}))
            process_strand(c.get("internal_stem", {}).get("strand3p", {}))

        parts = []
        for cid in sorted(chains.keys()):
            nts = chains[cid]
            sorted_s = sorted(nts.keys())
            if not sorted_s:
                continue
            curr, prev = "", None
            for s in sorted_s:
                if prev is not None and s != prev + 1:
                    parts.append(curr)
                    curr = nts[s]
                else:
                    curr += nts[s]
                prev = s
            parts.append(curr)
        return "-".join(parts)

    @staticmethod
    def get_closing_indices(motif: Dict) -> Set[int]:
        """Returns a set of serial indices for the closing base pairs of a motif."""
        indices = set()
        for s in motif.get("strands", []):
            f, l = Utils.to_int(s.get("first")), Utils.to_int(s.get("last"))
            if f is not None:
                indices.add(f)
            if l is not None:
                indices.add(l)
        return indices

    @staticmethod
    def enrich_motif_data(motif: Dict, full_data: Any) -> None:
        """Maps raw serial indices to detailed PDB auth numbering using the provided index."""
        if not isinstance(full_data, dict) or "bpseq_index" not in full_data:
            return
        mapping = full_data["bpseq_index"]
        strands_to_enrich = []
        if "strands" in motif:
            strands_to_enrich.extend(motif["strands"])
        elif "location" in motif and "strands" in motif["location"]:
            strands_to_enrich.extend(motif["location"]["strands"])
        if "strand5p" in motif:
            strands_to_enrich.append(motif["strand5p"])
        if "strand3p" in motif:
            strands_to_enrich.append(motif["strand3p"])

        def clean_icode(val: Any) -> Optional[str]:
            if val in [None, ".", "?", " "]:
                return None
            if isinstance(val, str) and not val.strip():
                return None
            return val

        for s in strands_to_enrich:
            for key_name in ["first", "last"]:
                raw = s.get(key_name)
                if raw is not None and not isinstance(raw, dict):
                    key = str(raw)
                    if key in mapping:
                        auth = mapping[key]["auth"]
                        s[key_name] = {
                            "chain": auth["chain"],
                            "number": auth["number"],
                            "icode": clean_icode(auth.get("icode")),
                            "name": auth["name"],
                            "serial": int(key),
                        }
            if "all_serials" in s:
                del s["all_serials"]

        interactions_details = (
            motif.get("modules", {}).get("interactions", {}).get("details", [])
        )
        for inter in interactions_details:
            for key in ["nt1", "nt2"]:
                raw_nt = inter.get(key)
                if raw_nt is not None and not isinstance(raw_nt, dict):
                    key_str = str(raw_nt)
                    if key_str in mapping:
                        auth = mapping[key_str]["auth"]
                        inter[key] = {
                            "chain": auth["chain"],
                            "number": auth["number"],
                            "icode": clean_icode(auth.get("icode")),
                            "name": auth["name"],
                            "serial": int(key_str),
                        }

    @staticmethod
    def get_stem_str(stem: Dict) -> str:
        """Returns a string representation of a stem's range for logging/debugging."""
        if not stem:
            return "None"
        if "strand5p" in stem and "strand3p" in stem:
            s5, s3 = stem["strand5p"], stem["strand3p"]
            f5, l5, f3, l3 = (
                s5.get("first", {}),
                s5.get("last", {}),
                s3.get("first", {}),
                s3.get("last", {}),
            )
            return f"[{f5.get('chain')}{f5.get('number')}-{l5.get('number')}, {f3.get('chain')}{f3.get('number')}-{l3.get('number')}]"
        elif "first" in stem:
            f, l = stem.get("first", {}), stem.get("last")
            if isinstance(f, dict):
                return f"[{f.get('chain')}{f.get('number')}-{l.get('number')}]"
            return f"[{f}-{l}]"
        return "None"

    @staticmethod
    def get_stem_nt_count(stem: Dict) -> int:
        """Calculates the total number of nucleotides in a given stem."""
        if not stem:
            return 0
        if "strand5p" in stem and "strand3p" in stem:
            l1 = len(stem["strand5p"].get("sequence", "")) or (
                abs(
                    Utils.to_int(stem["strand5p"].get("last", 0))
                    - Utils.to_int(stem["strand5p"].get("first", 0))
                )
                + 1
            )
            l2 = len(stem["strand3p"].get("sequence", "")) or (
                abs(
                    Utils.to_int(stem["strand3p"].get("last", 0))
                    - Utils.to_int(stem["strand3p"].get("first", 0))
                )
                + 1
            )
            return l1 + l2
        elif "first" in stem and "last" in stem:
            return (
                abs(Utils.to_int(stem["last"]) - Utils.to_int(stem["first"])) + 1
            ) * 2
        return 0

    @staticmethod
    def get_helix_unique_nt(motifs_in_helix: List[Dict]) -> int:
        """Counts unique nucleotides across all motifs and stems in a super-helix."""
        unique_nt = set()

        def add_range(f, l, chain):
            if f is not None and l is not None:
                for i in range(min(f, l), max(f, l) + 1):
                    unique_nt.add((chain, i))

        for m in motifs_in_helix:
            for s in m.get("location", {}).get("strands", []):
                add_range(
                    Utils.to_int(s.get("first")),
                    Utils.to_int(s.get("last")),
                    s.get("chain"),
                )
            context = m.get("location", {}).get("context", {})
            for key, stem in context.items():
                if "stem" in key or key in ["upstream", "downstream"]:
                    for skey in ["strand5p", "strand3p"]:
                        s = stem.get(skey, {})
                        add_range(
                            Utils.to_int(s.get("first")),
                            Utils.to_int(s.get("last")),
                            s.get("chain"),
                        )
        return len(unique_nt)


class Core:
    """Core RNA structural logic and classification."""
    @staticmethod
    def _has_pseudoknot_char(structure: str) -> bool:
        return False

    @staticmethod
    def loop_nt_count(structure: str) -> int:
        """Counts loop nucleotides (dots AND pseudoknot characters, i.e., anything except '(' and ')')."""
        return sum(1 for c in structure if c not in "()")

    @staticmethod
    def classify_motif(loop_strands: List[Dict]) -> Optional[str]:
        """Classifies an RNA loop into HAIRPIN, BULGE, INTERNAL_LOOP, or N-WAY_JUNCTION."""
        if len(loop_strands) == 1:
            s = loop_strands[0]
            if Core.loop_nt_count(s.get("structure", "")) > 0:
                return "HAIRPIN"
            return None

        if len(loop_strands) == 2:
            has_loop_nts = [
                (Core.loop_nt_count(s.get("structure", "")) > 0) for s in loop_strands
            ]
            if has_loop_nts[0] and has_loop_nts[1]:
                return "INTERNAL_LOOP"
            elif has_loop_nts[0] ^ has_loop_nts[1]:
                return "BULGE"
            return None

        if len(loop_strands) >= 3:
            return f"{len(loop_strands)}_WAY_JUNCTION"
        return None

    @staticmethod
    def build_pairs_index(data: Dict) -> Dict:
        """Creates a fast lookup index for base pairs based on serial indices."""
        mapping = data.get("bpseq_index", {})
        reverse_map = {}
        for s, l in mapping.items():
            if l and l.get("label"):
                label = l["label"]
                reverse_map[(label.get("chain"), label.get("number"))] = int(s)

        index = {}
        for bp in data.get("base_pairs", []):
            try:
                nt1, nt2 = bp.get("nt1"), bp.get("nt2")
                if not nt1 or not nt2:
                    continue
                l1, l2 = nt1.get("label"), nt2.get("label")
                if not l1 or not l2:
                    continue

                s1, s2 = (
                    reverse_map.get((l1.get("chain"), l1.get("number"))),
                    reverse_map.get((l2.get("chain"), l2.get("number"))),
                )
                if s1 is None or s2 is None:
                    continue
                key = tuple(sorted((s1, s2)))
                index[key] = {
                    "lw": bp.get("lw") or bp.get("LW") or "???",
                    "seq": f"{l1.get('name')}-{l2.get('name')}",
                }
            except Exception:
                continue
        return index

    @staticmethod
    def find_interactions(loop_strands: List[Dict], pairs_index: Dict) -> List[Dict]:
        """Finds non-canonical or tertiary interactions within a given set of loop strands."""
        loop_indices = set()
        for s in loop_strands:
            start, end = Utils.to_int(s.get("first")), Utils.to_int(s.get("last"))
            if start is not None and end is not None:
                loop_indices.update(range(min(start, end), max(start, end) + 1))
        interactions = []
        for i1, i2 in itertools.combinations(sorted(list(loop_indices)), 2):
            key = tuple(sorted((i1, i2)))
            if key in pairs_index:
                info = pairs_index[key]
                if info["lw"] == "cWW" and info["seq"] in [
                    "G-C",
                    "C-G",
                    "A-U",
                    "U-A",
                    "G-U",
                    "U-G",
                ]:
                    continue
                interactions.append(
                    {"nt1": i1, "nt2": i2, "seq": info["seq"], "lw": info["lw"]}
                )
        return interactions


class Analysis:
    """Higher-level categorization of RNA interactions."""

    @staticmethod
    def categorize_interaction(motif: Dict, interactions: List[Dict]) -> Dict[str, int]:
        """Categorizes interactions into TRIPLET, EXTERNAL, or INTERNAL for internal loops."""
        if not interactions:
            return {}
        closing_idxs, counts, nt_to_interactions, valid_interactions = (
            Utils.get_closing_indices(motif),
            {},
            {},
            [],
        )
        for i in interactions:
            n1, n2 = i["nt1"], i["nt2"]
            if Utils.to_int(n1) in closing_idxs and Utils.to_int(n2) in closing_idxs:
                continue
            valid_interactions.append(i)
            if Utils.to_int(n1) not in closing_idxs:
                nt_to_interactions.setdefault(Utils.to_int(n1), []).append(i)
            if Utils.to_int(n2) not in closing_idxs:
                nt_to_interactions.setdefault(Utils.to_int(n2), []).append(i)
        used_ids = set()
        for nt, bonds in nt_to_interactions.items():
            if len(bonds) >= 2:
                counts["TRIPLET"] = counts.get("TRIPLET", 0) + 1
                for b in bonds:
                    used_ids.add(id(b))
        for i in valid_interactions:
            if id(i) in used_ids:
                continue
            if (
                Utils.to_int(i["nt1"]) in closing_idxs
                or Utils.to_int(i["nt2"]) in closing_idxs
            ):
                counts["EXTERNAL"] = counts.get("EXTERNAL", 0) + 1
            else:
                counts["INTERNAL"] = counts.get("INTERNAL", 0) + 1
        return counts

    @staticmethod
    def categorize_junction_interaction(
        motif: Dict, interactions: List[Dict]
    ) -> Dict[str, int]:
        """Categorizes interactions for junction loops (shares logic with generic internal loops)."""
        return Analysis.categorize_interaction(motif, interactions)

    @staticmethod
    def categorize_bulge_interaction(
        motif: Dict, interactions: List[Dict]
    ) -> Dict[str, int]:
        """Categorizes interactions into TRIPLET, HELIX_ANCHOR, or INTRA_LOOP for bulges and hairpins."""
        if not interactions:
            return {}
        closing_idxs, counts = Utils.get_closing_indices(motif), {}
        valid_interactions = [
            i
            for i in interactions
            if not (
                Utils.to_int(i["nt1"]) in closing_idxs
                and Utils.to_int(i["nt2"]) in closing_idxs
            )
        ]
        nt_to_interactions = {}
        for i in valid_interactions:
            if Utils.to_int(i["nt1"]) not in closing_idxs:
                nt_to_interactions.setdefault(Utils.to_int(i["nt1"]), []).append(i)
            if Utils.to_int(i["nt2"]) not in closing_idxs:
                nt_to_interactions.setdefault(Utils.to_int(i["nt2"]), []).append(i)
        used_ids = set()
        for nt, bonds in nt_to_interactions.items():
            if len(bonds) >= 2:
                counts["TRIPLET"] = counts.get("TRIPLET", 0) + 1
                for b in bonds:
                    used_ids.add(id(b))
        for i in valid_interactions:
            if id(i) in used_ids:
                continue
            if (
                Utils.to_int(i["nt1"]) in closing_idxs
                or Utils.to_int(i["nt2"]) in closing_idxs
            ):
                counts["HELIX_ANCHOR"] = counts.get("HELIX_ANCHOR", 0) + 1
            else:
                counts["INTRA_LOOP"] = counts.get("INTRA_LOOP", 0) + 1
        return counts

    @staticmethod
    def categorize_hairpin_interaction(motif: dict, interactions: list) -> dict:
        """Categorizes interactions for hairpins (shares logic with bulges)."""
        return Analysis.categorize_bulge_interaction(motif, interactions)


class Stacking:
    """Analysis of coaxial stacking and stacking paths."""

    @staticmethod
    def build_stacking_index(stacking_list: list, full_data: dict = None) -> set:
        """Creates an index of stacked nucleotide pairs for fast graph traversal."""
        idx = set()
        reverse_map = {}
        if full_data and "bpseq_index" in full_data:
            mapping = full_data["bpseq_index"]
            for s_idx, l_data in mapping.items():
                if l_data and l_data.get("label"):
                    l = l_data["label"]
                    reverse_map[(l.get("chain"), l.get("number"))] = int(s_idx)

        for s in stacking_list:
            try:
                nt1 = s.get("nt1", {})
                nt2 = s.get("nt2", {})
                l1, l2 = nt1.get("label", {}), nt2.get("label", {})
                key1 = (l1.get("chain"), l1.get("number"))
                key2 = (l2.get("chain"), l2.get("number"))

                if reverse_map and key1 in reverse_map and key2 in reverse_map:
                    n1, n2 = reverse_map[key1], reverse_map[key2]
                else:
                    n1, n2 = l1.get("number"), l2.get("number")

                if n1 is not None and n2 is not None:
                    idx.add(frozenset([int(n1), int(n2)]))
            except Exception:
                continue
        return idx

    @staticmethod
    def _get_stacking_path_info(
        starts: List[int],
        targets: List[int],
        valid_nodes: Set[int],
        stacking_index: Set[frozenset],
    ) -> Tuple[bool, int, List[int]]:
        """Breadth-first search to find the shortest stacking path between sets of nucleotides."""
        if not starts or not targets:
            return False, 0, []
        graph = {}
        for edge in stacking_index:
            u, v = list(edge)
            if u in valid_nodes and v in valid_nodes:
                graph.setdefault(u, set()).add(v)
                graph.setdefault(v, set()).add(u)

        visited = {s: (None, 0) for s in starts}  # node: (parent, depth)
        queue = list(starts)
        target_set = set(targets)
        found_target = None

        while queue:
            curr = queue.pop(0)
            if curr in target_set:
                found_target = curr
                break
            for neighbor in graph.get(curr, []):
                if neighbor not in visited:
                    visited[neighbor] = (curr, visited[curr][1] + 1)
                    queue.append(neighbor)

        if found_target is None:
            return False, 0, []

        path_nodes = []
        curr = found_target
        while curr is not None:
            path_nodes.append(curr)
            curr = visited[curr][0]

        full_path = path_nodes[::-1]
        loop_nodes_in_path = [
            n for n in path_nodes if n not in starts and n not in targets
        ]
        return True, len(loop_nodes_in_path), full_path

    @staticmethod
    def _has_stacking_path(
        starts: List[int],
        targets: List[int],
        valid_nodes: Set[int],
        stacking_index: Set[frozenset],
    ) -> bool:
        """Convenience method to check if a stacking path exists."""
        has_path, _, _ = Stacking._get_stacking_path_info(
            starts, targets, valid_nodes, stacking_index
        )
        return has_path

    @staticmethod
    def check_general_coaxiality(motif: dict, stacking_index: set) -> dict:
        """Determines if a bulge motif is BULGE-IN or BULGE-OUT based on stacking paths."""
        s = next(
            (s for s in motif.get("strands", []) if Core.loop_nt_count(s.get("structure", "")) > 0), None
        )
        if not s:
            return {
                "status": "ERROR",
                "metrics": {"size_2d": [0, 0], "size_3d": [0, 0], "stacking_path": []},
            }

        f, l = Utils.to_int(s.get("first")), Utils.to_int(s.get("last"))
        n_2d = Core.loop_nt_count(s.get("structure", "")) 
        valid_nodes = set(range(min(f, l), max(f, l) + 1))

        has_path, stacked_count, path = Stacking._get_stacking_path_info(
            [f], [l], valid_nodes, stacking_index
        )

        if has_path:
            status = "BULGE-OUT" if stacked_count == 0 else "BULGE-IN"
        else:
            status = "NO"

        return {
            "status": status,
            "metrics": {
                "size_2d": [n_2d, 0],
                "size_3d": [stacked_count, 0],
                "stacking_path": path,
            },
        }

    @staticmethod
    def check_internal_coaxiality(motif: dict, stacking_index: set) -> dict:
        """Determines if an internal loop has FULL, PARTIAL, or NO coaxial stacking."""
        strands = sorted(
            motif.get("strands", []), key=lambda x: Utils.to_int(x["first"])
        )
        if len(strands) < 2:
            return {
                "status": "ERROR",
                "metrics": {"size_2d": [0, 0], "size_3d": [0, 0], "stacking_path": []},
            }

        f1, l1 = (
            Utils.to_int(strands[0].get("first")),
            Utils.to_int(strands[0].get("last")),
        )
        f2, l2 = (
            Utils.to_int(strands[1].get("first")),
            Utils.to_int(strands[1].get("last")),
        )

        interface_up = [n for n in [f1, l2] if n is not None]
        interface_down = [n for n in [l1, f2] if n is not None]

        valid_nodes = set()
        s2d = []
        for s in strands:
            f, l = Utils.to_int(s.get("first")), Utils.to_int(s.get("last"))
            s2d.append(Core.loop_nt_count(s.get("structure", "")))
            if f is not None and l is not None:
                valid_nodes.update(range(min(f, l), max(f, l) + 1))
        valid_nodes.update(interface_up)
        valid_nodes.update(interface_down)

        success1, _, path1 = Stacking._get_stacking_path_info(
            interface_up, interface_down, valid_nodes, stacking_index
        )
        if not success1:
            return {
                "status": "NO",
                "metrics": {"size_2d": s2d, "size_3d": [0, 0], "stacking_path": []},
            }

        start_used, end_used = path1[0], path1[-1]
        other_start = (
            interface_up[1] if start_used == interface_up[0] else interface_up[0]
        )
        other_target = (
            interface_down[1] if end_used == interface_down[0] else interface_down[0]
        )

        remaining_nodes = valid_nodes - set(path1[1:-1])
        success2, _, path2 = Stacking._get_stacking_path_info(
            [other_start], [other_target], remaining_nodes, stacking_index
        )

        s3d = [
            len(path1) - 2 if len(path1) > 2 else 0,
            len(path2) - 2 if len(path2) > 2 else 0,
        ]

        if success2:
            return {
                "status": "FULL",
                "metrics": {
                    "size_2d": s2d,
                    "size_3d": s3d,
                    "stacking_path": [path1, path2],
                },
            }
        else:
            return {
                "status": "PARTIAL",
                "metrics": {
                    "size_2d": s2d,
                    "size_3d": [s3d[0], 0],
                    "stacking_path": [path1, []],
                },
            }

    @staticmethod
    def check_junction_stacking(motif: dict, stacking_index: set) -> dict:
        """Identifies pairs of stems that are coaxially stacked in a junction."""
        strands = motif.get("location", {}).get("strands", [])
        n = len(strands)
        if n < 3:
            return {"status": "ERROR", "coaxial_pairs": [], "stacking_paths": {}}
        valid_nodes = set()
        for s in strands:
            f, l = Utils.to_int(s.get("first")), Utils.to_int(s.get("last"))
            if f is not None and l is not None:
                valid_nodes.update(range(min(f, l), max(f, l) + 1))

        stem_ends = {}
        for i in range(1, n + 1):
            e1 = Utils.to_int(strands[i - 1].get("last"))
            e2 = Utils.to_int(strands[i % n].get("first"))
            stem_number = (i % n) + 1
            stem_ends[stem_number] = [x for x in (e1, e2) if x is not None]
            valid_nodes.update(stem_ends[stem_number])

        pairs = []
        stacking_paths = {}
        for i, j in itertools.combinations(range(1, n + 1), 2):
            has_path1, _, path1 = Stacking._get_stacking_path_info(
                stem_ends[i], stem_ends[j], valid_nodes, stacking_index
            )
            if not has_path1:
                continue

            found = [path1]
            remaining_nodes = valid_nodes - set(path1[1:-1])
            start_used, end_used = path1[0], path1[-1]
            other_starts = [x for x in stem_ends[i] if x != start_used]
            other_targets = [x for x in stem_ends[j] if x != end_used]
            if other_starts and other_targets:
                has_path2, _, path2 = Stacking._get_stacking_path_info(
                    other_starts, other_targets, remaining_nodes, stacking_index
                )
                if has_path2:
                    found.append(path2)

            pairs.append([f"stem_{i}", f"stem_{j}"])
            stacking_paths[f"stem_{i}_stem_{j}"] = found

        return {
            "status": "FULL" if pairs else "NO",
            "coaxial_pairs": pairs,
            "stacking_paths": stacking_paths,
        }


    @staticmethod
    def check_hairpin_stacking(motif: dict, stacking_index: set) -> dict:
        """Checks whether the hairpin loop has a genuine stacking path
        connecting its two flanking ends (f to l) through the loop's
        interior. FULL means such a path exists (the loop is stacked/bridged,
        whether directly or via intermediate - possibly non-canonical -
        nucleotides). NO means no such path exists (the stacking is broken)."""
        s = motif.get("strands", [{}])[0]
        f, l = Utils.to_int(s.get("first")), Utils.to_int(s.get("last"))
        n2d = Core.loop_nt_count(s.get("structure", ""))
        if f is None or l is None:
            return {
                "status": "ERROR",
                "metrics": {"size_2d": [0, 0], "size_3d": [0, 0], "stacking_path": []},
            }

        valid_nodes = set(range(min(f, l), max(f, l) + 1))
        has_path, stacked_count, path = Stacking._get_stacking_path_info(
            [f], [l], valid_nodes, stacking_index
        )

        return {
            "status": "FULL" if has_path else "NO",
            "metrics": {
                "size_2d": [n2d, 0],
                "size_3d": [stacked_count, 0],
                "stacking_path": path,
            },
        }


class GeometryCalculator:
    """Handles 3D geometric calculations using PyMOL and X3DNA."""

    _loaded_structure = None
    _axis_cache = {}
    _pymol_initialized = False

    @staticmethod
    def _init_pymol():
        """Initializes PyMOL in headless mode."""
        if not GeometryCalculator._pymol_initialized:
            pymol.finish_launching(["pymol", "-cq"])
            GeometryCalculator._pymol_initialized = True

    @staticmethod
    def load_structure(path: str):
        """Loads a CIF/PDB structure into PyMOL if not already loaded."""
        GeometryCalculator._init_pymol()
        if GeometryCalculator._loaded_structure == path:
            return
        cmd.reinitialize()
        cmd.load(path, "full_struct")
        GeometryCalculator._loaded_structure = path

    @staticmethod
    def extract_stem(pdb_path: str, selection_query: str, out_pdb_path: str):
        """Saves a specific structural selection to a temporary PDB file."""
        GeometryCalculator.load_structure(pdb_path)
        cmd.save(out_pdb_path, f"full_struct and ({selection_query})")

    @staticmethod
    def get_stem_selection(stem: Dict) -> Optional[str]:
        """Constructs a PyMOL selection string for an RNA stem."""
        if not stem or not stem.get("strand5p") or not stem.get("strand3p"):
            return None
        s5p_f, s5p_l = (
            stem["strand5p"].get("first", {}),
            stem["strand5p"].get("last", {}),
        )
        s3p_f, s3p_l = (
            stem["strand3p"].get("first", {}),
            stem["strand3p"].get("last", {}),
        )

        sel_5p = f"(chain {s5p_f.get('chain')} and resi {s5p_f.get('number')}-{s5p_l.get('number')})"
        sel_3p = f"(chain {s3p_f.get('chain')} and resi {s3p_f.get('number')}-{s3p_l.get('number')})"
        return f"{sel_5p} or {sel_3p}"

    @staticmethod
    def get_helical_axis(pdb_path: str) -> list:
        """Computes the helical axis vector using X3DNA find_pair and analyze."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_pdb = os.path.join(tmpdir, os.path.basename(pdb_path))
            with open(pdb_path, "r") as s, open(tmp_pdb, "w") as d:
                d.write(s.read())
            try:
                p1 = subprocess.Popen(
                    f"find_pair {tmp_pdb} stdout".split(),
                    stdout=subprocess.PIPE,
                    cwd=tmpdir,
                )
                p2 = subprocess.Popen(
                    "analyze".split(),
                    stdin=p1.stdout,
                    stdout=subprocess.PIPE,
                    cwd=tmpdir,
                )
                p1.stdout.close()
                p2.communicate()
            except Exception:
                return None
            out = tmp_pdb.replace(".pdb", ".out")
            if not os.path.exists(out):
                return None
            with open(out, "r") as f:
                for line in f:
                    if line.startswith("Helix:"):
                        p = line.split()
                        if len(p) >= 4:
                            return [float(p[1]), float(p[2]), float(p[3])]
        return None

    @staticmethod
    def get_stem_axis_cached(
        cif_path: str, stem_data: Dict, label: str = "tmp"
    ) -> Optional[List[float]]:
        """Retrieves helical axis with internal caching and orienting it 5' -> 3'."""
        sel = GeometryCalculator.get_stem_selection(stem_data)
        if not sel:
            return None
        cache_key = (cif_path, sel)
        if cache_key in GeometryCalculator._axis_cache:
            return GeometryCalculator._axis_cache[cache_key]

        tmp_pdb = f"{label}.pdb"
        try:
            GeometryCalculator.extract_stem(cif_path, sel, tmp_pdb)
            axis = GeometryCalculator.get_helical_axis(tmp_pdb)
            if os.path.exists(tmp_pdb):
                os.remove(tmp_pdb)

            if axis:
                # Orient axis 5' -> 3' based on C1' atoms
                s5p_f = stem_data["strand5p"].get("first", {})
                s5p_l = stem_data["strand5p"].get("last", {})
                sel_f = f"full_struct and chain {s5p_f.get('chain')} and resi {s5p_f.get('number')} and name C1'"
                sel_l = f"full_struct and chain {s5p_l.get('chain')} and resi {s5p_l.get('number')} and name C1'"

                GeometryCalculator.load_structure(cif_path)
                coords_f = cmd.get_coords(sel_f)
                coords_l = cmd.get_coords(sel_l)

                if coords_f is not None and len(coords_f) > 0 and coords_l is not None and len(coords_l) > 0:
                    v_ref = [
                        coords_l[0][0] - coords_f[0][0],
                        coords_l[0][1] - coords_f[0][1],
                        coords_l[0][2] - coords_f[0][2],
                    ]
                    dot_ref = sum(a * b for a, b in zip(axis, v_ref))
                    if dot_ref < 0:
                        axis = [-x for x in axis]

            GeometryCalculator._axis_cache[cache_key] = axis
            return axis
        except Exception:
            if os.path.exists(tmp_pdb):
                os.remove(tmp_pdb)
            return None

    @staticmethod
    def calculate_bend_angle(v1: list, v2: list) -> float:
        """Calculates the angle in degrees between two 3D vectors."""
        if not v1 or not v2:
            return None

        dot = max(-1.0, min(1.0, sum(a * b for a, b in zip(v1, v2))))
        return round(math.degrees(math.acos(dot)), 2)

    @staticmethod
    def validate_junction_compactness(loop_strands: List[Dict], mapping: Dict) -> bool:
        if not mapping:
            return True
        for s in loop_strands:
            f_ser, l_ser = Utils.to_int(s.get("first")), Utils.to_int(s.get("last"))
            if f_ser is None or l_ser is None:
                continue

            step = 1 if f_ser <= l_ser else -1
            expected_len = abs(l_ser - f_ser) + 1
            if expected_len <= 1:
                continue

            # Detect numbering direction (ascending or descending)
            try:
                n_start = Utils.to_int(
                    mapping.get(str(f_ser), {}).get("auth", {}).get("number")
                )
                n_end = Utils.to_int(
                    mapping.get(str(l_ser), {}).get("auth", {}).get("number")
                )
                auth_step = (
                    1
                    if (n_end is not None and n_start is not None and n_end >= n_start)
                    else -1
                )
            except Exception:
                auth_step = 1

            for i in range(expected_len - 1):
                curr = mapping.get(str(f_ser + i * step))
                nxt = mapping.get(str(f_ser + (i + 1) * step))
                if not curr or not nxt:
                    return False

                ac, an = curr["auth"], nxt["auth"]
                if ac["chain"] != an["chain"]:
                    return False

                n_curr, n_next = Utils.to_int(ac["number"]), Utils.to_int(an["number"])
                if n_curr is None or n_next is None:
                    return False

                # Strict DNA check: Reject if any nucleotide is DNA (DA, DC, DG, DT)
                if (ac.get("name") or "").upper() in ["DA", "DC", "DG", "DT"] or \
                   (an.get("name") or "").upper() in ["DA", "DC", "DG", "DT"]:
                    return False

                diff = n_next - n_curr
                if diff == 0:
                    ic, in_c = ac.get("icode") or "", an.get("icode") or ""
                    if ic == " ":
                        ic = ""
                    if in_c == " ":
                        in_c = ""
                    if auth_step == 1 and in_c <= ic:
                        return False
                    if auth_step == -1 and in_c >= ic:
                        return False
                elif diff == auth_step:
                    if auth_step == 1 and (an.get("icode") or "").strip() not in [
                        "",
                        ".",
                        "?",
                        " ",
                    ]:
                        return False
                    if auth_step == -1 and (ac.get("icode") or "").strip() not in [
                        "",
                        ".",
                        "?",
                        " ",
                    ]:
                        return False
                else:
                    return False
        return True

    @staticmethod
    def get_junction_bend_angles(
        cif_path: str, stems_data: Dict[str, Dict]
    ) -> Dict[str, float]:
        """Calculates all pair-wise bend angles between stems connected to a junction."""
        if not cif_path:
            return {}
        angles, axes = {}, []
        # Sort keys based on numeric part (stem1, stem2, ...)
        keys = sorted(
            stems_data.keys(), key=lambda x: int("".join(filter(str.isdigit, x)))
        )
        for k in keys:
            axes.append(
                GeometryCalculator.get_stem_axis_cached(
                    cif_path, stems_data[k], label=f"tmp_{k}"
                )
            )

        for i, j in itertools.combinations(range(len(keys)), 2):
            v1, v2 = axes[i], axes[j]
            # Use format stem_1_stem_2
            s1_id = "".join(filter(str.isdigit, keys[i]))
            s2_id = "".join(filter(str.isdigit, keys[j]))
            angles[f"stem_{s1_id}_stem_{s2_id}"] = (
                GeometryCalculator.calculate_bend_angle(v1, v2)
            )
        return angles


class Visualizer:
    """Handles 2D and 3D visualization exports (PyMOL scripts and VARNA JSON)."""

    # High-contrast palette for Helices
    HELIX_COLORS = [
        "#1890FF",
        "#A0D911",
        "#722ED1",
        "#13C2C2",
        "#EB2F96",
        "#237804",
        "#FAAD14",
        "#2F54EB",
    ]
    # Vibrant Warm Palette for Junction Core
    JUNCTION_COLORS = [
        "#FA8C16",
        "#FADB14",
        "#EB2F96",
        "#FA541C",
        "#A0D911",
        "#D4380D",
        "#C41D7F",
        "#722ED1",
    ]

    PYMOL_HELIX = [
        "blue",
        "lime",
        "purple",
        "cyan",
        "hotpink",
        "darkgreen",
        "gold",
        "royal",
    ]
    PYMOL_JUNCTION = [
        "orange",
        "yellow",
        "hotpink",
        "red",
        "lime",
        "darkorange",
        "magenta",
        "violet",
    ]

    @staticmethod
    def hex_to_rgb(hex_color: str) -> Dict[str, int]:
        """Converts hex color string to RGB dictionary."""
        hex_color = hex_color.lstrip("#")
        return {
            "r": int(hex_color[0:2], 16),
            "g": int(hex_color[2:4], 16),
            "b": int(hex_color[4:6], 16),
        }

    @staticmethod
    def prepare_3d_components(data_obj: Dict, is_junction: bool = False) -> List[Dict]:
        """Prepares a simplified 3D component list for frontend visualization."""
        comp3d = []

        def add_res(strand):
            if not strand or not strand.get("first"):
                return []
            f, l = strand["first"], strand["last"]
            return [
                {
                    "start": f["number"],
                    "end": l["number"],
                    "chain": f.get("chain", "A"),
                    "start_icode": f.get("icode"),
                    "end_icode": l.get("icode"),
                }
            ]

        def process_stems(target):
            if not is_junction:
                strands = target.get("strands", {})
                for s_t in ["upstream", "downstream"]:
                    s_d = strands.get(s_t, {})
                    res = add_res(s_d.get("strand5p")) + add_res(s_d.get("strand3p"))
                    if res:
                        comp3d.append(
                            {"id": s_t.upper(), "type": "STEM", "residues": res}
                        )
            else:
                for s_n, s_d in target.get("context", {}).items():
                    if "stem" in s_n:
                        res = add_res(s_d.get("strand5p")) + add_res(
                            s_d.get("strand3p")
                        )
                        if res:
                            comp3d.append(
                                {"id": s_n.upper(), "type": "STEM", "residues": res}
                            )

        def process_components(target, palette, prefix=""):
            for i, c in enumerate(target.get("components", [])):
                loc_res = []
                for loc in c.get("location", []):
                    loc_res.extend(add_res(loc))
                if loc_res:
                    c_type = c.get("type", "MOTIF").upper()
                    comp_id = f"{prefix}{c.get('type')}_{c.get('m_id', i + 1)}"
                    color = (
                        Visualizer.hex_to_rgb(palette[i % 8])
                        if "STEM" not in c_type
                        else None
                    )
                    comp3d.append(
                        {
                            "id": comp_id,
                            "type": c_type,
                            "residues": loc_res,
                            "color": color,
                        }
                    )
                i_stem = c.get("internal_stem", {})
                ires = add_res(i_stem.get("strand5p")) + add_res(i_stem.get("strand3p"))
                if ires:
                    comp3d.append(
                        {
                            "id": f"{prefix}INTERNAL_STEM_{i + 1}",
                            "type": "STEM",
                            "residues": ires,
                        }
                    )

        if is_junction:
            for i, loc in enumerate(data_obj.get("location", [])):
                res = add_res(loc)
                if res:
                    comp3d.append(
                        {
                            "id": f"LOOP_{i + 1}",
                            "type": "LOOP",
                            "residues": res,
                            "color": Visualizer.hex_to_rgb(
                                Visualizer.JUNCTION_COLORS[i % 8]
                            ),
                        }
                    )
            process_stems(data_obj)
            for ext in data_obj.get("extended_helices", []):
                ext_id = ext.get("h_id", "EXT")
                for s_t in ["upstream", "downstream"]:
                    s_d = ext.get("strands", {}).get(s_t, {})
                    res = add_res(s_d.get("strand5p")) + add_res(s_d.get("strand3p"))
                    if res:
                        comp3d.append(
                            {
                                "id": f"EXT_{ext_id}_{s_t.upper()}",
                                "type": "EXTENDED_STEM",
                                "residues": res,
                            }
                        )
                process_components(
                    ext, Visualizer.HELIX_COLORS, prefix=f"EXT_{ext_id}_"
                )
        else:
            process_stems(data_obj)
            process_components(data_obj, Visualizer.HELIX_COLORS)
        return comp3d

    @staticmethod
    def _get_single_sel(s: Dict) -> Optional[str]:
        """Creates a PyMOL selection for a single RNA strand."""
        if not s:
            return None
        f_obj, l_obj = s.get("first", {}), s.get("last", {})
        seq = s.get("sequence", "")
        c1 = f_obj.get("chain")
        c2 = l_obj.get("chain")
        f_ser, l_ser = f_obj.get("serial"), l_obj.get("serial")

        if c1 is None or f_ser is None or c2 is None or l_ser is None:
            return None

        # Check continuity using serial indices
        is_continuous = (c1 == c2) and (abs(f_ser - l_ser) + 1 == len(seq))

        def fmt_res(obj):
            num = obj.get("number")
            icode = obj.get("icode") or ""
            if icode in [None, ".", "?", " "]:
                icode = ""
            return f"{num}{icode}"

        if is_continuous:
            # resi 1406-1406C works in PyMOL if they are on the same chain
            return f"(chain {c1} and resi {fmt_res(f_obj)}-{fmt_res(l_obj)})"
        else:
            return f"(chain {c1} and resi {fmt_res(f_obj)}) or (chain {c2} and resi {fmt_res(l_obj)})"

    @staticmethod
    def get_continuous_selection(data: Dict, is_junction: bool = False) -> str:
        """Generates a complete PyMOL selection string for a complex motif or helix."""
        parts = []

        def add_s(s):
            sel = Visualizer._get_single_sel(s)
            if sel:
                parts.append(sel)

        if is_junction:
            for s in data.get("location", []):
                add_s(s)
            for stem in [v for k, v in data.get("context", {}).items() if "stem" in k]:
                add_s(stem.get("strand5p"))
                add_s(stem.get("strand3p"))
            for ext_h in data.get("extended_helices", []):
                for k in ["upstream", "downstream"]:
                    add_s(ext_h.get("strands", {}).get(k, {}).get("strand5p"))
                    add_s(ext_h.get("strands", {}).get(k, {}).get("strand3p"))
                for comp in ext_h.get("components", []):
                    for loc in comp.get("location", []):
                        add_s(loc)
                    add_s(comp.get("internal_stem", {}).get("strand5p"))
                    add_s(comp.get("internal_stem", {}).get("strand3p"))
        else:
            for k in ["upstream", "downstream"]:
                add_s(data.get("strands", {}).get(k, {}).get("strand5p"))
                add_s(data.get("strands", {}).get(k, {}).get("strand3p"))
            for comp in data.get("components", []):
                for loc in comp.get("location", []):
                    add_s(loc)
                add_s(comp.get("internal_stem", {}).get("strand5p"))
                add_s(comp.get("internal_stem", {}).get("strand3p"))
        return " or ".join(parts) if parts else "none"

    @staticmethod
    def generate_pml_script(
        pml_path: str,
        cif_path: str,
        data: Dict,
        full_sel: str,
        is_junction: bool = False,
    ):
        """Creates a .pml script for high-quality PyMOL rendering of the motif."""
        rel_cif = os.path.relpath(cif_path, os.path.dirname(pml_path)).replace(
            "\\", "/"
        )
        with open(pml_path, "w", encoding="UTF-8") as out:
            out.write(
                "reinitialize\nload "
                + rel_cif
                + ", rna\nbg_color white\nhide everything\n\n"
            )
            sname = "junc" if is_junction else "helix"
            out.write(
                f"select {sname}, {full_sel}\nshow sticks, {sname}\nshow cartoon, {sname}\nset cartoon_sampling, 20\ncolor gray80, {sname}\n\n"
            )
            if is_junction:
                # Core Loop
                for i, strand in enumerate(data.get("location", [])):
                    sel = Visualizer._get_single_sel(strand)
                    if sel:
                        out.write(
                            f"select s_{i}, {sel}\ncolor {Visualizer.PYMOL_JUNCTION[i % 8]}, s_{i}\n"
                        )
                # Extended Helices
                for h_idx, ext in enumerate(data.get("extended_helices", [])):
                    for c_idx, comp in enumerate(ext.get("components", [])):
                        m_sel = " or ".join(
                            [
                                Visualizer._get_single_sel(loc)
                                for loc in comp.get("location", [])
                                if Visualizer._get_single_sel(loc)
                            ]
                        )
                        if m_sel:
                            out.write(
                                f"select ext_{h_idx}_m_{c_idx}, {m_sel}\ncolor {Visualizer.PYMOL_HELIX[c_idx % 8]}, ext_{h_idx}_m_{c_idx}\n"
                            )
            else:
                for i, comp in enumerate(data.get("components", [])):
                    m_sel = " or ".join(
                        [
                            Visualizer._get_single_sel(loc)
                            for loc in comp.get("location", [])
                            if Visualizer._get_single_sel(loc)
                        ]
                    )
                    if m_sel:
                        out.write(
                            f"select m_{i}, {m_sel}\ncolor {Visualizer.PYMOL_HELIX[i % 8]}, m_{i}\n"
                        )
            out.write(f"zoom {sname}\norient {sname}\n")

    @staticmethod
    def export_varna_json(
        data: Dict, out_path: str, mapping: Dict = None, is_junction: bool = False
    ):
        """Exports the RNA secondary structure to a JSON format compatible with VARNA visualization."""
        nt_data = {}

        def add_s(s, color):
            if not s:
                return
            f, l, seq, struct = (
                Utils.to_int(s.get("first")),
                Utils.to_int(s.get("last")),
                s.get("sequence", ""),
                s.get("structure", ""),
            )

            cleaned_struct = "".join(
                c if c in "()" else "." for c in struct
            )

            if f is None or l is None:
                return
            step = 1 if f <= l else -1
            for i in range(len(seq)):
                ser = f + i * step
                label = str(ser)
                if mapping and str(ser) in mapping:
                    auth = mapping[str(ser)].get("auth", {})
                    num, icode = (
                        auth.get("number"),
                        auth.get("icode") or auth.get("insertion_code"),
                    )
                    if num is not None:
                        label = (
                            f"{num}{icode}"
                            if icode and str(icode).strip() not in [".", "?", ""]
                            else str(num)
                        )
                nt_data[ser] = {
                    "char": seq[i],
                    "struct": cleaned_struct[i] if i < len(cleaned_struct) else ".",
                    "orig_struct": struct[i] if i < len(struct) else ".",
                    "color": color,
                    "label": label,
                }

        if is_junction:
            stacked = {
                s
                for p in data.get("modules", {})
                .get("stacking", {})
                .get("coaxial_pairs", [])
                for s in p
            }
            for ext_h in data.get("extended_helices", []):
                for k in ["upstream", "downstream"]:
                    add_s(
                        ext_h.get("strands", {}).get(k, {}).get("strand5p"), "#D3D3D3"
                    )
                    add_s(
                        ext_h.get("strands", {}).get(k, {}).get("strand3p"), "#D3D3D3"
                    )
                for comp in ext_h.get("components", []):
                    add_s(comp.get("internal_stem", {}).get("strand5p"), "#D3D3D3")
                    add_s(comp.get("internal_stem", {}).get("strand3p"), "#D3D3D3")
            for k, stem in data.get("context", {}).items():
                if "stem" in k:
                    c = "#555555" if k in stacked else "#D3D3D3"
                    add_s(stem.get("strand5p"), c)
                    add_s(stem.get("strand3p"), c)
            for i, s in enumerate(data.get("location", [])):
                add_s(s, Visualizer.JUNCTION_COLORS[i % 8])
            for ext_h in data.get("extended_helices", []):
                for i, comp in enumerate(ext_h.get("components", [])):
                    for loc in comp.get("location", []):
                        add_s(loc, Visualizer.HELIX_COLORS[i % 8])
        else:
            for k in ["upstream", "downstream"]:
                add_s(data.get("strands", {}).get(k, {}).get("strand5p"), "#D3D3D3")
                add_s(data.get("strands", {}).get(k, {}).get("strand3p"), "#D3D3D3")
            for comp in data.get("components", []):
                add_s(comp.get("internal_stem", {}).get("strand5p"), "#D3D3D3")
                add_s(comp.get("internal_stem", {}).get("strand3p"), "#D3D3D3")
            for i, comp in enumerate(data.get("components", [])):
                for loc in comp.get("location", []):
                    add_s(loc, Visualizer.HELIX_COLORS[i % 8])

        sorted_serials = sorted(nt_data.keys())
        counts = Counter(nt_data[s]["label"] for s in nt_data)
        nucleotides, ser_to_v_id = [], {}
        for i, ser in enumerate(sorted_serials):
            v_id = i + 1
            ser_to_v_id[ser] = v_id
            lbl = nt_data[ser]["label"]
            nucleotides.append(
                {
                    "id": v_id,
                    "number": lbl if counts[lbl] == 1 else str(ser),
                    "char": nt_data[ser]["char"],
                    "innerColor": nt_data[ser]["color"],
                    "outlineColor": "black",
                }
            )

        base_pairs = []
        stacks = {}
        close_to_open = {')': '(', ']': '[', '}': '{', '>': '<'}
        for i in range(26):
            close_to_open[chr(ord('a') + i)] = chr(ord('A') + i)

        for ser in sorted_serials:
            c = nt_data[ser].get("orig_struct", ".")
            
            if c in "([{<" or (c.isalpha() and c.isupper()):
                stacks.setdefault(c, []).append(ser)
            elif c in ")]}>" or (c.isalpha() and c.islower()):
                open_c = close_to_open.get(c)
                if open_c and stacks.get(open_c):
                    p_ser = stacks[open_c].pop()
                    base_pairs.append({
                        "id1": ser_to_v_id[p_ser],
                        "id2": ser_to_v_id[ser],
                        "edge5": "WC",
                        "edge3": "WC",
                        "stericity": "CIS",
                        "canonical": True,
                    })


        lw_s, lw_e = (
            {"c": "CIS", "t": "TRANS"},
            {"w": "WC", "h": "HOOGSTEEN", "s": "SUGAR"},
        )
        interactions = []
        if is_junction:
            interactions.extend(
                data.get("modules", {}).get("interactions", {}).get("details", [])
            )
            for ext_h in data.get("extended_helices", []):
                for comp in ext_h.get("components", []):
                    interactions.extend(comp.get("interactions", []))
        else:
            for comp in data.get("components", []):
                interactions.extend(comp.get("interactions", []))

        for inter in interactions:
            s1, s2, lw = (
                Utils.to_int(inter.get("nt1")),
                Utils.to_int(inter.get("nt2")),
                inter.get("lw", "").strip(),
            )
            id1, id2 = ser_to_v_id.get(s1), ser_to_v_id.get(s2)
            if id1 and id2 and len(lw) >= 3:
                base_pairs.append(
                    {
                        "id1": id1,
                        "id2": id2,
                        "stericity": lw_s.get(lw[0].lower(), "CIS"),
                        "edge5": lw_e.get(lw[1].lower(), "WC"),
                        "edge3": lw_e.get(lw[2].lower(), "WC"),
                        "canonical": False,
                    }
                )

        # --- DRAW STACKING PATHS (JUMPS) ---
        v_stackings = []
        core_stack_nodes = set()

        def add_jump_stacks(p_list, track=True):
            if not p_list or not isinstance(p_list, list):
                return
            if p_list and isinstance(p_list[0], list):
                for sub_path in p_list:
                    add_jump_stacks(sub_path, track)
                return
            for i in range(len(p_list) - 1):
                n1, n2 = Utils.to_int(p_list[i]), Utils.to_int(p_list[i + 1])
                if n1 is None or n2 is None:
                    continue
                if track:
                    core_stack_nodes.add(n1)
                    core_stack_nodes.add(n2)
                if abs(n1 - n2) > 1:
                    vid1, vid2 = ser_to_v_id.get(n1), ser_to_v_id.get(n2)
                    if vid1 and vid2:
                        v_stackings.append(
                            {"id1": vid1, "id2": vid2, "color": "255,0,0"}
                        )

        if is_junction:
            for path in (
                data.get("modules", {})
                .get("stacking", {})
                .get("stacking_paths", {})
                .values()
            ):
                add_jump_stacks(path, track=True)
            for ext_h in data.get("extended_helices", []):
                for comp in ext_h.get("components", []):
                    if comp.get("type") in ["BULGE", "INTERNAL_LOOP"]:
                        add_jump_stacks(
                            comp.get("metrics", {}).get("stacking_path"), track=False
                        )
            for ser in core_stack_nodes:
                v_id = ser_to_v_id.get(ser)
                if v_id:
                    nucleotides[v_id - 1]["innerColor"] = "255,0,0"
        else:
            for comp in data.get("components", []):
                if comp.get("type") in ["BULGE", "INTERNAL_LOOP"]:
                    add_jump_stacks(comp.get("metrics", {}).get("stacking_path"), track=False)

        with open(out_path, "w", encoding="UTF-8") as f:
            json.dump(
                {
                    "drawingAlgorithm": "NAVIEW",
                    "stackingArrowPlacement": "opposing",
                    "stackingArrowGap": 8.0,
                    "nucleotides": nucleotides,
                    "basePairs": base_pairs,
                    "stackings": v_stackings,
                },
                f,
                indent=4,
            )


class HelicesBuilder:
    """Builds super-helices from motifs and validates global geometry."""

    @staticmethod
    def get_global_angle(
        cif_path: str, motif_start: Dict, motif_end: Dict
    ) -> Optional[float]:
        """Calculates the bend angle between the first and last stem of a super-helix."""
        ctx_start = (
            motif_start.get("location", {}).get("context", {}).get("upstream", {})
        )
        if motif_end.get("meta", {}).get("type") == "HAIRPIN":
            ctx_end = (
                motif_end.get("location", {}).get("context", {}).get("upstream", {})
            )
        else:
            ctx_end = (
                motif_end.get("location", {}).get("context", {}).get("downstream", {})
            )

        def stem_key(stem):
            f = stem.get("strand5p", {}).get("first", {})
            return f.get("serial") if isinstance(f, dict) else f

        k1, k2 = stem_key(ctx_start), stem_key(ctx_end)
        if k1 is not None and k1 == k2:
            return 0.0

        v1 = GeometryCalculator.get_stem_axis_cached(
            cif_path, ctx_start, label="tmp_g1"
        )
        v2 = GeometryCalculator.get_stem_axis_cached(cif_path, ctx_end, label="tmp_g2")
        return GeometryCalculator.calculate_bend_angle(v1, v2)

    @staticmethod
    def stems_match(d: Dict, u: Dict) -> bool:
        """Checks if two stems are consecutive based on their serial indices."""
        try:
            return Utils.to_int(d.get("strand5p", {}).get("first")) == Utils.to_int(
                u.get("strand5p", {}).get("first")
            ) and Utils.to_int(d.get("strand5p", {}).get("last")) == Utils.to_int(
                u.get("strand5p", {}).get("last")
            )
        except Exception:
            return False

    @staticmethod
    def extract_helices(motifs: List[Dict]) -> List[Dict]:
        """Groups motifs into raw helix sequences based on stem connectivity."""
        v = [
            m
            for m in motifs
            if "WAY_JUNCTION" not in (m.get("meta", {}).get("type") or "")
            and (
                m.get("modules", {}).get("stacking", {}).get("status")
                in ["FULL", "BULGE-IN", "BULGE-OUT"]
                or (
                    m.get("meta", {}).get("type") == "HAIRPIN"
                    and m.get("modules", {}).get("stacking", {}).get("status") == "FULL"
                )
            )
        ]
        if not v:
            return []
        v.sort(
            key=lambda x: (
                Utils.to_int(x.get("location", {}).get("strands", [{}])[0].get("first"))
                or 0
            )
        )
        h, cur = [], [v[0]]
        for i in range(1, len(v)):
            prev_down = (
                cur[-1].get("location", {}).get("context", {}).get("downstream", {})
            )
            next_up = v[i].get("location", {}).get("context", {}).get("upstream", {})
            if cur[-1].get("meta", {}).get(
                "type"
            ) == "HAIRPIN" or not HelicesBuilder.stems_match(prev_down, next_up):
                h.append(cur)
                cur = [v[i]]
            else:
                cur.append(v[i])
        h.append(cur)
        return HelicesBuilder._format_helices(h)

    @staticmethod
    def _format_helices(groups: List[List[Dict]]) -> List[Dict]:
        """Formats grouped motifs into the final helix JSON structure."""
        res = []
        for idx, g in enumerate(groups, 1):
            if len(g) < 1:
                continue
            f, l = g[0], g[-1]
            s1f = (
                f.get("location", {})
                .get("context", {})
                .get("upstream", {})
                .get("strand5p", {})
                .get("first")
            )
            s1l = (
                l.get("location", {})
                .get("context", {})
                .get("downstream", {})
                .get("strand5p", {})
                .get("last")
            )
            s2f = (
                l.get("location", {})
                .get("context", {})
                .get("downstream", {})
                .get("strand3p", {})
                .get("first")
            )
            s2l = (
                f.get("location", {})
                .get("context", {})
                .get("upstream", {})
                .get("strand3p", {})
                .get("last")
            )
            total = (
                (abs(Utils.to_int(s1l) - Utils.to_int(s1f)) + 1)
                + (abs(Utils.to_int(s2l) - Utils.to_int(s2f)) + 1)
                if all(
                    x is not None
                    for x in [
                        Utils.to_int(s1f),
                        Utils.to_int(s1l),
                        Utils.to_int(s2f),
                        Utils.to_int(s2l),
                    ]
                )
                else 0
            )
            comps = []
            for m in g:
                m_metrics = m.get("modules", {}).get("stacking", {}).get("metrics", {})
                comps.append(
                    {
                        "m_id": m.get("meta", {}).get("id"),
                        "type": m.get("meta", {}).get("type"),
                        "location": m.get("location", {}).get("strands", []),
                        "metrics": {
                            "size_2d": m_metrics.get("size_2d"),
                            "size_3d": m_metrics.get("size_3d"),
                            "stacking_path": m_metrics.get("stacking_path", []),
                            "stacking": m.get("modules", {})
                            .get("stacking", {})
                            .get("status"),
                            "bend_angle": m.get("modules", {})
                            .get("geometry", {})
                            .get("bend_angle"),
                        },
                        "interactions": m.get("modules", {})
                        .get("interactions", {})
                        .get("details", []),
                        "internal_stem": m.get("location", {})
                        .get("context", {})
                        .get("downstream", {}),
                    }
                )
            res.append(
                {
                    "h_id": idx,
                    "total_nt": total,
                    "strands": {
                        "upstream": {"first": s1f, "last": s1l},
                        "downstream": {"first": s2f, "last": s2l},
                    },
                    "components": comps,
                }
            )
        return res
