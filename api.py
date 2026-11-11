import time
import urllib.request
import os
import json
import csv
import io
import zipfile
from threading import Lock
from fastapi import FastAPI, Query, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, PlainTextResponse, StreamingResponse
from sqlalchemy import create_engine, or_, and_, func, String, Integer
from sqlalchemy.orm import sessionmaker
from sqlalchemy import tuple_
from database import Helix, Segment, Junction
from pathlib import Path
from config import DATABASE_URL, OUTPUT_DIR, BASE_DIR
from rnabridge import Visualizer
from typing import List, Optional, Dict, Any, Set

app = FastAPI()

# Enable CORS for frontend interaction
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Database initialization
engine = (
    create_engine(DATABASE_URL, pool_pre_ping=True)
    if DATABASE_URL.startswith("postgresql")
    else create_engine(DATABASE_URL)
)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

OUTPUT_ROOTS = {
    Path(OUTPUT_DIR).resolve(),
    (BASE_DIR / "result").resolve(),
    (BASE_DIR / "rnasolo" / "result").resolve(),
}


def resolve_output_file(file_path: str) -> Path:
    """Resolve a path only if it stays within an allowed results directory."""
    candidate = Path(file_path)
    if not candidate.is_absolute():
        candidate = Path.cwd() / candidate

    resolved = candidate.resolve()

    for root in OUTPUT_ROOTS:
        try:
            resolved.relative_to(root)
            return resolved
        except ValueError:
            continue

    raise HTTPException(status_code=404)

NR_CACHE_FILE = "non_redundant_cache.json"
_nr_memory_cache = {"data": {}, "timestamp": 0}

_NR_ALLOWED_PAIRS_CACHE = {"map_id": None, "pairs": frozenset()}
_NR_ALLOWED_PAIRS_LOCK = Lock()


def get_nr_allowed_pairs(session, nr_map: Dict[str, List[Set[str]]]) -> Set[tuple]:
    map_id = id(nr_map)

    with _NR_ALLOWED_PAIRS_LOCK:
        if _NR_ALLOWED_PAIRS_CACHE["map_id"] == map_id:
            return _NR_ALLOWED_PAIRS_CACHE["pairs"]

        allowed_pairs = set()

        for model in (Helix, Junction):
            rows = (
                session.query(model.pdb_id, model.chains)
                .filter(model.pdb_id.isnot(None), model.chains.isnot(None))
                .distinct()
                .all()
            )

            for pdb_id, chains_str in rows:
                motif_chains = {
                    chain.strip()
                    for chain in str(chains_str).split("-")
                    if chain.strip()
                }
                if not motif_chains:
                    continue

                representative_groups = nr_map.get(str(pdb_id).upper(), [])
                if any(
                    motif_chains.issubset(rep_group)
                    for rep_group in representative_groups
                ):
                    allowed_pairs.add((str(pdb_id).upper(), str(chains_str)))

        result = frozenset(allowed_pairs)
        _NR_ALLOWED_PAIRS_CACHE["map_id"] = map_id
        _NR_ALLOWED_PAIRS_CACHE["pairs"] = result
        return result

def get_non_redundant_map() -> Dict[str, List[Set[str]]]:
    """
    Fetches the BGSU representative set from RNAsolo and returns a dict:
    { 'PDB_ID': [set_of_chains_1, set_of_chains_2, ...] }
    Preserves exact chain case (e.g., 'c' vs 'C') and multi-chain complexes (e.g., 'A-B-C-D-E-F-G-H').
    """
    now = time.time()
    if _nr_memory_cache["data"] and (now - _nr_memory_cache["timestamp"] < 86400):
        return _nr_memory_cache["data"]

    try:
        ver_url = "https://rnasolo.cs.put.poznan.pl/api/statistics/version/"
        with urllib.request.urlopen(ver_url, timeout=10) as resp:
            ver_data = json.loads(resp.read().decode("utf-8"))
        bgsu_ver = ver_data.get("bgsu_ver", "").replace(".", "_")
        
        if bgsu_ver and os.path.exists(NR_CACHE_FILE):
            try:
                with open(NR_CACHE_FILE, "r") as f:
                    cached = json.load(f)

                if cached.get("bgsu_ver") == bgsu_ver:
                    raw_map = cached.get("map", {})
                    nr_map = {
                        k: [set(chains) for chains in groups]
                        for k, groups in raw_map.items()
                    }
                    if nr_map:
                        _nr_memory_cache["data"] = nr_map
                        _nr_memory_cache["timestamp"] = now
                        return nr_map
            except Exception:
                pass

        if bgsu_ver:
            txt_url = f"https://rnasolo.cs.put.poznan.pl/media/files/bunches/BGSU__R__All__All__1000_0__txt_{bgsu_ver}.txt"
            with urllib.request.urlopen(txt_url, timeout=15) as resp:
                lines = resp.read().decode("utf-8").splitlines()

            nr_map: Dict[str, List[Set[str]]] = {}
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("_")
                if parts[0] == "PDB" and len(parts) >= 4:
                    pdb_id = parts[1][-4:].upper()
                    chain_part = "_".join(parts[3:])
                elif len(parts) >= 3:
                    pdb_id = parts[0][-4:].upper()
                    chain_part = "_".join(parts[2:])
                else:
                    continue

                chain_set = {c for c in chain_part.split("-") if c}
                if chain_set:
                    nr_map.setdefault(pdb_id, []).append(chain_set)

            if nr_map:
                _nr_memory_cache["data"] = nr_map
                _nr_memory_cache["timestamp"] = now
                try:
                    serializable = {k: [list(s) for s in v] for k, v in nr_map.items()}
                    with open(NR_CACHE_FILE, "w") as f:
                        json.dump({"bgsu_ver": bgsu_ver, "map": serializable}, f)
                except Exception:
                    pass
                return nr_map
    except Exception as e:
        print(f"[WARNING] Could not fetch non-redundant list from RNAsolo: {e}")

    if os.path.exists(NR_CACHE_FILE):
        try:
            with open(NR_CACHE_FILE, "r") as f:
                raw_map = json.load(f).get("map", {})
                nr_map = {k: [set(s) for s in v] for k, v in raw_map.items()}
                _nr_memory_cache["data"] = nr_map
                _nr_memory_cache["timestamp"] = now
                return nr_map
        except Exception:
            pass

    return _nr_memory_cache["data"]

def format_helix_response(h: Helix, session) -> Dict[str, Any]:
    """Formats a Helix database object into a detailed API response with cleaned motif sequences."""
    db_segments = (
        session.query(Segment).filter(Segment.helix_id == h.id).order_by(Segment.s_id).all()
    )

    clean_sequences = {}
    comp3d = []

    if h.path_json and os.path.exists(h.path_json):
        try:
            with open(h.path_json, "r") as f:
                data = json.load(f).get("helices", [])
                if data:
                    h_data = data[0]
                    comp3d = Visualizer.prepare_3d_components(h_data, is_junction=False)

                    curr_gid = 1
                    if h_data.get("strands", {}).get("upstream", {}):
                        curr_gid += 1

                    for c in h_data.get("components", []):
                        c_seq_parts = []
                        for loc in c.get("location", []):
                            seq, struct = loc.get("sequence", ""), loc.get("structure", "")
                            if seq and struct and len(seq) == len(struct):
                                clean_seq = "".join(
                                    [seq[idx] for idx in range(len(seq)) if struct[idx] not in "()"]
                                )
                                if clean_seq:
                                    c_seq_parts.append(clean_seq)
                            else:
                                if seq:
                                    c_seq_parts.append(seq)

                        clean_sequences[curr_gid] = "-".join(c_seq_parts)
                        curr_gid += 1
                        # Skip internal stem of motif if exists
                        i_stem = c.get("internal_stem", {})
                        if i_stem and (
                            i_stem.get("strand5p", {}).get("sequence")
                            or i_stem.get("strand3p", {}).get("sequence")
                        ):
                            curr_gid += 1
        except Exception:
            pass

    segments_list = []
    for s in db_segments:
        try:
            path = json.loads(s.stacking_path) if s.stacking_path else []
        except Exception:
            path = []

        display_seq = s.sequence_full
        if s.type.upper() != "STEM" and s.s_id in clean_sequences:
            display_seq = clean_sequences[s.s_id]

        segments_list.append(
            {
                "s_id": s.s_id,
                "type": s.type,
                "size_2d": f"{s.size_2d_left}x{s.size_2d_right}",
                "size_3d": f"{s.size_3d_left}x{s.size_3d_right}",
                "bend_angle": s.bend_angle,
                "stacking": s.stacking,
                "stacking_path": path,
                "sequence": display_seq,
            }
        )

    return {
        "id": h.id,
        "pdb_id": h.pdb_id,
        "molecule": h.molecule,
        "organism": h.organism,
        "resolution": h.resolution,
        "method": h.method,
        "total_nt": h.total_nt,
        "global_bend_angle": h.global_bend_angle,
        "segment_count_folder": h.segment_count_folder,
        "path_json": h.path_json,
        "path_svg": h.path_svg,
        "path_cif": h.path_cif,
        "path_pml": h.path_pml,
        "type": "helix",
        "details": {"db_segments": segments_list, "components": comp3d},
    }


def format_junction_response(j: Junction) -> Dict[str, Any]:
    """Formats a Junction database object into a detailed API response."""
    try:
        raw_angles = json.loads(j.bend_angles) if j.bend_angles else {}
        angles = {k.replace("_stem_", "_"): v for k, v in raw_angles.items()}
    except Exception:
        angles = {}
    try:
        pairs = json.loads(j.coaxial_pairs) if j.coaxial_pairs else []
    except Exception:
        pairs = []
    try:
        raw_paths = json.loads(j.stacking_paths) if j.stacking_paths else {}
        paths = {k.replace("_stem_", "_"): v for k, v in raw_paths.items()}
    except Exception:
        paths = {}

    stems_info, comp3d = {}, []
    if j.path_json and os.path.exists(j.path_json):
        try:
            with open(j.path_json, "r") as f:
                data = json.load(f).get("junctions", [])
                if data:
                    target = next(
                        (jd for jd in data if int(jd.get("j_id", 0)) == (j.j_id or 1)),
                        data[0],
                    )
                    for s_n, s_d in target.get("context", {}).items():
                        if "stem" in s_n:
                            stems_info[s_n.lower()] = {
                                "strand5p": s_d.get("strand5p"),
                                "strand3p": s_d.get("strand3p"),
                            }
                    comp3d = Visualizer.prepare_3d_components(target, is_junction=True)
        except Exception:
            pass

    return {
        "id": j.id,
        "pdb_id": j.pdb_id,
        "molecule": j.molecule,
        "organism": j.organism,
        "resolution": j.resolution,
        "method": j.method,
        "total_nt": j.total_nt,
        "global_bend_angle": j.global_bend_angle,
        "segment_count_folder": j.segment_count_folder,
        "path_json": j.path_json,
        "path_svg": j.path_svg,
        "path_cif": j.path_cif,
        "path_pml": j.path_pml,
        "type": "junction",
        "details": {
            "stacking_status": j.stacking_status,
            "coaxial_pairs": pairs,
            "all_angles": angles,
            "stacking_paths": paths,
            "stems": stems_info,
            "components": comp3d,
        },
    }


def apply_filters(h_q, j_q, filters: Dict[str, Any]):
    """Applies a consistent set of filters to both Helix and Junction queries."""
    
    methods = filters.get("method")
    if methods:
        h_q = h_q.filter(Helix.method.in_(methods))
        j_q = j_q.filter(Junction.method.in_(methods))

    actual = filters.get("segment_type") or filters.get("motif_type") or []
    if actual:
        h_q = h_q.filter(or_(*[Helix.segment_count_folder.contains(t) for t in actual]))
        j_f = []
        for t in actual:
            if t == "8plus-junctions":
                j_f.append(
                    and_(
                        ~Junction.segment_count_folder.contains("3-way"),
                        ~Junction.segment_count_folder.contains("4-way"),
                        ~Junction.segment_count_folder.contains("5-way"),
                        ~Junction.segment_count_folder.contains("6-way"),
                        ~Junction.segment_count_folder.contains("7-way"),
                        Junction.segment_count_folder.contains("way-junction"),
                    )
                )
            else:
                j_f.append(Junction.segment_count_folder.contains(t))
        j_q = j_q.filter(or_(*j_f))

    min_angle = filters.get("min_angle")
    if min_angle is not None:
        h_q, j_q = (
            h_q.filter(Helix.global_bend_angle >= min_angle),
            j_q.filter(Junction.global_bend_angle >= min_angle),
        )

    max_angle = filters.get("max_angle")
    if max_angle is not None:
        h_q, j_q = (
            h_q.filter(Helix.global_bend_angle <= max_angle),
            j_q.filter(Junction.global_bend_angle <= max_angle),
        )

    min_nt = filters.get("min_nt")
    if min_nt is not None:
        h_q, j_q = (
            h_q.filter(Helix.total_nt >= min_nt),
            j_q.filter(Junction.total_nt >= min_nt),
        )

    max_nt = filters.get("max_nt")
    if max_nt is not None:
        h_q, j_q = (
            h_q.filter(Helix.total_nt <= max_nt),
            j_q.filter(Junction.total_nt <= max_nt),
        )

    search_pdb = filters.get("search_pdb")
    if search_pdb:
        pdb_id = search_pdb.strip().upper()
        h_q = h_q.filter(func.upper(Helix.pdb_id) == pdb_id)
        j_q = j_q.filter(func.upper(Junction.pdb_id) == pdb_id)

    sequence = filters.get("sequence")
    if sequence:
        s_up, s_rev = sequence.upper(), sequence.upper()[::-1]
        h_q, j_q = (
            h_q.filter(
                or_(
                    Helix.sequence_full.contains(s_up),
                    Helix.sequence_full.contains(s_rev),
                )
            ),
            j_q.filter(
                or_(
                    Junction.sequence.contains(s_up),
                    Junction.sequence.contains(s_rev),
                )
            ),
        )

    s1, s2 = filters.get("stacking_stem1"), filters.get("stacking_stem2")
    if s1 or s2:
        h_q = h_q.filter(Helix.id < 0)

        stem1 = s1.strip().lower() if s1 else None
        stem2 = s2.strip().lower() if s2 else None

        matching_ids = []
        rows = j_q.with_entities(Junction.id, Junction.coaxial_pairs).all()

        for junction_id, raw_pairs in rows:
            try:
                pairs = (
                    json.loads(raw_pairs)
                    if isinstance(raw_pairs, str)
                    else raw_pairs
                )
            except (TypeError, json.JSONDecodeError):
                continue

            if not isinstance(pairs, list):
                continue

            normalized_pairs = [
                [str(stem).strip().lower() for stem in pair]
                for pair in pairs
                if isinstance(pair, list)
            ]

            if stem1 and stem2:
                matches = any(
                    len(pair) == 2 and set(pair) == {stem1, stem2}
                    for pair in normalized_pairs
                )
            else:
                selected_stem = stem1 or stem2
                matches = any(
                    selected_stem in pair for pair in normalized_pairs
                )

            if matches:
                matching_ids.append(junction_id)

        j_q = j_q.filter(Junction.id.in_(matching_ids))

    dataset = "rnasolo" if filters.get("non_redundant") else "pdb"
    h_q = h_q.filter(Helix.dataset == dataset)
    j_q = j_q.filter(Junction.dataset == dataset)

    return h_q, j_q


@app.get("/api/search")
async def search(
    segment_type: Optional[List[str]] = Query(None),
    motif_type: Optional[List[str]] = Query(None),
    method: Optional[List[str]] = Query(None),
    min_angle: Optional[float] = Query(None),
    max_angle: Optional[float] = Query(None),
    min_nt: Optional[int] = Query(None),
    max_nt: Optional[int] = Query(None),
    search_pdb: Optional[str] = Query(None),
    sequence: Optional[str] = Query(None),
    stacking_stem1: Optional[str] = Query(None),
    stacking_stem2: Optional[str] = Query(None),
    sort_by: Optional[str] = Query(None),
    sort_order: str = Query("asc"),
    page: int = Query(1, ge=1),
    limit: int = Query(10, ge=1, le=100),
    non_redundant: Optional[bool] = Query(False),
):
    """General search endpoint for helices and junctions."""
    with SessionLocal() as session:
        offset = (page - 1) * limit
        h_q, j_q = apply_filters(
            session.query(Helix),
            session.query(Junction),
            {
                "segment_type": segment_type,
                "motif_type": motif_type,
                "method": method,
                "min_angle": min_angle,
                "max_angle": max_angle,
                "min_nt": min_nt,
                "max_nt": max_nt,
                "search_pdb": search_pdb,
                "sequence": sequence,
                "stacking_stem1": stacking_stem1,
                "stacking_stem2": stacking_stem2,
                "non_redundant": non_redundant,
            },
        )

        # Combined sorting and pagination
        s_col = sort_by or "pdb_id"
        a_h, a_j = (
            getattr(Helix, s_col, Helix.pdb_id),
            getattr(Junction, s_col, Junction.pdb_id),
        )
        h_ids = h_q.with_entities(Helix.id, a_h).all()
        j_ids = j_q.with_entities(Junction.id, a_j).all()
        combined = [(i[0], i[1], "h") for i in h_ids] + [(i[0], i[1], "j") for i in j_ids]
        combined.sort(
            key=lambda x: (
                x[1]
                if x[1] is not None
                else ("" if isinstance(a_h.type, String) else -1.0)
            ),
            reverse=(sort_order == "desc"),
        )

        total_count = len(combined)
        page_items = combined[offset : offset + limit]
        final_results = []
        for item in page_items:
            if item[2] == "h":
                h_obj = session.query(Helix).filter(Helix.id == item[0]).first()
                if h_obj:
                    final_results.append(format_helix_response(h_obj, session))
            else:
                j_obj = session.query(Junction).filter(Junction.id == item[0]).first()
                if j_obj:
                    final_results.append(format_junction_response(j_obj))

        # Statistics for the dashboard
        pie_stats = {}
        for f, c in (
            h_q.with_entities(Helix.segment_count_folder, func.count(Helix.id))
            .group_by(Helix.segment_count_folder)
            .all()
        ):
            pie_stats[f] = pie_stats.get(f, 0) + c
        for f, c in (
            j_q.with_entities(Junction.segment_count_folder, func.count(Junction.id))
            .group_by(Junction.segment_count_folder)
            .all()
        ):
            pie_stats[f] = pie_stats.get(f, 0) + c

        angle_stats = []
        h_ang = (
            session.query(
                Helix.segment_count_folder,
                func.cast(func.floor(Helix.global_bend_angle / 5) * 5, Integer).label("bin"),
                func.count(Helix.id),
            )
            .filter(Helix.id.in_(h_q.with_entities(Helix.id)))
            .group_by(Helix.segment_count_folder, "bin")
            .all()
        )
        for f, b, c in h_ang:
            angle_stats.append({"folder": f, "bin": b, "count": c})
        j_ang = (
            session.query(
                Junction.segment_count_folder,
                func.cast(func.floor(Junction.global_bend_angle / 5) * 5, Integer).label("bin"),
                func.count(Junction.id),
            )
            .filter(Junction.id.in_(j_q.with_entities(Junction.id)))
            .group_by(Junction.segment_count_folder, "bin")
            .all()
        )
        for f, b, c in j_ang:
            angle_stats.append({"folder": f, "bin": b, "count": c})

        return {
            "results": final_results,
            "total": total_count,
            "page": page,
            "limit": limit,
            "stats": {"pie": pie_stats, "angles": angle_stats},
        }


@app.get("/api/export-csv")
async def export_csv(
    segment_type: Optional[List[str]] = Query(None),
    motif_type: Optional[List[str]] = Query(None),
    method: Optional[List[str]] = Query(None),
    min_angle: Optional[float] = Query(None),
    max_angle: Optional[float] = Query(None),
    min_nt: Optional[int] = Query(None),
    max_nt: Optional[int] = Query(None),
    search_pdb: Optional[str] = Query(None),
    sequence: Optional[str] = Query(None),
    stacking_stem1: Optional[str] = Query(None),
    stacking_stem2: Optional[str] = Query(None),
    sort_by: Optional[str] = Query(None),
    sort_order: str = Query("asc"),
    non_redundant: Optional[bool] = Query(False),
):
    """Exports all filtered results as a CSV file."""
    with SessionLocal() as session:
        h_q, j_q = apply_filters(
            session.query(Helix),
            session.query(Junction),
            {
                "segment_type": segment_type,
                "motif_type": motif_type,
                "method": method,
                "min_angle": min_angle,
                "max_angle": max_angle,
                "min_nt": min_nt,
                "max_nt": max_nt,
                "search_pdb": search_pdb,
                "sequence": sequence,
                "stacking_stem1": stacking_stem1,
                "stacking_stem2": stacking_stem2,
                "non_redundant": non_redundant,
            },
        )

        # Fetch results
        s_col = sort_by or "pdb_id"
        a_h, a_j = (
            getattr(Helix, s_col, Helix.pdb_id),
            getattr(Junction, s_col, Junction.pdb_id),
        )
        h_res = h_q.with_entities(
            Helix.pdb_id,
            Helix.organism,
            Helix.method,
            Helix.resolution,
            Helix.segment_count_folder,
            Helix.total_nt,
            Helix.global_bend_angle,
            a_h,
        ).all()
        j_res = j_q.with_entities(
            Junction.pdb_id,
            Junction.organism,
            Junction.method,
            Junction.resolution,
            Junction.segment_count_folder,
            Junction.total_nt,
            Junction.global_bend_angle,
            a_j,
        ).all()

        combined = []
        for r in h_res:
            combined.append([r[0], r[1], r[2], r[3], "HELIX", r[4], r[5], r[6], r[7]])
        for r in j_res:
            combined.append([r[0], r[1], r[2], r[3], "JUNCTION", r[4], r[5], r[6], r[7]])

        is_str_col = isinstance(a_h.type, String)
        
        def sort_key(x):
            val = x[8]
            if val is None:
                return "" if is_str_col else -1.0
            return val

        combined.sort(key=sort_key, reverse=(sort_order == "desc"))

        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(
            [
                "CIF ID",
                "Source/Molecule",
                "Method",
                "Res. (A)",
                "Nts Count",
                "Bend Angle (deg)",
                "Type",
                "Segment Count",
            ]
        )
        for r in combined:
            writer.writerow(
                [
                    r[0],
                    r[1] or "Unknown",
                    r[2] or "N/A",
                    f"{r[3]:.2f}" if r[3] is not None else "-",
                    r[6],
                    f"{r[7]:.1f}" if r[7] is not None else "-",
                    r[4],
                    r[5],
                ]
            )

        return PlainTextResponse(output.getvalue(), media_type="text/csv")


@app.post("/api/download-zip")
async def download_zip(file_paths: List[str]):
    """Creates a ZIP archive from the provided file paths and returns it as a stream."""
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
        for path in file_paths:
            try:
                full_path = resolve_output_file(path)
            except HTTPException:
                raise HTTPException(
                    status_code=400,
                    detail="Files must be inside the configured results directory.",
                )

            if full_path.is_file():
                zip_file.write(full_path, full_path.name)

    zip_buffer.seek(0)
    return StreamingResponse(
        zip_buffer,
        media_type="application/x-zip-compressed",
        headers={"Content-Disposition": "attachment; filename=RNAbridge_results.zip"},
    )


@app.get("/api/export-zip")
async def export_zip(
    segment_type: Optional[List[str]] = Query(None),
    motif_type: Optional[List[str]] = Query(None),
    method: Optional[List[str]] = Query(None),
    min_angle: Optional[float] = Query(None),
    max_angle: Optional[float] = Query(None),
    min_nt: Optional[int] = Query(None),
    max_nt: Optional[int] = Query(None),
    search_pdb: Optional[str] = Query(None),
    sequence: Optional[str] = Query(None),
    stacking_stem1: Optional[str] = Query(None),
    stacking_stem2: Optional[str] = Query(None),
    non_redundant: Optional[bool] = Query(False),
):
    """Exports all filtered results as a ZIP file."""
    with SessionLocal() as session:
        h_q, j_q = apply_filters(
            session.query(Helix),
            session.query(Junction),
            {
                "segment_type": segment_type,
                "motif_type": motif_type,
                "method": method,
                "min_angle": min_angle,
                "max_angle": max_angle,
                "min_nt": min_nt,
                "max_nt": max_nt,
                "search_pdb": search_pdb,
                "sequence": sequence,
                "stacking_stem1": stacking_stem1,
                "stacking_stem2": stacking_stem2,
                "non_redundant": non_redundant,
            },
        )

        results = (
            h_q.with_entities(Helix.path_cif, Helix.path_pml).all()
            + j_q.with_entities(Junction.path_cif, Junction.path_pml).all()
        )

        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
            for r in results:
                for path in r:
                    if path:
                        try:
                            full_path = resolve_output_file(path)
                        except HTTPException:
                            continue

                        if full_path.is_file():
                            zip_file.write(full_path, full_path.name)

        zip_buffer.seek(0)
        return StreamingResponse(
            zip_buffer,
            media_type="application/x-zip-compressed",
            headers={
                "Content-Disposition": "attachment; filename=RNAbridge_all_results.zip"
            },
        )



@app.get("/api/stats")
async def get_global_stats(
    non_redundant: bool = Query(False),
):
    """Returns statistics for the selected dataset."""
    dataset = "rnasolo" if non_redundant else "pdb"

    with SessionLocal() as session:
        h_q = session.query(Helix).filter(Helix.dataset == dataset)
        j_q = session.query(Junction).filter(Junction.dataset == dataset)

        max_nt_h = h_q.with_entities(func.max(Helix.total_nt)).scalar() or 0
        max_nt_j = j_q.with_entities(func.max(Junction.total_nt)).scalar() or 0

        h_folders = h_q.with_entities(Helix.segment_count_folder).distinct().all()
        j_folders = j_q.with_entities(Junction.segment_count_folder).distinct().all()

        all_folders = set()
        for (folder,) in h_folders + j_folders:
            if folder:
                all_folders.add(folder)

        max_way = 0
        for folder in all_folders:
            try:
                if "way" in folder:
                    val = int(folder.split("-")[0])
                    max_way = max(max_way, val)
            except Exception:
                pass

        h_methods = h_q.with_entities(Helix.method).distinct().all()
        j_methods = j_q.with_entities(Junction.method).distinct().all()
        all_methods = {
            method
            for (method,) in h_methods + j_methods
            if method
        }

        return {
            "max_nt": max(max_nt_h, max_nt_j),
            "max_way": max_way,
            "folders": sorted(all_folders),
            "methods": sorted(all_methods),
        }

@app.get("/api/ids")
async def get_all_ids(
    non_redundant: bool = Query(False),
):
    """Returns IDs from the selected dataset."""
    dataset = "rnasolo" if non_redundant else "pdb"

    with SessionLocal() as session:
        helix_ids = (
            session.query(Helix.pdb_id)
            .filter(Helix.dataset == dataset)
            .distinct()
            .all()
        )
        junction_ids = (
            session.query(Junction.pdb_id)
            .filter(Junction.dataset == dataset)
            .distinct()
            .all()
        )

        return sorted({
            row[0]
            for row in helix_ids + junction_ids
            if row[0]
        })

@app.get("/api/files/{file_path:path}")
async def get_file(file_path: str):
    path = resolve_output_file(file_path)
    if path.is_file():
        mime = (
            "application/json"
            if path.suffix.lower() == ".json"
            else "image/svg+xml"
            if path.suffix.lower() == ".svg"
            else "application/octet-stream"
        )
        return FileResponse(str(path), media_type=mime)

    raise HTTPException(status_code=404)

if os.path.exists("dist"):
    app.mount("/", StaticFiles(directory="dist", html=True), name="static")

    @app.exception_handler(404)
    async def custom_404_handler(request, __):
        return FileResponse("dist/index.html")


if __name__ == "__main__":
    import uvicorn
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "export":
        from fastapi.openapi.utils import get_openapi

        openapi_schema = get_openapi(
            title=app.title,
            version=app.version,
            openapi_version=app.openapi_version,
            description=app.description,
            routes=app.routes,
        )
        with open("openapi.json", "w") as f:
            json.dump(openapi_schema, f, indent=4)
        print("Successfully exported OpenAPI schema to openapi.json")
    else:
        uvicorn.run(app, host="0.0.0.0", port=8000)
