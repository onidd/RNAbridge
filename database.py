import os
import json
import re
from pathlib import Path
from sqlalchemy import create_engine, Column, Integer, String, Float, ForeignKey, Text
from sqlalchemy.orm import declarative_base, sessionmaker, relationship
from rnabridge import Utils, MetadataProvider
from config import CIF_DIR, JSON_DIR, OUTPUT_DIR, DATABASE_URL, BASE_DIR

Base = declarative_base()

CACHE_FILE = BASE_DIR / "metadata_cache.json"

def load_metadata_cache():
    if os.path.exists(CACHE_FILE):
        try:
            with open(CACHE_FILE, "r") as f:
                return json.load(f)
        except Exception:
            return {}
    return {}

def save_metadata_cache(cache):
    with open(CACHE_FILE, "w") as f:
        json.dump(cache, f, indent=4)


class Helix(Base):
    """Database model for RNA Super-Helices."""

    __tablename__ = "helices"
    id = Column(Integer, primary_key=True, autoincrement=True)
    pdb_id = Column(String, index=True)
    dataset = Column(String, index=True, nullable=False, default="pdb")
    molecule = Column(String)
    organism = Column(String)
    resolution = Column(Float)
    method = Column(String)
    segment_count_folder = Column(String)
    chains = Column(String, index=True)
    total_nt = Column(Integer)
    global_bend_angle = Column(Float)
    path_json = Column(String, index=True, unique=True)
    path_cif = Column(String)
    path_pml = Column(String)
    path_svg = Column(String)
    sequence_full = Column(Text)
    segments = relationship(
        "Segment", back_populates="helix", cascade="all, delete-orphan"
    )


class Segment(Base):
    """Database model for components (motifs/stems) within a helix."""

    __tablename__ = "segments"
    id = Column(Integer, primary_key=True, autoincrement=True)
    helix_id = Column(Integer, ForeignKey("helices.id"), index=True)
    s_id = Column(Integer)
    type = Column(String)
    size_2d_left = Column(Integer)
    size_2d_right = Column(Integer)
    size_3d_left = Column(Integer)
    size_3d_right = Column(Integer)
    bend_angle = Column(Float)
    stacking = Column(String)
    stacking_path = Column(Text)
    sequence_full = Column(Text)
    helix = relationship("Helix", back_populates="segments")


class Junction(Base):
    """Database model for RNA Junctions."""

    __tablename__ = "junctions"
    id = Column(Integer, primary_key=True, autoincrement=True)
    pdb_id = Column(String, index=True)
    dataset = Column(String, index=True, nullable=False, default="pdb")
    molecule = Column(String)
    organism = Column(String)
    resolution = Column(Float)
    method = Column(String)
    segment_count_folder = Column(String)
    chains = Column(String, index=True)
    total_nt = Column(Integer)
    extended_helices_count = Column(Integer, default=0)
    j_id = Column(Integer)
    stacking_status = Column(String)
    coaxial_pairs = Column(String)
    bend_angles = Column(String)
    stacking_paths = Column(Text)
    global_bend_angle = Column(Float)
    sequence = Column(Text)
    path_json = Column(String, index=True, unique=True)
    path_cif = Column(String)
    path_pml = Column(String)
    path_svg = Column(String)


def get_engine():
    """Initializes the database engine using the URL from config."""
    return (
        create_engine(DATABASE_URL, pool_pre_ping=True)
        if DATABASE_URL.startswith("postgresql")
        else create_engine(DATABASE_URL)
    )


def setup_database(engine):
    """Creates tables and performs necessary migrations."""
    from sqlalchemy import inspect, text

    Base.metadata.create_all(engine)

    inspector = inspect(engine)
    cols_h = {c["name"] for c in inspector.get_columns("helices")}
    cols_j = {c["name"] for c in inspector.get_columns("junctions")}

    with engine.begin() as conn:
        if "dataset" not in cols_h:
            conn.execute(
                text(
                    "ALTER TABLE helices "
                    "ADD COLUMN dataset VARCHAR NOT NULL DEFAULT 'pdb'"
                )
            )

        if "dataset" not in cols_j:
            conn.execute(
                text(
                    "ALTER TABLE junctions "
                    "ADD COLUMN dataset VARCHAR NOT NULL DEFAULT 'pdb'"
                )
            )

        if "method" not in cols_h:
            conn.execute(text("ALTER TABLE helices ADD COLUMN method VARCHAR"))
        if "chains" not in cols_h:
            conn.execute(text("ALTER TABLE helices ADD COLUMN chains VARCHAR"))

        if "method" not in cols_j:
            conn.execute(text("ALTER TABLE junctions ADD COLUMN method VARCHAR"))
        if "chains" not in cols_j:
            conn.execute(text("ALTER TABLE junctions ADD COLUMN chains VARCHAR"))
        if "global_bend_angle" not in cols_j:
            conn.execute(
                text("ALTER TABLE junctions ADD COLUMN global_bend_angle FLOAT")
            )
        if "stacking_paths" not in cols_j:
            conn.execute(text("ALTER TABLE junctions ADD COLUMN stacking_paths TEXT"))

    return sessionmaker(bind=engine)()

def get_auth_label_from_cache(bp_cache, pdb_id, ser):
    """Converts a raw serial number to its PDB auth-number label (e.g. 'A1005')."""
    entry = bp_cache.get(pdb_id, {}).get(str(ser))
    if entry and entry.get("auth"):
        auth = entry["auth"]
        num, icode = (
            auth.get("number", ""),
            auth.get("icode") or auth.get("insertion_code"),
        )
        icode_str = (
            icode if icode and str(icode).strip() not in [".", "?", ""] else ""
        )
        return f"{auth.get('chain', '')}{num}{icode_str}"
    return str(ser)


def convert_junction_stacking_paths(raw_paths, bp_cache, pdb_id):
    """Converts every serial number inside a junction's stacking_paths dict to
    its auth label, matching what's shown on the VARNA 2D rendering."""
    converted = {}
    for k, paths in raw_paths.items():
        if paths and isinstance(paths[0], list):
            converted[k] = [
                [get_auth_label_from_cache(bp_cache, pdb_id, s) for s in p]
                for p in paths
            ]
        else:
            converted[k] = [get_auth_label_from_cache(bp_cache, pdb_id, s) for s in paths]
    return converted



def add_helix_segments(session, h_obj, h_data, pdb_id, bp_cache, bpseq_key=None):
    """Helper to populate the segments table for a specific helix."""
    cache_key = bpseq_key or pdb_id

    if cache_key not in bp_cache:
        orig = Path(JSON_DIR) / f"{cache_key}.json"
        if not orig.exists() and cache_key != pdb_id:
            orig = Path(JSON_DIR) / f"{pdb_id}.json"

        bp_cache[cache_key] = (
            json.load(open(orig, "r")).get("bpseq_index", {}) if orig.exists() else {}
        )

    def get_auth_label(ser):
        entry = bp_cache[cache_key].get(str(ser))        
        if entry and entry.get("auth"):
            auth = entry["auth"]
            num, icode = (
                auth.get("number", ""),
                auth.get("icode") or auth.get("insertion_code"),
            )
            icode_str = (
                icode if icode and str(icode).strip() not in [".", "?", ""] else ""
            )
            return f"{auth.get('chain', '')}{num}{icode_str}"
        return str(ser)

    gid = 1
    up = h_data.get("strands", {}).get("upstream", {})
    if up:
        s5, s3 = up.get("strand5p", {}), up.get("strand3p", {})
        session.add(
            Segment(
                helix_id=h_obj.id,
                s_id=gid,
                type="STEM",
                size_2d_left=s5.get("structure", "").count("("),
                size_2d_right=s3.get("structure", "").count(")"),
                size_3d_left=s5.get("structure", "").count("("),
                size_3d_right=s3.get("structure", "").count(")"),
                bend_angle=0.0,
                stacking="FULL",
                stacking_path="[]",
                sequence_full=f"{s5.get('sequence', '')}-{s3.get('sequence', '')}",
            )
        )
        gid += 1

    for comp in h_data.get("components", []):
        m = comp.get("metrics", {})
        s2, s3d = m.get("size_2d", [0, 0]), m.get("size_3d", [0, 0])
        rsp = (
            [[get_auth_label(x) for x in p] for p in m.get("stacking_path", [])]
            if m.get("stacking_path") and isinstance(m.get("stacking_path")[0], list)
            else [get_auth_label(x) for x in m.get("stacking_path", [])]
        )
        c_seq = "-".join(
            [
                loc.get("sequence", "")
                for loc in comp.get("location", [])
                if loc.get("sequence")
            ]
        )
        session.add(
            Segment(
                helix_id=h_obj.id,
                s_id=gid,
                type=comp.get("type", "").upper(),
                size_2d_left=s2[0] if isinstance(s2, list) else 0,
                size_2d_right=s2[1] if isinstance(s2, list) and len(s2) > 1 else 0,
                size_3d_left=s3d[0] if isinstance(s3d, list) else 0,
                size_3d_right=s3d[1] if isinstance(s3d, list) and len(s3d) > 1 else 0,
                bend_angle=m.get("bend_angle"),
                stacking=m.get("stacking"),
                stacking_path=json.dumps(rsp),
                sequence_full=c_seq,
            )
        )
        gid += 1

        i_stem = comp.get("internal_stem", {})
        if i_stem and (
            i_stem.get("strand5p", {}).get("sequence")
            or i_stem.get("strand3p", {}).get("sequence")
        ):
            is5, is3 = i_stem.get("strand5p", {}), i_stem.get("strand3p", {})
            session.add(
                Segment(
                    helix_id=h_obj.id,
                    s_id=gid,
                    type="STEM",
                    size_2d_left=is5.get("structure", "").count("("),
                    size_2d_right=is3.get("structure", "").count(")"),
                    size_3d_left=is5.get("structure", "").count("("),
                    size_3d_right=is3.get("structure", "").count(")"),
                    bend_angle=0.0,
                    stacking="FULL",
                    stacking_path="[]",
                    sequence_full=f"{is5.get('sequence', '')}-{is3.get('sequence', '')}",
                )
            )
            gid += 1


def extract_motif_chains(data_obj: dict) -> str:
    """Extracts sorted unique chain IDs (case-sensitive, e.g. 'c' or 'A-B') from a motif JSON object."""
    chains = set()

    def add_s(s):
        if isinstance(s, dict):
            c = s.get("first", {}).get("chain") or s.get("last", {}).get("chain")
            if c is not None and str(c).strip() != "":
                chains.add(str(c).strip())

    for loc in data_obj.get("location", []):
        add_s(loc)
    for stem in data_obj.get("context", {}).values():
        if isinstance(stem, dict):
            add_s(stem.get("strand5p"))
            add_s(stem.get("strand3p"))
    for s in ["upstream", "downstream"]:
        add_s(data_obj.get("strands", {}).get(s, {}).get("strand5p"))
        add_s(data_obj.get("strands", {}).get(s, {}).get("strand3p"))
    for c in data_obj.get("components", []):
        for loc in c.get("location", []):
            add_s(loc)
        add_s(c.get("internal_stem", {}).get("strand5p"))
        add_s(c.get("internal_stem", {}).get("strand3p"))

    return "-".join(sorted(chains))

def update_database(session, root_dir):
    """Scans classification_result and updates the database with new entries."""
    print("--- STARTING DATABASE UPDATE ---")

    def to_rel(p):
        try:
            return os.path.relpath(str(p), os.getcwd())
        except Exception:
            return str(p)

    svg_cache = {}
    for f in Path(root_dir).glob("**/*.svg"):
        key = (
            f.name.replace("varna-tz-", "")
            .replace("-clean.svg", "")
            .replace("_varna", "")
        )
        svg_cache[key] = to_rel(f)

    json_files = [
        f
        for f in Path(root_dir).glob("**/*.json")
        if not f.name.endswith("_varna.json")
    ]

    # Reload existing with variants
    existing_h, existing_j = {}, {}
    for h in session.query(Helix).all():
        existing_h[h.path_json] = h
        existing_h[os.path.abspath(h.path_json)] = h
    for j in session.query(Junction).all():
        existing_j[j.path_json] = j
        existing_j[os.path.abspath(j.path_json)] = j

    metadata_cache = load_metadata_cache()
    bpseq_cache, batch_size = {}, 200
    cache_modified = False

    for i in range(0, len(json_files), batch_size):
        batch = json_files[i : i + batch_size]
        for path in batch:
            ps_rel = to_rel(path)
            ps_abs = str(path.absolute())
            name_parts = path.stem.split("_")
            if len(name_parts) >= 2 and name_parts[0].upper() == "PDB":
                pdb_id = name_parts[1][-4:].upper()
            else:
                pdb_id = name_parts[0].upper()

            dataset_name = os.environ.get("RNABRIDGE_DATASET", "pdb").lower()

            bpseq_key = pdb_id
            if dataset_name == "rnasolo":
                match = re.match(
                    r"^(.*)_(?:h\d+_\d+-segment|j\d+_\d+way)$",
                    path.stem,
                )
                if match:
                    bpseq_key = match.group(1)

            if bpseq_key not in bpseq_cache:
                orig = Path(JSON_DIR) / f"{bpseq_key}.json"
                if not orig.exists() and bpseq_key != pdb_id:
                    orig = Path(JSON_DIR) / f"{pdb_id}.json"

                try:
                    bpseq_cache[bpseq_key] = (
                        json.load(open(orig, "r", encoding="utf-8")).get("bpseq_index", {})
                        if orig.exists()
                        else {}
                    )
                except Exception:
                    bpseq_cache[bpseq_key] = {}

            svg = svg_cache.get(path.stem)

            # Check if record exists but is missing SVG
            target = (
                existing_h.get(ps_rel)
                or existing_h.get(ps_abs)
                or existing_j.get(ps_rel)
                or existing_j.get(ps_abs)
            )
            if target:
              if not target.path_svg and svg:
                  target.path_svg = svg
                  session.flush()
              if not target.chains:
                  try:
                      d_tmp = json.load(open(path, "r", encoding="utf-8"))
                      items = d_tmp.get("helices") or d_tmp.get("junctions") or []
                      if items:
                          target.chains = extract_motif_chains(items[0])
                          session.flush()
                  except Exception:
                      pass
              continue

            if pdb_id not in metadata_cache:
                print(f"   -> Fetching metadata for {pdb_id} from PDB...")
                metadata_cache[pdb_id] = MetadataProvider.fetch_metadata(pdb_id)
                cache_modified = True
            
            mol, org, res, met = metadata_cache[pdb_id]

            try:
                data = json.load(open(path, "r", encoding="utf-8"))
            except Exception:
                continue

            if "helices" in data:
                for h_d in data["helices"]:
                    h = Helix(
                        pdb_id=pdb_id,
                        dataset=dataset_name,
                        molecule=mol,
                        organism=org,
                        resolution=res,
                        method=met,
                        segment_count_folder=path.parent.name,
                        chains=extract_motif_chains(h_d),
                        total_nt=h_d.get("total_nt"),
                        global_bend_angle=h_d.get("global_measures", {}).get("angle") if path.parent.name != '1-segment-helis' else (h_d.get("global_measures", {}).get("angle") or 0.0),
                        path_json=ps_rel,
                        path_cif=to_rel(path.with_suffix(".cif")),
                        path_pml=to_rel(path.with_suffix(".pml")),
                        path_svg=svg,
                        sequence_full=Utils.get_full_sequence_from_data(h_d),
                    )
                    session.add(h)
                    session.flush()
                    add_helix_segments(session, h, h_d, pdb_id, bpseq_cache, bpseq_key)
            elif "junctions" in data:
                for j_d in data["junctions"]:
                    angs = j_d.get("modules", {}).get("geometry", {}).get("bend_angles", {})
                    coaxial_pairs = j_d.get("modules", {}).get("stacking", {}).get("coaxial_pairs", [])
                    
                    coaxial_angles = []
                    for p in coaxial_pairs:
                        if len(p) >= 2:
                            s1, s2 = str(p[0]), str(p[1])
                            keys_to_check = [f"{s1}_{s2}", f"{s2}_{s1}", f"stem_{s1}_stem_{s2}", f"stem_{s2}_stem_{s1}"]
                            for k in keys_to_check:
                                if k in angs and angs[k] is not None:
                                    coaxial_angles.append(angs[k])
                                    break
                    best_angle = min(coaxial_angles) if coaxial_angles else None

                    conv_stacking_paths = convert_junction_stacking_paths(
                        j_d.get("modules", {}).get("stacking", {}).get("stacking_paths", {}),
                        bpseq_cache,
                        bpseq_key,
                    )
                    j = Junction(                        
                        pdb_id=pdb_id,
                        dataset=dataset_name,
                        molecule=mol,
                        organism=org,
                        resolution=res,
                        method=met,
                        segment_count_folder=path.parent.name,
                        chains=extract_motif_chains(j_d),
                        total_nt=j_d.get("total_nt"),
                        j_id=j_d.get("j_id"),
                        stacking_status=j_d.get("modules", {})
                        .get("stacking", {})
                        .get("status"),
                        coaxial_pairs=json.dumps(
                            j_d.get("modules", {})
                            .get("stacking", {})
                            .get("coaxial_pairs", [])
                        ),
                        bend_angles=json.dumps(angs),
                        stacking_paths=json.dumps(conv_stacking_paths),
                        global_bend_angle=best_angle,    
                        sequence=Utils.get_full_sequence_from_data(j_d),
                        path_json=ps_rel,
                        path_cif=to_rel(path.with_suffix(".cif")),
                        path_pml=to_rel(path.with_suffix(".pml")),
                        path_svg=svg,
                    )
                    session.add(j)
        session.commit()
        if cache_modified:
            save_metadata_cache(metadata_cache)
            cache_modified = False
        print(
            f"Database progress: {min(i + batch_size, len(json_files))}/{len(json_files)}"
        )


def cleanup_database(session):
    """Removes records for CIF files that no longer exist locally."""
    cif_p = Path(CIF_DIR)
    if not cif_p.exists():
        return
    active = {f.stem.upper() for f in cif_p.glob("*.cif*")}

    if len(active) < 5000:
        print(f"[SAFETY GUARD] Only {len(active)} CIF files found in {cif_p}. Skipping database cleanup!")
        return

    td_h = [
        row.id
        for row in session.query(Helix.id, Helix.pdb_id)
        .filter(Helix.dataset == "pdb")
        .all()
        if row.pdb_id.upper() not in active
    ]
    td_j = [
        row.id
        for row in session.query(Junction.id, Junction.pdb_id)
        .filter(Junction.dataset == "pdb")
        .all()
        if row.pdb_id.upper() not in active
    ]

    if len(td_h) + len(td_j) > 500:
        print(
            f"[SAFETY GUARD] Too many database records marked for deletion "
            f"(Helices: {len(td_h)}, Junctions: {len(td_j)}). Skipping database cleanup!"
        )
        return

    if td_h:
        session.query(Segment).filter(Segment.helix_id.in_(td_h)).delete(
            synchronize_session=False
        )
        session.query(Helix).filter(Helix.id.in_(td_h)).delete(
            synchronize_session=False
        )
    if td_j:
        session.query(Junction).filter(Junction.id.in_(td_j)).delete(
            synchronize_session=False
        )
    session.commit()

if __name__ == "__main__":
    eng = get_engine()
    sess = setup_database(eng)
    dataset_name = os.environ.get("RNABRIDGE_DATASET", "pdb").lower()
    if dataset_name == "pdb":
        cleanup_database(sess)
    else:
        print(f"Skipping CIF cleanup for dataset: {dataset_name}")

    update_database(sess, OUTPUT_DIR)
