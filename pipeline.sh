#!/bin/bash

# --- ENVIRONMENT CONFIGURATION ---
PROJECT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$PROJECT_DIR"

# Python Path setup
if [ -d "$PROJECT_DIR/venv" ]; then
    export PATH="$PROJECT_DIR/venv/bin:$PATH"
fi
PY=python3
export X3DNA="$PROJECT_DIR/x3dna-v2.4"
export PATH="$X3DNA/bin:$PATH"

# Get parameters from config.py
get_config() {
    $PY -c "import config; print(config.$1)" | tr -d '\r' | xargs
}

input_cif_dir=$(get_config "CIF_DIR")
output_json=$(get_config "JSON_DIR")
output_helices=$(get_config "RESULTS_DIR")
MAX_BEND_ANGLE=$(get_config "MAX_BEND_ANGLE")

echo "Configuration: CIF=$input_cif_dir, JSON=$output_json, HELICES=$output_helices, ANGLE=$MAX_BEND_ANGLE"

# --- STEP 0: WAIT FOR DATABASE ---
if [[ "$DATABASE_URL" == *"postgresql"* ]]; then
    echo "--- STEP 0: Waiting for PostgreSQL... ---"
    # Extract host and port from URL
    db_host=$(echo $DATABASE_URL | sed -e 's|.*@||' -e 's|/.*||' -e 's|:.*||')
    db_port=$(echo $DATABASE_URL | sed -e 's|.*:||' -e 's|/.*||')
    [[ -z "$db_port" ]] && db_port=5432
    
    until printf "" 2>>/dev/null >>/dev/tcp/$db_host/$db_port; do
        echo "Postgres is unavailable - sleeping"
        sleep 2
    done
    echo "Postgres is up!"
fi

# --- STEP 1: PDB SYNCHRONIZATION ---
echo "--- STEP 1: Synchronizing files with PDB ---"
#$PY fetch_pdb.py
if [[ $? -ne 0 ]]; then
    echo "ERROR: PDB synchronization failed. Aborting pipeline!"
    exit 1
fi

# --- STEP 2: FILE ANALYSIS ---
echo "--- STEP 2: Processing CIF files ---"
mkdir -p "$output_json" "$output_helices"

if [[ ! -d "$input_cif_dir" ]]; then
    echo "ERROR: Input directory $input_cif_dir does not exist."
    exit 1
fi

find "$input_cif_dir" -maxdepth 1 -name "*.cif" | while read cif_file; do  
    filename=$(basename "$cif_file")
    id="${filename%.cif}" 
    json_file="${output_json}/${id}.json"
    helices_file="${output_helices}/${id}.json" 

    # Skip if final result already exists
    if [[ -f "$helices_file" ]]; then
        continue
    fi

    echo "Processing ID: $id"

    # RNA Annotation
    if [[ ! -f "$json_file" ]]; then
        echo "   -> Running annotator (with 10m timeout)..."
        timeout 600s annotator --decompose-pseudoknot-free --json "$json_file" "$cif_file"
	status=$?
        if [[ $status -ne 0 ]]; then
            echo "   -> ERROR: annotator failed for $id. Skipping."
            continue
        fi
    fi

    # Analysis and Helix Building (Combined Stages)
    echo "   -> Analyzing and building helices..."
    $PY analyze.py "$json_file" "$helices_file" "$cif_file" "$MAX_BEND_ANGLE"

    echo "--------------------------------------------------------"
done

# --- STEP 3: CATEGORIZATION AND VISUALIZATION ---
echo "--- STEP 3: Generating visualizations and categorization ---"
$PY categorize.py

# --- STEP 4: DATABASE UPDATE ---
echo "--- STEP 4: Updating database ---"
$PY database.py

echo "--- PIPELINE FINISHED: $(date) ---"
