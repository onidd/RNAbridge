# RNAbridge 🧬

RNAbridge is a professional platform for the analysis of RNA structural motifs. It automates the extraction, categorization, and visualization of complex super-helices and junctions from CIF data, providing a searchable database with interactive 2D and 3D perspectives.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22693019.svg)](https://doi.org/10.5281/zenodo.22693019)
[![Python](https://img.shields.io/badge/Python-3.12+-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![FastAPI](https://img.shields.io/badge/API-FastAPI-009688?style=flat-square&logo=fastapi&logoColor=white)](https://fastapi.tiangolo.com/)
[![React](https://img.shields.io/badge/Frontend-React-61DAFB?style=flat-square&logo=react&logoColor=black)](https://reactjs.org/)
[![Docker](https://img.shields.io/badge/Docker-Ready-2496ED?style=flat-square&logo=docker&logoColor=white)](https://www.docker.com/)

---

## 🖼️ Visual Overview

| Integrated Structural Preview | Statistical Insights |
| :---: | :---: |
| ![Combined Preview](frontend/public/guide/preview_combined.png) | ![Stats Histogram](frontend/public/guide/stats_histogram.png) |
| **Synchronized 2D & 3D Viewer** | **Automated Geometry Statistics** |

---

## ✨ Key Features

- **Automated Pipeline**: Single-pass analysis from raw CIF files to a structured database.
- **Advanced Geometry**: Detailed calculation of bend angles, stacking interactions, and coaxial pairings.
- **Rich Visuals**: Automatic generation of SVG diagrams (via VARNA) and PyMOL scripts.
- **Modern Interface**: React-based dashboard with Mol* integration for seamless 3D viewing.
- **Flexible Storage**: Native support for PostgreSQL (production) and SQLite (portable).

---

## 🚀 Quick Start (Docker - Recommended)

The easiest way to run the full stack (Backend, Frontend, and Database) with all system dependencies pre-configured.

### 1. Prepare X3DNA
Due to licensing, you must provide your own X3DNA distribution:
1. Download **X3DNA v2.4** from [x3dna.org](http://x3dna.org/).
2. Extract it into the `x3dna-v2.4/` directory in the project root.

### 2. Launch
```bash
docker-compose up --build
```
- **Web Interface**: [http://localhost:8000](http://localhost:8000)
- **API Documentation**: [http://localhost:8000/docs](http://localhost:8000/docs)

---

## 🛠️ Installation & Local Development

If you prefer to run RNAbridge locally (e.g., for development), ensure the following prerequisites are met:

### Prerequisites

| Tool | Purpose | Note |
| :--- | :--- | :--- |
| **Python 3.12+** | Main Logic | Required |
| **Docker** | Visualizations | Required for `cli2rest-bio` (VARNA) |
| **X3DNA (v2.4)** | Helix Analysis | Must be in `x3dna-v2.4/` |
| **Node.js** | Frontend | Only for `npm run dev` |
| **PyMOL** | 3D Processing | Recommended |

### Setup Options

#### Option A: Conda (Recommended for Dev)
```bash
conda env create -f environment.yml
conda activate rnabridge
# Build frontend once
cd frontend && npm install && npm run build && cd ..
# Start server
python api.py
```

#### Option B: Pip
```bash
pip install -r requirements.txt
python api.py
```

---

## 📖 Usage Guide

### 1. Run Analysis Pipeline
Place your `.cif` files in the `cif/` directory, then trigger the processing:

**Using Docker:**
```bash
docker-compose exec app bash pipeline.sh
```

**Using Local Python:**
```bash
./pipeline.sh
```

### 2. Accessing the Interface
- **Web Interface**: [http://localhost:8000](http://localhost:8000)
- **API Documentation**: 
  - **Swagger UI**: [http://localhost:8000/docs](http://localhost:8000/docs)
  - **ReDoc**: [http://localhost:8000/redoc](http://localhost:8000/redoc)

- **Development Mode (Local)**: 
  - **Backend**: `python api.py` (serves API at `:8000`)
  - **Frontend (Live Reload)**: `cd frontend && npm run dev` (available at `:5173`)

---

## ⚙️ Configuration

RNAbridge uses a centralized `config.py`. Key environment variables:

- `DATABASE_URL`: Connection string. Defaults to SQLite (`sqlite:///rnabridge.db`).
  - **PostgreSQL Example**: `export DATABASE_URL="postgresql://user:password@localhost:5432/rnabridge"`
- `X3DNA`: Path to X3DNA installation (default: `./x3dna-v2.4`).

---

## 📄 Data Sources & Licensing

- **Project License**: This project is licensed under the [MIT License](LICENSE).
- **X3DNA**: Users must obtain their own license from [x3dna.org](http://x3dna.org/).
- **RCSB PDB**: Structural data is fetched from the [Protein Data Bank](https://www.rcsb.org/).

## 👥 Authors

- **Damian Zakrzewski**
- **Tomasz Zok**
- **Maciej Antczak**
