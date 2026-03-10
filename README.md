# Epitope Pipeline

Structural bioinformatics pipeline for identifying druggable VHH epitope space on human membrane protein targets.

## Overview

This pipeline takes a list of membrane protein targets and systematically identifies surface regions suitable for VHH antibody binding. It filters for regions that are:

1. **Structurally resolved** — uses AlphaFold DB full-length predictions by default, or experimental PDB structures via opt-in
2. **Distal from the membrane** — ≥80 Angstroms from the membrane surface (default; configurable)
3. **Surface accessible** — solvent-exposed patches large enough for a VHH CDR footprint (≥600 A²)
4. **Cynomolgus-conserved** — max 2 cyno mismatches per ~600 A² VHH footprint (3D sliding window)
5. **Target-specific** — max 2 non-specific residues per ~600 A² VHH footprint (3D sliding window, ≥75% off-target BLAST identity)

## Quick Start

```python
from epitope_pipeline.run import run_pipeline

# Single target
results = run_pipeline(["ERBB2"])

# Multiple targets (gene names or UniProt IDs, can be mixed)
results = run_pipeline(["ERBB2", "EGFR", "P04626"])

# With custom thresholds
results = run_pipeline(
    ["ERBB2"],
    min_distance_a=60.0,            # Override 80A default
    cyno_max_mismatches=3,          # Override 2 mismatches per 600A² default
    specificity_threshold=0.90,     # Override 0.75 default
    force_experimental=True,        # Use experimental PDB instead of AlphaFold DB
)
```

### Command Line

```bash
python -m epitope_pipeline.run ERBB2 EGFR
```

## Pipeline Steps

| Step | Module | Description |
|------|--------|-------------|
| 1 | `target_input.py` | Resolve UniProt IDs or gene names to full protein metadata + InterPro domains |
| 2 | `structure.py` | Acquire structure: AlphaFold DB (default), experimental PDB (opt-in), or Tamarind AlphaFold (fallback) |
| 3 | `membrane.py` | Define membrane plane (OPM, TM annotations, or GPI-anchor) |
| 4 | `spatial.py` | Filter to ectodomain residues ≥80A from membrane surface (default; bispecific uses 60A/40A) |
| 5 | `surface.py` | Calculate SASA, cluster into surface patches ≥600 A² |
| 6 | `conservation.py` | Align with cyno ortholog, sliding window mismatch check + contiguous sub-patch |
| 7 | `specificity.py` | Full-sequence BLAST, sliding window non-specific density check |
| 8 | `scoring.py` | Composite score: area + distance + conservation + specificity + accessibility |
| 9 | `export.py` | Generate all output files (CSV, XLSX, PDB, FASTA, JSON) |
| 10 | `visualize.py` | 6-track epitope map with domains, distance, SASA, conservation, specificity |

## Output Files

Each run creates a date-stamped directory under `runs/`:

```
runs/2026-02-24_erbb2_egfr/
├── epitope_candidates.csv          # Main results table
├── epitope_candidates.xlsx         # Color-formatted Excel workbook
├── input_manifest.json             # Parameters for reproducibility
├── log.txt                         # Full pipeline log
├── annotated_pdbs/
│   ├── erbb2_epitope.pdb           # B-factor = epitope score (PyMOL: spectrum b)
│   └── egfr_epitope.pdb
├── annotated_sequences/
│   ├── erbb2_epitope.fasta         # FASTA with epitope annotations in header
│   ├── erbb2_residue_table.csv     # Per-residue: topology, SASA, conservation, patch
│   └── ...
├── annotations/
│   ├── erbb2_annotation.json       # Full intermediate results (JSON)
│   └── ...
├── figures/
│   ├── erbb2_epitope_map.png       # Linear epitope map
│   ├── egfr_epitope_map.png
│   └── scoring_summary.png         # Multi-target comparison
└── structures/
    ├── 1n8z.pdb                    # Downloaded/predicted structures
    └── ...
```

### Annotated PDB Scoring Scheme

The B-factor column in annotated PDBs encodes the epitope analysis:

| B-factor | Meaning |
|----------|---------|
| > 50 | Epitope patch residue (composite_score × 100) |
| 25 | Extracellular, not in qualifying patch |
| 0 | Transmembrane |
| -25 | Intracellular |

**PyMOL**: `spectrum b, blue_white_red`
**ChimeraX**: `color bfactor palette cyan:white:green`

Each annotated PDB comes with a companion `.pml` script for PyMOL with clickable named selections for each pipeline tier. Cytoplasmic residues are hidden; TM helices remain visible for structural context. For multi-pass TM proteins, the PML shows the full structure with topology-based coloring (TM helices in purple, ECD loops in teal).

### Scoring Components

| Component | Weight | Description |
|-----------|--------|-------------|
| Area | 25% | Patch size relative to VHH footprint (900 A²) |
| Distance | 20% | Distance from membrane (normalized to 300A) |
| Conservation | 25% | Cyno sequence identity within the patch |
| Specificity | 20% | Binary: 1.0 if patch passes 75% screen, 0.0 if rejected |
| Accessibility | 10% | Average relative SASA (solvent exposure) |

## Configuration

All thresholds are defined in `config.py` and can be overridden via `run_pipeline()` arguments:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ECTODOMAIN_MIN_DISTANCE_A` | 80.0 | Min Angstroms from membrane |
| `RESOLUTION_THRESHOLD_A` | 3.5 | Max PDB resolution accepted |
| `VHH_FOOTPRINT_MIN_A2` | 600.0 | Min surface patch area (A²) |
| `MAX_CYNO_MISMATCHES_PER_600A2` | 2 | Max cyno mismatches per ~600A² VHH footprint |
| `MAX_NONSPECIFIC_PER_600A2` | 2 | Max non-specific residues per ~600A² footprint |
| `SPECIFICITY_IDENTITY_THRESHOLD` | 0.75 | Off-target BLAST identity for binary screening |
| `PATCH_CLUSTERING_DISTANCE_A` | 8.0 | Cα-Cα distance for patch connectivity |
| `MIN_ECTODOMAIN_COVERAGE` | 0.80 | Min PDB ectodomain coverage (applies when `force_experimental=True`) |

## Structure Acquisition

The pipeline defaults to AlphaFold DB for full-length, gap-free structures. Experimental PDB structures are available as an opt-in override.

1. **AlphaFold DB** (default) — Pre-computed full-length structures from EBI's AlphaFold Protein Structure Database. Queries the API to discover the latest model version (currently v6). No API key required. Full-length models provide complete SASA coverage without gaps from disordered loops or uncovered regions.
2. **Experimental PDB** (opt-in via `force_experimental=True`) — RCSB search for X-ray/cryo-EM structures ≤3.5A resolution, ranked by ectodomain coverage × (1/resolution). Requires ≥80% coverage of the ectodomain (not full-length protein). PDB entity sequences are aligned to the UniProt sequence via BioPython PairwiseAligner for accurate coverage mapping.
3. **Tamarind AlphaFold** — De novo AlphaFold prediction via the Tamarind API (requires API key). Used only for proteins not in AlphaFold DB.

## Dependencies

```
# Core
biopython>=1.81          # PDB parsing, alignment, BLAST
freesasa>=2.2.0          # Solvent-accessible surface area
scipy>=1.10              # Spatial clustering
numpy>=1.24,<2.0         # Array operations
requests>=2.28           # API calls
python-dotenv>=0.21      # Environment config

# Visualization & export
matplotlib>=3.7.1,<3.8   # Figures
openpyxl>=3.0            # Excel output
```

### Installation

```bash
pip install -r requirements.txt
```

**Note**: Most human proteins have pre-computed structures in AlphaFold DB (no API key needed). The Tamarind API is only required as a last resort for proteins not in AlphaFold DB. Set your API key in `tamarind/.env`:
```
TAMARIND_API_KEY=your_key_here
```

### Local BLAST Database Setup (Recommended)

For faster specificity screening (10-20x speedup), set up a local BLAST database:

**1. Install NCBI BLAST+ tools:**
```bash
# macOS
brew install blast

# Ubuntu/Debian
sudo apt-get install ncbi-blast+
```

**2. Run the automated setup script:**
```bash
bash epitope_pipeline/scripts/setup_blast_db.sh
```

This will:
- Download the SwissProt database (~500 MB)
- Filter to human proteins only
- Build a local BLAST database (~300 MB)
- Takes ~10-15 minutes

**3. Enable local BLAST in config.py:**
```python
USE_LOCAL_BLAST = True  # Set to False to use remote NCBI BLAST API
```

**4. Verify setup:**
```bash
blastp -db epitope_pipeline/blast_db/swissprot/swissprot_human \
       -query <test.fasta> -outfmt 6
```

**Benefits:**
- **Performance**: ~100ms per query (vs 1-2 seconds with remote API)
- **Reliability**: No network dependency or NCBI service downtime
- **No rate limiting**: Can run queries in parallel

**Note**: SwissProt is updated quarterly by UniProt. Re-run the setup script periodically to refresh the database.

## Epitope Map Figure

Each target gets a 6-track linear sequence map:

| Track | Content | Details |
|-------|---------|---------|
| **Domains** | Topology brackets + domain blocks | Extracellular/Cytoplasmic brackets above; InterPro domains (teal), UniProt domains (blue), TM (purple), SP (gray), Disordered (gray) blocks |
| **Distance** | Distance from membrane (A) | Continuous line through all residues; ECD regions in teal fill, TM/cytoplasmic in light gray fill; dashed threshold line at the configured min distance |
| **SASA** | Per-residue solvent-accessible surface area (A²) | Bar chart (teal) |
| **Cyno Conservation** | Per-residue identity with cynomolgus ortholog | Mint = conserved, red = mismatch |
| **Human Specificity** | Binary per-patch off-target screening (BLAST) | Mint = specific, red = non-specific (≥75% identity off-target HSP on ≥600 A² ectodomain patch) |
| **Target Epitope** | Final qualifying patches | Green = passing all filters. Label: "Patch N — total A² \| score X.XX" |

Domain annotations are sourced from both UniProt (TM, Signal peptide, Domain, Region) and the InterPro API (ECD subdomains like Receptor L-domain, Furin-like cysteine-rich, etc.). InterPro entries that overlap significantly with UniProt features are automatically filtered out.

All font sizes across single-target and bispecific figures were calibrated for readability at 150 DPI (titles 16pt, track labels 10-11pt, tick labels 9-10pt, domain text 7-9pt, patch labels 8.5pt).

## Supported Membrane Protein Types

| Type | Detection | Membrane Plane Method |
|------|-----------|----------------------|
| **Single-pass TM** | UniProt `Transmembrane` annotation | PCA on TM helix Cα atoms, or OPM lookup |
| **Multi-pass TM** | Multiple `Transmembrane` annotations | Plane through TM helix midpoints, or OPM |
| **GPI-anchored** | UniProt `Lipidation` with "GPI" in description | PCA on ectodomain + 15A extrapolation from C-terminus |

For ectodomain-only crystal structures (no TM residues resolved), the membrane plane is estimated from the C-terminal boundary of the resolved structure with a 30A extrapolation.

### Distance Measurement

Distances are measured as the projected distance along the membrane normal from the **bilayer surface** (top of the lipid bilayer = membrane center + half-thickness), not from the membrane center. This is equivalent to "height above the cell surface."

The membrane surface is defined geometrically from the membrane half-thickness (15A default), which is more robust than using the most proximal ectodomain residue as a proxy — AlphaFold models predict in vacuum without membrane context, so juxtamembrane residues can artifactually collapse to membrane-plane height (observed for ERBB2 and EGFR).

## Cynomolgus Ortholog Resolution

The cyno conservation step requires the cynomolgus macaque ortholog sequence. Resolution follows a fallback chain:

1. **UniProt reviewed** — search `gene:{name} AND organism_id:9541 AND reviewed:true`
2. **UniProt unreviewed** — same query without `reviewed:true`
3. **Ensembl ortholog** — REST API homology lookup (`/homology/symbol/homo_sapiens/{gene}`)
4. **None** — if no ortholog found, conservation step is skipped for that target

The Ensembl fallback handles genes with divergent cyno orthologs (e.g., CEACAM5, where the CEACAM family is rearranged in macaques and the ortholog is only ~57% identical).

## API Caching

All external API calls are cached to `cache/` for efficient re-runs:

| Cache | Location | Contents |
|-------|----------|----------|
| UniProt entries | `cache/uniprot/` | Full JSON entries by accession |
| Cyno orthologs | `cache/cyno_sequences/` | Ortholog accession + sequence |
| InterPro domains | `cache/interpro/` | Domain-type annotations by accession |
| PDB metadata | `cache/pdb/` | RCSB entry metadata |
| OPM orientations | `cache/opm/` | Pre-oriented PDB files |
| BLAST results | `cache/blast/` | Per-patch + full-sequence hits keyed by sequence hash |

To force fresh API calls, delete the relevant cache files or directories.

## Architecture

```
epitope_pipeline/
├── config.py           # Central configuration
├── utils.py            # Shared helpers (extract_ca_coords, get_chain, setup_logging)
├── run.py              # Pipeline orchestrator + CLI
├── target_input.py     # Step 1: UniProt resolution
├── structure.py        # Step 2: PDB/AlphaFold acquisition
├── membrane.py         # Step 3: Membrane plane definition
├── spatial.py          # Step 4: Distance filtering
├── surface.py          # Step 5: SASA + patch clustering
├── conservation.py     # Step 6: Cyno conservation
├── specificity.py      # Step 7: BLAST specificity
├── scoring.py          # Step 8: Composite scoring
├── export.py           # Step 9: Output generation + membrane CGO
├── visualize.py        # Step 10: Figures
├── bispecific.py       # Bispecific dual-targeting orchestrator
├── export_bispecific.py    # Bispecific output generation
└── visualize_bispecific.py # Bispecific figures (stacked maps, combined orientation, pair summary)
```

Each module is independent and can be used individually:

```python
from epitope_pipeline.target_input import resolve_targets
from epitope_pipeline.membrane import annotate_membrane

targets = resolve_targets(["ERBB2"])
# ... use individual modules as needed
```

## Bispecific Dual-Targeting Mode

Evaluates pairs of membrane protein targets for complementary epitope space suitable for bispecific antibody design. One target provides a "distal" epitope (>=60A from membrane surface) and the other provides a "proximal" epitope (<=40A from membrane surface). Both orientations are tested for each pair.

### Usage

```bash
# Single pair
python -m epitope_pipeline.bispecific ERBB2:NECTIN4

# Multiple pairs
python -m epitope_pipeline.bispecific ERBB2:NECTIN4 ERBB2:MSLN
```

### How it works

1. **Shared steps**: Target resolution, structure acquisition, and membrane annotation are performed once per unique target
2. **Zone analysis**: Each target is analyzed in both distal (>=60A) and proximal (<=40A) modes — Steps 4-8 run independently per zone
3. **Pair evaluation**: Both orientations evaluated per pair (A-distal/B-proximal and B-distal/A-proximal)
4. **Scoring**: Geometric mean of best distal and proximal composite scores. 20% flexibility bonus when both orientations are valid.

### Outputs

```
runs/YYMMDD_HHMM_bispecific_erbb2_nectin/
├── bispecific_pairs.csv               # Pair scores, orientations, component scores
├── bispecific_pairs.xlsx              # Color-formatted pair scores
├── input_manifest.json
├── log.txt
├── figures/
│   ├── bispecific_combined_*.png      # Side-by-side orientation comparison (3 tracks per target)
│   ├── bispecific_*_distal__*.png     # Stacked dual-panel 6-track maps
│   ├── bispecific_pair_summary.png    # Bar chart comparing pair scores
│   └── zone_details/                  # Per-zone deep-dive maps
│       ├── erbb2_distal_epitope_map.png
│       └── ...                        # 4 zone-specific 6-track single-target maps
├── pymol/                             # All 3D artifacts
│   ├── *_distal__*_proximal.pml       # Dual PML: both structures side-by-side
│   │                                  #   with shared membrane bilayer CGO
│   └── zone_sessions/                 # Per-zone deep-dive PDB + PML pairs
│       ├── erbb2_distal_epitope.pdb
│       ├── erbb2_distal_epitope.pml
│       └── ...
└── structures/
```

### Dual PML Alignment

The dual PML scripts place two structures side-by-side on a shared membrane bilayer. Three alignment steps ensure a clean presentation regardless of membrane normal orientation:

1. **Membrane-aligned view**: A `set_view` rotation matrix is computed so screen-Y = membrane normal (membrane horizontal) and screen-X = side-by-side direction. An ectodomain-up check flips the view if needed.
2. **Screen-right translation**: The proximal structure is translated along the `right` vector (perpendicular to the membrane normal in the viewing plane), not along model-X. This prevents vertical offset when the membrane normal is tilted relative to the model axes.
3. **Vertical membrane alignment**: Both proteins' membrane centers are aligned along the normal axis so their TM regions sit at the same height in the shared bilayer. Without this, averaging two different membrane centers would bury one protein.

Epitope patches are shown as green surfaces (`show surface`) with 30% transparency. **Note**: PyMOL `ray_trace_mode 3` cannot render molecular surfaces — use `ray_trace_mode 0` or `1` for ray-traced images.

### Bispecific Test Set

| Pair | Score | Validity | Best Orientation | Notes |
|------|-------|----------|------------------|-------|
| ERBB2 x NECTIN4 | 1.016 | DUAL-VALID | NECTIN4(distal) x ERBB2(proximal) | Gold standard. ERBB2 prox trimmed: 83-res sub-patch salvaged from 129-res parent |
| MSLN x ADAM8 | 0.833 | Single | MSLN(distal) x ADAM8(proximal) | Only one valid orientation |
| ITGB6 x ADAM8 | 0.814 | Single | ITGB6(distal) x ADAM8(proximal) | Only one valid orientation |
| CEACAM5 x MUCL3 | 0.000 | None | — | Cyno + specificity wipeout, no viable patches |

## Test Set (22 targets)

### Targets with qualifying epitope patches

| Target | Gene | UniProt | Type | Structure | Patches | Score | Notes |
|--------|------|---------|------|-----------|---------|-------|-------|
| EGFR | EGFR | P00533 | Single-pass | 1YY9 (X-ray 2.6A) | 1 | 0.835 | Tested at 60A override |
| NECTIN4 | NECTIN4 | Q96NY8 | Single-pass | AlphaFold DB v6 | 2 | 0.818 | 99.4% cyno, 100% patch conservation |
| GPNMB | GPNMB | Q14956 | Single-pass | AlphaFold DB v6 | 3 | 0.803 | 1/4 patches rejected (85% cyno) |
| ITGAV | ITGAV | P06756 | Single-pass | 3IJE (X-ray 2.9A) | 6 | 0.788 | 99.1% cyno, all specific |
| FOLH1 | FOLH1 | Q04609 | Single-pass | 3BI1 (X-ray 1.5A) | 1 | 0.653 | 97.6% overall cyno, patch 100% |

### Targets with 0 qualifying patches

| Target | Gene | UniProt | Type | Structure | Rejection Reason |
|--------|------|---------|------|-----------|------------------|
| HER2 | ERBB2 | P04626 | Single-pass | AlphaFold DB v6 | Needs re-test with sliding window threshold (was rejected at 97.2% under old 98% threshold) |
| ITGB6 | ITGB6 | P18564 | Single-pass | 4UM8 (X-ray 2.9A) | 33 res ≥80A, too few for ≥600A² patch |
| ADAM8 | ADAM8 | P78325 | Single-pass | AlphaFold DB v6 | 4DD8 33% ECD → AF; max 79.5A |
| MSLN | MSLN | Q13421 | GPI-anchored | AlphaFold DB v6 | 7UED 55% ECD → AF; cyno 85.1% |
| MELTF | MELTF | P08582 | GPI-anchored | 6XR0 (X-ray 3.1A) | Max projected dist 80.2A (GPI PCA) |
| CEA | CEACAM5 | P06731 | GPI-anchored | AlphaFold DB v6 | 8BW0 28% ECD → AF; cyno 62%, CEACAM cross-reactive |
| CLDN18 | CLDN18 | P56856 | Multi-pass (4TM) | AlphaFold DB v6 | ECD loops only 36.6A |
| CD19 | CD19 | P15391 | Single-pass | 6AL5 (X-ray 3.0A) | Patch cyno 82.4% |
| CD20 | MS4A1 | P11836 | Multi-pass (4TM) | 8VGN (cryo-EM 2.5A) | Max 15.2A (short ECD loops) |
| PD-L1 | CD274 | Q9NZQ7 | Single-pass | AlphaFold DB v6 | 5O45 53% ECD → AF; patch cyno 93.9% |
| TROP2 | TACSTD2 | P09758 | Single-pass | 7PEE (X-ray 2.8A) | Max 62.9A |
| HER3 | ERBB3 | P21860 | Single-pass | 1M6B (X-ray 2.6A) | Needs re-test (was rejected at 96.2% under old threshold) |
| DLL3 | DLL3 | Q9NYJ7 | Single-pass | AlphaFold DB v6 | Max 67.1A |
| CLDN6 | CLDN6 | P56747 | Multi-pass (4TM) | AlphaFold DB v6 | Max 34.8A (short ECD loops) |
| FLT3 | FLT3 | P36888 | Single-pass | AlphaFold DB v6 | 1RJB 0% ECD → AF; needs re-test (was rejected at 97.9% under old threshold) |
| BCMA | TNFRSF17 | Q02223 | Single-pass | 4ZFO (X-ray 1.9A) | Max 48.2A (tiny 54-res ECD) |
| GPC3 | GPC3 | P51654 | GPI-anchored | 7ZAW (X-ray 2.6A) | Max 69.6A |

## Known Limitations

| Limitation | Affected Targets | Details |
|------------|-----------------|---------|
| **PCA membrane normal for GPI proteins** | MELTF, GPC3 | For multi-lobe GPI-anchored proteins, PCA principal axis doesn't align with the max protein dimension. MELTF: 80A projected vs 107A Euclidean. GPC3: 69.6A projected. |
| **Multi-pass TM — tiny ECD** | CLDN18, CLDN6, MS4A1 | Proteins with 4+ TM helices have short extracellular loops (CLDN18: 36.6A, CLDN6: 34.8A, MS4A1: 15.2A). Always 0 patches at 80A. |
| **Monomer-only analysis** | ITGB6/ITGAV | Heterodimer interfaces (e.g., integrin αVβ6) are not captured — each chain is analyzed independently. |
| **Low cyno conservation** | MSLN, CD19, CD274, CEACAM5 | Multiple ADC targets have clustered cyno mismatches that fail the sliding window check. Borderline targets (ERBB2, ERBB3, FLT3) need re-test with new threshold. |
| **Small ectodomains** | TNFRSF17, TACSTD2, DLL3 | Ectodomains too small to reach 80A: BCMA (48A, 54 res), TROP2 (63A), DLL3 (67A). |
