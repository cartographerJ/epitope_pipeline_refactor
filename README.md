# Epitope Pipeline

Structural bioinformatics pipeline for identifying druggable VHH epitope space on human membrane protein targets.

## Overview

This pipeline takes a list of membrane protein targets and systematically identifies surface regions suitable for VHH antibody binding. It filters for regions that are:

1. **Structurally resolved** — uses AlphaFold DB full-length predictions by default, or experimental PDB structures via opt-in
2. **Distal from the membrane** — ≥80 Angstroms from the membrane surface (default; configurable)
3. **Surface accessible** — solvent-exposed patches large enough for a VHH CDR footprint (≥600 A²)
4. **Cynomolgus-conserved** — cyno mismatch threshold scales with patch size: `min(base% * sqrt(n_residues/20), 30%)`. Default base 15% → a 20-residue patch allows 15% (3 mismatches), a 120-residue patch allows up to 30% (36 mismatches)
5. **Target-specific** — non-specific residue threshold scales identically with patch size (same sqrt formula, 30% ceiling). ≥70% off-target BLAST identity defines non-specificity

## Quick Start

```python
from epitope_pipeline.run import run_pipeline

# Single target
results = run_pipeline(["ERBB2"])

# Multiple targets (gene names or UniProt IDs, can be mixed)
results = run_pipeline(["ERBB2", "EGFR", "P04626"])

# Isoform support (GENE.X notation or UniProt isoform IDs)
results = run_pipeline(["CLDN18.2", "P56856-2"])  # Both formats supported

# Proximal mode — find epitopes close to membrane
results = run_pipeline(
    ["KISS1R"],
    max_distance_a=40.0,                # Epitopes within 40A of membrane
)

# Distal mode with custom thresholds
results = run_pipeline(
    ["ERBB2"],
    min_distance_a=60.0,                # Override 80A default
    max_cyno_mismatch_percent=20.0,     # Override 15% default (more permissive)
    max_nonspecific_percent=20.0,       # Override 15% default (more permissive)
    specificity_threshold=0.75,         # Override 0.70 default (stricter off-target identity)
    force_experimental=True,            # Use experimental PDB instead of AlphaFold DB
)
```

### Command Line

```bash
python -m epitope_pipeline.run ERBB2 EGFR
```

### Web Interface

```bash
streamlit run epitope_pipeline/app.py
```

Interactive Streamlit app supporting single-target and bispecific modes with configurable distance thresholds, cyno mismatch %, and non-specific %. Results are displayed inline with epitope maps and downloadable outputs.

## Pipeline Steps

| Step | Module | Description |
|------|--------|-------------|
| 1 | `target_input.py` | Resolve UniProt IDs or gene names to full protein metadata + InterPro domains |
| 2 | `structure.py` | Acquire structure: AlphaFold DB (default), experimental PDB (opt-in), or Tamarind AlphaFold (fallback) |
| 3 | `membrane.py` | Define membrane plane (OPM, TM annotations, or GPI-anchor) |
| 4 | `spatial.py` | Filter to ectodomain residues ≥80A from membrane surface (default; bispecific uses 60A/40A) |
| 5 | `surface.py` | Calculate SASA, cluster into surface patches ≥600 A² |
| 6 | `conservation.py` | Align with cyno ortholog, size-scaled patch evaluation (base 15%, sqrt scaling, 30% cap) |
| 7 | `specificity.py` | Full-sequence BLAST, size-scaled patch evaluation (same formula) + merge adjacent patches |
| 8 | `scoring.py` | Composite score: area (60%) + conservation (25%) + specificity (15%) |
| 9 | `export.py` | Generate all output files (CSV, XLSX, PDB, FASTA, JSON) |
| 10 | `visualize.py` | 6-track epitope map with domains, distance, SASA, conservation, specificity |

## Filtering Strategy: Whole-Patch Evaluation

**March 2026 Update**: The pipeline uses a fundamentally different approach for conservation and specificity filtering compared to earlier versions.

### Philosophy: VHH Binding Units

VHH antibodies bind contiguous ~600 A² surface patches as atomic units. The CDR loops contact ~15-25 residues simultaneously, and not all contacted residues are equally critical for binding. This biological reality informs our filtering strategy:

**Whole-patch evaluation** (current approach):
- Each ≥600 A² patch is treated as an indivisible VHH binding surface
- Accept/reject based on **percentage of problematic residues**
- Simple pass/fail: "This patch has 8% cyno mismatches → PASS (under threshold)" — threshold scales with patch size via `min(base% * sqrt(n/20), 30%)`
- Patches remain intact throughout filtering — no trimming, no size changes
- Adjacent patches that pass filtering are merged into larger contiguous epitope regions (15Å centroid distance threshold)

**Sliding window trim-and-recluster** (deprecated, pre-March 2026):
- Scanned each residue's ~600 A² neighborhood for mismatches
- If any window exceeded threshold, trimmed failing residues and re-clustered survivors
- Could fragment a 149-residue patch into a 19-residue sub-patch
- Biologically questionable: a VHH that binds the 149-residue surface won't bind the 19-residue fragment

### Default Thresholds

- **Cyno conservation**: size-scaled threshold `min(base% * sqrt(n/20), 30%)`, base=`MAX_CYNO_MISMATCH_PERCENT` (default 15%)
- **Human specificity (paralog homology)**: flat threshold — `MAX_NONSPECIFIC_PERCENT` (default 85%, i.e. up to 85% of patch residues may be shared with off-target paralogs; min 15% unique)
- **Post-filtering merge**: 15Å centroid distance (`MERGE_DISTANCE_THRESHOLD_A`)

**April 2026 update — paralog specificity:**
The specificity filter was reframed from a "15% max non-specific, size-scaled" rule to a flat **85% max shared with off-targets** (or equivalently, min 15% unique). Size-scaling was removed for specificity because it inflated the effective threshold for large patches beyond what the literature supports. The 85% default is supported by empirical antibody data: an antibody footprint of ~15-25 residues with just 2-4 unique hot-spot contacts (Clackson & Wells, 1995) can be selective. The webapp slider now reads **"Max % shared with off-targets"** (50-95%, default 85%). Higher = more permissive.

A 20-residue patch with 3 cyno mismatches = 15.0% → **PASS** (at threshold)
A 25-residue patch with 4 cyno mismatches = 16.0% → **REJECT** (above threshold)

### Tuning Guidance

**If too restrictive** (rejecting viable targets):
- Increase thresholds to 20-25%
- Biological justification: VHH binding is forgiving; not all surface residues are contacted, and not all contacted residues are critical

**If too permissive** (accepting risky patches):
- Decrease thresholds to 10%
- Use per-residue diagnostic visualizations to inspect exact problem locations

**Parameter sweep**: Run the same target at multiple thresholds (10%, 15%, 20%, 25%) and compare patch counts to find the right balance for your portfolio.

### Diagnostic Visualizations

All per-residue conservation and specificity data is preserved for diagnostic purposes. The epitope map tracks show raw per-residue data (mint = pass, red = fail), allowing you to see exactly where problematic residues are located, even when a patch passes the threshold. The threshold scales with patch size so that larger patches (where the VHH can find a clean ~20-residue footprint) tolerate more total mismatches.

## Output Files

Each run creates a date-stamped directory under `runs/`:

```
runs/YYMMDD_HHMM_erbb2_egfr/
├── .epitope_candidates.csv         # Main results table (hidden)
├── Figures/
│   ├── erbb2_epitope_map.png       # Linear epitope map
│   ├── egfr_epitope_map.png
│   ├── scoring_summary.png         # Multi-target druggability comparison
│   └── BLAST/
│       ├── erbb2_blast_offtargets.png
│       └── egfr_blast_offtargets.png
├── Structures/
│   ├── erbb2_epitope.pml           # PyMOL session script (visible)
│   ├── .erbb2_epitope.pdb          # Annotated PDB (hidden, loaded by PML)
│   └── ...
└── Supplementary Files/
    ├── epitope_candidates.xlsx     # Color-formatted Excel workbook
    ├── Annotated Sequences/
    │   ├── erbb2_epitope.fasta     # FASTA with epitope annotations
    │   ├── erbb2_residue_table.csv # Per-residue: topology, SASA, conservation
    │   └── ...
    ├── Annotations/
    │   ├── erbb2_annotation.json   # Full intermediate results (JSON)
    │   └── ...
    ├── BLAST/
    │   ├── erbb2_blast_hsps.csv    # All BLAST HSPs with identity, range, and call
    │   ├── erbb2_blast_specificity.csv
    │   └── ...
    └── Logs/
        ├── log.txt                 # Full pipeline log
        └── input_manifest.json     # Parameters for reproducibility
```

### BLAST Detail Files

The `Supplementary Files/BLAST/` directory contains two files per target for full transparency into the specificity analysis:

**`{gene}_blast_hsps.csv`** — One row per BLAST HSP (high-scoring segment pair). Each row shows the off-target protein, alignment range, identity percentage, and call:
- `non-specific` — identity >= 70%, covered residues marked non-specific
- `specific` — identity < 70% but passes pre-filter (>= 40% identity, >= 30aa)
- `filtered` — below pre-filter thresholds, not used for scoring

**`{gene}_blast_specificity.csv`** — One row per residue in the target sequence. Shows the max off-target identity at each position, which protein drives that score, and the final binary call (specific/non-specific/not_assessed) against the 70% threshold.

### Annotated PDB Scoring Scheme

The B-factor column in annotated PDBs encodes a tiered epitope analysis:

| B-factor | Meaning |
|----------|---------|
| -20 | Non-target chain (antibody, etc.) |
| 0 | Intracellular |
| 10 | Transmembrane |
| 25 | Extracellular (below distance cutoff) |
| 40 | Distance-qualified (≥ threshold) |
| 55 | + Surface exposed (rel SASA > 25%) |
| 70 | + Cyno conserved |
| 80-95 | Epitope patch (composite_score × 100) |
| 99 | Membrane reference (3 most proximal) |

Each annotated PDB comes with a companion `.pml` script for PyMOL. The full structure is shown as white cartoon with green epitope patch surfaces (30% transparency). Selections are organized into collapsible groups: `Epitopes` (per-patch selections + combined) and `Filters` (pipeline tier selections, clickable but no default coloring). PDB files are dot-prefixed (hidden from file browsers) and loaded automatically by the PML script.

**Cleaved region handling**: Signal peptide residues (N-terminal) and GPI anchor signal residues (C-terminal, after omega site) are hidden in PyMOL and excluded from the distance plot, since these regions are cleaved in the mature protein and never present on the cell surface.

**PDB orientation**: All topology types have their PDB coordinates pre-rotated for consistent PyMOL orientation (membrane horizontal, ectodomain up):
- **Single-pass TM / GPI**: Ecto-axis (anchor → farthest ECD residue) aligned to Y-up, anchor at origin
- **Multi-pass TM**: Membrane normal aligned to Y-up, membrane center at origin

### Scoring Components

| Component | Weight | Description |
|-----------|--------|-------------|
| Area | 60% | Log-scaled patch size (0.0 at 600 A², 1.0 at 5000 A²) |
| Conservation | 25% | Cyno sequence identity within the patch |
| Specificity | 15% | Per-residue max off-target identity across patch |

## Configuration

All thresholds are defined in `config.py` and can be overridden via `run_pipeline()` arguments:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ECTODOMAIN_MIN_DISTANCE_A` | 80.0 | Min Angstroms from membrane |
| `RESOLUTION_THRESHOLD_A` | 3.5 | Max PDB resolution accepted |
| `VHH_FOOTPRINT_MIN_A2` | 600.0 | Min surface patch area (A²) |
| `MAX_CYNO_MISMATCH_PERCENT` | 15.0 | Base % cyno mismatches (scales with patch size: `min(base * sqrt(n/20), 30%)`) |
| `MAX_NONSPECIFIC_PERCENT` | 15.0 | Base % non-specific residues (same scaling formula) |
| `MERGE_DISTANCE_THRESHOLD_A` | 15.0 | Merge patches with centroids within 15Å (post-filtering) |
| `SPECIFICITY_IDENTITY_THRESHOLD` | 0.70 | Off-target BLAST identity for marking residues as non-specific |
| `PATCH_CLUSTERING_DISTANCE_A` | 8.0 | Cα-Cα distance for patch connectivity |
| `MIN_ECTODOMAIN_COVERAGE` | 0.80 | Min PDB ectodomain coverage (applies when `force_experimental=True`) |
| `MAX_CYNO_MISMATCHES_PER_600A2` | 2 | **DEPRECATED** — Legacy sliding window threshold (not used) |
| `MAX_NONSPECIFIC_PER_600A2` | 2 | **DEPRECATED** — Legacy sliding window threshold (not used) |

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

**Path with spaces**: The local BLAST database path is resolved automatically at runtime. If the path contains spaces (e.g., Dropbox directories), the pipeline creates a temporary symlink under `/tmp/` so that `blastp` can parse it correctly. No manual symlink setup is required.

**Note**: SwissProt is updated quarterly by UniProt. Re-run the setup script periodically to refresh the database.

## Epitope Map Figure

Each target gets a 6-track linear sequence map:

| Track | Content | Details |
|-------|---------|---------|
| **Domains** | Topology brackets + domain blocks | Extracellular/Cytoplasmic brackets above; InterPro domains (teal), UniProt domains (blue), TM (purple), SP (gray), Disordered (gray) blocks |
| **Distance** | Distance from membrane (A) | Continuous line through all residues; ECD regions in teal fill, TM/cytoplasmic in light gray fill; dashed threshold line at the configured min distance |
| **SASA** | Per-residue solvent-accessible surface area (A²) | Bar chart (teal) |
| **Cyno Conservation** | Per-residue identity with cynomolgus ortholog | Mint = conserved, red = mismatch. Threshold scales with patch size (`min(base% * sqrt(n/20), 30%)`) |
| **Human Specificity** | Per-residue off-target screening (BLAST) | Mint = specific, red = non-specific (≥70% identity off-target HSP). Same size-scaled threshold |
| **Target Epitope** | Final qualifying patches (post-merge) | Green = passing all filters. Adjacent patches merged if centroids within 15Å. Label: "Patch N — total A² \| score X.XX" |

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

For **single-pass TM** proteins, distances are measured as the **Euclidean (3D) distance from the TM/ectodomain boundary Cα atom** (the anchor residue). This captures the true spatial extent of kinked or bent ectodomains that project sideways rather than straight up — the normal-projected distance underestimates these (e.g., ROR1 projects 54A along the membrane normal but extends 125A in 3D from its TM anchor).

For **multi-pass TM** and **GPI-anchored** proteins, the projected distance along the membrane normal from the bilayer surface is used (no single TM anchor available).

The membrane surface is defined geometrically from the membrane half-thickness (15A default), which is more robust than using the most proximal ectodomain residue as a proxy — AlphaFold models predict in vacuum without membrane context, so juxtamembrane residues can artifactually collapse to membrane-plane height (observed for ERBB2 and EGFR).

The `SpatialFilter` result stores both metrics: `residue_distances` (primary — Euclidean for single-pass TM, projected for others) and `residue_distances_projected` (always projected, for diagnostics). It also records `anchor_resnum` (TM/ecto boundary) and `farthest_resnum` (most distal ECD residue) when available.

## Target Resolution

### Gene Names and UniProt IDs

The pipeline accepts three input formats:

1. **Gene symbols** (e.g., `ERBB2`) — searched against reviewed human UniProt entries
2. **UniProt accessions** (e.g., `P04626`) — used directly
3. **Isoform notation** (e.g., `CLDN18.2` or `P56856-2`) — auto-detected and resolved to isoform-specific entries

**Gene name validation** (March 2026): When searching by gene name, the pipeline now prioritizes results where the requested name is the **primary gene name**, not just an alias. This prevents incorrect mappings like `LY6G6D` → `Q5SQ64` (LY6G6F, where LY6G6D is only an alias) instead of the correct `LY6G6D` → `O95868`. If no exact primary name match is found, the pipeline falls back to the first search result with a warning.

**Isoform support** (March 2026): The pipeline detects isoform notation in two formats:
- **Dot notation**: `CLDN18.2` → searches for gene `CLDN18`, appends isoform `-2` to get `P56856-2`
- **Dash notation**: `P56856-2` → recognized as UniProt isoform ID and used directly

**Known limitation**: UniProt isoform entries (e.g., P56856-2) may lack feature annotations like transmembrane regions that are present in the canonical entry. For membrane proteins specified via isoform notation, TM annotations must be manually verified.

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

# Multiple pairs (all results in one run directory)
python -m epitope_pipeline.bispecific ERBB2:NECTIN4 ERBB2:MSLN

# Custom distance thresholds
python -m epitope_pipeline.bispecific ERBB2:NECTIN4 --distal 80 --proximal 30
```

### How it works

1. **Shared steps**: Target resolution, structure acquisition, and membrane annotation are performed once per unique target
2. **Zone analysis**: Each target is analyzed in both distal (>=60A) and proximal (<=40A) modes — Steps 4-8 run independently per zone
3. **Pair evaluation**: Both orientations evaluated per pair (A-distal/B-proximal and B-distal/A-proximal)
4. **Scoring**: Geometric mean of best distal and proximal composite scores. 20% flexibility bonus when both orientations are valid.

### Outputs

```
runs/YYMMDD_HHMM_bispecific_erbb2_nectin/
├── Figures/
│   ├── bispecific_combined_*.png      # Side-by-side orientation comparison
│   ├── bispecific_*_distal__*.png     # Stacked dual-panel 6-track maps
│   ├── bispecific_pair_summary.png    # Bar chart comparing pair scores
│   └── zone_details/                  # Per-zone deep-dive maps
│       ├── erbb2_distal_epitope_map.png
│       └── ...
├── Structures/
│   ├── *_distal__*_proximal.pml       # Dual PML scripts (shared membrane CGO)
│   └── zone_sessions/                 # Per-zone PML + hidden PDB pairs
│       ├── erbb2_distal_epitope.pml
│       ├── .erbb2_distal_epitope.pdb  # Hidden — loaded by PML
│       └── ...
└── Supplementary Files/
    ├── bispecific_pairs.xlsx           # Color-formatted pair scores
    ├── BLAST/                         # Per unique target
    │   ├── ERBB2_blast_hsps.csv
    │   ├── ERBB2_blast_specificity.csv
    │   └── ...
    └── Logs/
        ├── log.txt
        └── input_manifest.json
```

### Dual PML Alignment

All topology types have PDB coordinates pre-rotated for consistent orientation (membrane horizontal, ectodomain up):
- **Single-pass TM / GPI**: Ecto-axis (anchor → farthest) aligned to Y-up, anchor at origin. Bilayer CGO at Y = -half_thickness.
- **Multi-pass TM**: Membrane normal aligned to Y-up, membrane center at origin. Bilayer CGO at Y = 0.

When both targets in a bispecific pair are rotated, translation is purely along X and the shared bilayer is horizontal.

Epitope patches are shown as green surfaces (`show surface`) with 30% transparency. **Note**: PyMOL `ray_trace_mode 3` cannot render molecular surfaces — use `ray_trace_mode 0` or `1` for ray-traced images.

## Known Limitations

| Limitation | Affected Protein Types | Details |
|------------|-----------------|---------|
| **PCA membrane normal for GPI proteins** | Multi-lobe GPI-anchored proteins | For GPI-anchored proteins with multiple structural domains, PCA principal axis may not align with the max protein dimension, leading to underestimated membrane distances. |
| **Multi-pass TM — tiny ECD** | 4+ TM helix proteins | Proteins with 4+ TM helices typically have short extracellular loops (<40A from membrane), making them incompatible with the default 80A distance threshold. |
| **Monomer-only analysis** | Hetero-oligomeric complexes | Heterodimer interfaces (e.g., integrin αVβ6) are not captured — each chain is analyzed independently. |
| **Low cyno conservation** | Targets with clustered mismatches | Proteins with low cynomolgus conservation may fail the size-scaled threshold. The sqrt scaling helps large patches (a 120-residue patch tolerates up to 30% mismatches) but small patches with >15% mismatches will still fail. Consider increasing the base threshold for challenging targets. |
| **Small ectodomains** | Compact membrane proteins | Ectodomains too small to reach the default 80A threshold. Consider lowering `min_distance_a` for targets with known compact ectodomains (<70A). |

## Version Control

This project is tracked in a **private GitHub repo**: `cartographybio/epitope_pipeline`, authenticated as `amartinko-cartography`.

The git repository lives **inside** `epitope_pipeline/` (not the parent `Claude_Code/` directory).

### Commit Guidelines

Every commit message must include specifics about what changed. Do not use vague messages like "update code" or "fix bugs." Each commit should describe:

1. **What was added, changed, or fixed** — name the files and functions affected
2. **Why** — the motivation or user request that prompted the change
3. **Behavioral impact** — what the pipeline does differently now (new outputs, changed thresholds, fixed bugs, etc.)

Example:
```
Add BLAST detail export for specificity transparency

Export per-run blast/ directory with two reference files:
- {gene}_blast_hsps.csv: all HSPs with identity, range, and call
- {gene}_blast_specificity.csv: per-residue identity scores with top hit

Files modified: specificity.py (SpecificityResult fields), export.py
(export_blast_details), export_bispecific.py (wiring)
```
