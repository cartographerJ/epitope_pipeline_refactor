# Epitope Pipeline

Structural bioinformatics pipeline for identifying druggable VHH epitope space on human membrane protein targets.

## Overview

This pipeline takes a list of membrane protein targets and systematically identifies surface regions suitable for VHH antibody binding. It filters for regions that are:

1. **Structurally resolved** — uses AlphaFold DB full-length predictions by default, or experimental PDB structures via opt-in
2. **Distal from the membrane** — ≥80 Angstroms from the membrane surface (default; configurable)
3. **Surface accessible** — solvent-exposed patches large enough for a VHH CDR footprint (≥600 A²)
4. **Cynomolgus-conserved** — ≤15% cyno mismatches per patch (whole-patch evaluation; configurable threshold)
5. **Target-specific** — ≤15% non-specific residues per patch (whole-patch evaluation; ≥70% off-target BLAST identity threshold)

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

## Pipeline Steps

| Step | Module | Description |
|------|--------|-------------|
| 1 | `target_input.py` | Resolve UniProt IDs or gene names to full protein metadata + InterPro domains |
| 2 | `structure.py` | Acquire structure: AlphaFold DB (default), experimental PDB (opt-in), or Tamarind AlphaFold (fallback) |
| 3 | `membrane.py` | Define membrane plane (OPM, TM annotations, or GPI-anchor) |
| 4 | `spatial.py` | Filter to ectodomain residues ≥80A from membrane surface (default; bispecific uses 60A/40A) |
| 5 | `surface.py` | Calculate SASA, cluster into surface patches ≥600 A² |
| 6 | `conservation.py` | Align with cyno ortholog, whole-patch evaluation (≤15% mismatches per patch) |
| 7 | `specificity.py` | Full-sequence BLAST, whole-patch evaluation (≤15% non-specific per patch) + merge adjacent patches |
| 8 | `scoring.py` | Composite score: area + distance + conservation + specificity + accessibility |
| 9 | `export.py` | Generate all output files (CSV, XLSX, PDB, FASTA, JSON) |
| 10 | `visualize.py` | 6-track epitope map with domains, distance, SASA, conservation, specificity |

## Filtering Strategy: Whole-Patch Evaluation

**March 2026 Update**: The pipeline uses a fundamentally different approach for conservation and specificity filtering compared to earlier versions.

### Philosophy: VHH Binding Units

VHH antibodies bind contiguous ~600 A² surface patches as atomic units. The CDR loops contact ~15-25 residues simultaneously, and not all contacted residues are equally critical for binding. This biological reality informs our filtering strategy:

**Whole-patch evaluation** (current approach):
- Each ≥600 A² patch is treated as an indivisible VHH binding surface
- Accept/reject based on **percentage of problematic residues**
- Simple pass/fail: "This patch has 8% cyno mismatches → PASS (under 15% threshold)"
- Patches remain intact throughout filtering — no trimming, no size changes
- Adjacent patches that pass filtering are merged into larger contiguous epitope regions (15Å centroid distance threshold)

**Sliding window trim-and-recluster** (deprecated, pre-March 2026):
- Scanned each residue's ~600 A² neighborhood for mismatches
- If any window exceeded threshold, trimmed failing residues and re-clustered survivors
- Could fragment a 149-residue patch into a 19-residue sub-patch
- Biologically questionable: a VHH that binds the 149-residue surface won't bind the 19-residue fragment

### Default Thresholds

- **Cyno conservation**: ≤15% mismatches per patch (`MAX_CYNO_MISMATCH_PERCENT`)
- **Human specificity**: ≤15% non-specific residues per patch (`MAX_NONSPECIFIC_PERCENT`)
- **Post-filtering merge**: 15Å centroid distance (`MERGE_DISTANCE_THRESHOLD_A`)

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

All per-residue conservation and specificity data is preserved for diagnostic purposes, even though filtering uses whole-patch evaluation. The epitope map tracks show raw per-residue data (mint = pass, red = fail), allowing you to see exactly where problematic residues are located, even when a patch passes the percentage threshold. For example, a patch with 12% mismatches (passing 15% threshold) will still show those specific mismatched residues in red on the conservation track.

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
| `MAX_CYNO_MISMATCH_PERCENT` | 15.0 | Max % cyno mismatches to accept patch (whole-patch evaluation) |
| `MAX_NONSPECIFIC_PERCENT` | 15.0 | Max % non-specific residues to accept patch (whole-patch evaluation) |
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

**Note**: SwissProt is updated quarterly by UniProt. Re-run the setup script periodically to refresh the database.

## Epitope Map Figure

Each target gets a 6-track linear sequence map:

| Track | Content | Details |
|-------|---------|---------|
| **Domains** | Topology brackets + domain blocks | Extracellular/Cytoplasmic brackets above; InterPro domains (teal), UniProt domains (blue), TM (purple), SP (gray), Disordered (gray) blocks |
| **Distance** | Distance from membrane (A) | Continuous line through all residues; ECD regions in teal fill, TM/cytoplasmic in light gray fill; dashed threshold line at the configured min distance |
| **SASA** | Per-residue solvent-accessible surface area (A²) | Bar chart (teal) |
| **Cyno Conservation** | Per-residue identity with cynomolgus ortholog | Mint = conserved, red = mismatch. Shows raw per-residue data for diagnostic purposes; filtering uses whole-patch evaluation (≤15% mismatches per patch) |
| **Human Specificity** | Per-residue off-target screening (BLAST) | Mint = specific, red = non-specific (≥70% identity off-target HSP). Shows raw per-residue data; filtering uses whole-patch evaluation (≤15% non-specific per patch) |
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

## Known Limitations

| Limitation | Affected Protein Types | Details |
|------------|-----------------|---------|
| **PCA membrane normal for GPI proteins** | Multi-lobe GPI-anchored proteins | For GPI-anchored proteins with multiple structural domains, PCA principal axis may not align with the max protein dimension, leading to underestimated membrane distances. |
| **Multi-pass TM — tiny ECD** | 4+ TM helix proteins | Proteins with 4+ TM helices typically have short extracellular loops (<40A from membrane), making them incompatible with the default 80A distance threshold. |
| **Monomer-only analysis** | Hetero-oligomeric complexes | Heterodimer interfaces (e.g., integrin αVβ6) are not captured — each chain is analyzed independently. |
| **Low cyno conservation** | Targets with clustered mismatches | Proteins with regions of low cynomolgus conservation (e.g., divergent gene families) may fail the 15% whole-patch threshold. Consider relaxing to 20-25% for challenging targets. |
| **Small ectodomains** | Compact membrane proteins | Ectodomains too small to reach the default 80A threshold. Consider lowering `min_distance_a` for targets with known compact ectodomains (<70A). |
