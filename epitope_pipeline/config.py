"""
Epitope pipeline configuration.

All tunable parameters, API URLs, path conventions, and thresholds.
"""

from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# Repo root (one level above the package directory). Used as the default base
# for run output, cache, and the local BLAST DB. Works for editable installs
# from a checkout; if this package is ever site-installed in a venv where
# writing into the install location isn't desirable, override these via env
# vars or a user config.
PIPELINE_ROOT = Path(__file__).resolve().parent.parent

# Output directories
RUNS_DIR = PIPELINE_ROOT / "runs"
CACHE_DIR = PIPELINE_ROOT / "cache"

# ---------------------------------------------------------------------------
# External API URLs
# ---------------------------------------------------------------------------
UNIPROT_API = "https://rest.uniprot.org"
OPM_API = "https://opm.phar.umich.edu"
RCSB_API = "https://data.rcsb.org/rest/v1"
RCSB_SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_FILES = "https://files.rcsb.org/download"
NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
ALPHAFOLD_DB_URL = "https://alphafold.ebi.ac.uk/files"

# ---------------------------------------------------------------------------
# Taxonomy
# ---------------------------------------------------------------------------
HUMAN_TAXID = 9606
CYNO_TAXID = 9541          # Macaca fascicularis (cynomolgus macaque)
RHESUS_TAXID = 9544        # Macaca mulatta (rhesus macaque) — fallback for cyno

# ---------------------------------------------------------------------------
# Structure acquisition
# ---------------------------------------------------------------------------
RESOLUTION_THRESHOLD_A = 3.5        # Max resolution (Angstroms) to accept
ALLOWED_METHODS = [                 # No NMR
    "X-RAY DIFFRACTION",
    "ELECTRON MICROSCOPY",
]
MIN_COVERAGE_FRACTION = 0.30        # Prefer AlphaFold if best PDB covers <30% of protein
MIN_ECTODOMAIN_COVERAGE = 0.80      # Require >=80% of ectodomain covered by PDB

# ---------------------------------------------------------------------------
# Membrane / spatial thresholds
# ---------------------------------------------------------------------------
ECTODOMAIN_MIN_DISTANCE_A = 80.0    # Min Angstroms from membrane plane
MEMBRANE_SLAB_HALF_THICKNESS = 15.0 # Default half-thickness of membrane (A)

# ---------------------------------------------------------------------------
# Surface analysis
# ---------------------------------------------------------------------------
SASA_PROBE_RADIUS = 1.4             # Standard water probe radius (A)
SURFACE_EXPOSURE_THRESHOLD = 0.25   # Fraction of max SASA for "exposed"
VHH_FOOTPRINT_MIN_A2 = 600.0       # Minimum patch area (A^2) for VHH CDR
VHH_FOOTPRINT_MAX_A2 = 900.0       # Typical max VHH footprint (A^2)
PATCH_CLUSTERING_DISTANCE_A = 8.0   # Max CA-CA distance for same surface patch

# Reference max SASA values per amino acid (Gly-X-Gly tripeptide, Lee & Richards)
MAX_SASA = {
    "ALA": 113.0, "ARG": 241.0, "ASN": 158.0, "ASP": 151.0, "CYS": 140.0,
    "GLU": 183.0, "GLN": 189.0, "GLY":  85.0, "HIS": 194.0, "ILE": 182.0,
    "LEU": 180.0, "LYS": 211.0, "MET": 204.0, "PHE": 218.0, "PRO": 143.0,
    "SER": 122.0, "THR": 146.0, "TRP": 259.0, "TYR": 229.0, "VAL": 160.0,
}
# Single-letter to three-letter mapping for SASA lookup
AA_1TO3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}

# ---------------------------------------------------------------------------
# Conservation
# ---------------------------------------------------------------------------
# Whole-patch evaluation (NEW)
MAX_CYNO_MISMATCH_PERCENT = 15.0        # Max % cyno mismatches to accept patch

# DEPRECATED: Legacy sliding window threshold (not used in whole-patch mode)
MAX_CYNO_MISMATCHES_PER_600A2 = 2       # Kept for backward compatibility only

# ---------------------------------------------------------------------------
# Specificity (BLAST)
# ---------------------------------------------------------------------------
# Per-paralog patch rule: a patch fails if any single paralog matches
# more than MAX_NONSPECIFIC_PERCENT of the patch residues. HSPs below
# 40% identity or 30 aa alignment are filtered out before counting.
MAX_NONSPECIFIC_PERCENT = 85.0          # Max % of patch that can match any single paralog

# Post-filtering merge
MERGE_DISTANCE_THRESHOLD_A = 15.0       # Merge patches with centroids within 15Å
BLAST_EVALUE_CUTOFF = 0.001
BLAST_WORD_SIZE = 3
BLAST_DATABASE = "swissprot"
BLAST_ORGANISM_FILTER = ""              # Post-filter to human (entrez_query broken on NCBI)
BLAST_HITLIST_SIZE = 50
BLAST_DELAY_SECONDS = 3                 # NCBI rate limit

# Local BLAST database configuration
USE_LOCAL_BLAST = True                  # Toggle: True = local blastp, False = remote NCBI API
# Spaces in path are handled automatically (temp symlink created at runtime)
LOCAL_BLAST_DB_PATH = PIPELINE_ROOT / "blast_db" / "swissprot" / "swissprot_human"

# ---------------------------------------------------------------------------
# Scoring weights — DEPRECATED: not used by scoring.py (hardcoded 60/25/15)
# Kept for backward compatibility with external callers only.
# ---------------------------------------------------------------------------
SCORE_WEIGHTS = {
    "area":          0.60,
    "conservation":  0.25,
    "specificity":   0.15,
}

# ---------------------------------------------------------------------------
# Bispecific / dual-targeting thresholds
# ---------------------------------------------------------------------------
DISTAL_MIN_DISTANCE_A = 60.0         # Distal mode: ECD residues >= 60A
PROXIMAL_MAX_DISTANCE_A = 50.0       # Proximal mode: ECD residues <= 50A
FLEXIBILITY_BONUS = 1.20             # 20% bonus when BOTH orientations are valid
DUAL_PML_GAP_A = 120.0              # X-axis translation gap between structures in dual PML

# ---------------------------------------------------------------------------
# Visualization palette (Cartography)
# ---------------------------------------------------------------------------
PALETTE = {
    "light_gray":  "#D3D3D3",
    "mint":        "#E0F3DB",
    "green":       "#6BC291",
    "teal":        "#18B5CB",
    "blue":        "#2E95D2",
    "dark_purple": "#28154C",
}

# Semantic color assignments
COLOR_EXTRACELLULAR = PALETTE["teal"]
COLOR_TRANSMEMBRANE = PALETTE["dark_purple"]
COLOR_INTRACELLULAR = PALETTE["light_gray"]
COLOR_EPITOPE_PATCH = PALETTE["green"]
COLOR_CONSERVED = PALETTE["mint"]
COLOR_MISMATCH = "#E74C3C"  # Red for conservation mismatches

# ---------------------------------------------------------------------------
# API retry settings
# ---------------------------------------------------------------------------
RETRY_UNIPROT = {"attempts": 3, "backoff_seconds": 5}
RETRY_RCSB = {"attempts": 3, "backoff_seconds": 2}
RETRY_BLAST = {"attempts": 2, "backoff_seconds": 30}
