"""
Epitope pipeline configuration.

All tunable parameters, API URLs, path conventions, and thresholds.
Tamarind API key is loaded from the existing tamarind/.env file.
"""

import os
from pathlib import Path
from dotenv import load_dotenv

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PIPELINE_ROOT = Path(__file__).parent
PROJECT_ROOT = PIPELINE_ROOT.parent
TAMARIND_ROOT = PROJECT_ROOT / "tamarind"

# Load Tamarind API key from existing .env
load_dotenv(TAMARIND_ROOT / ".env")
TAMARIND_API_KEY = os.environ.get("TAMARIND_API_KEY", "")

# Output directories
RUNS_DIR = PIPELINE_ROOT / "runs"
RAW_JOBS_DIR = PIPELINE_ROOT / "raw_jobs"
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
FOLDING_TOOL = "alphafold"          # Tamarind folding service name
FOLDING_NUM_MODELS = 1
FOLDING_NUM_RECYCLES = 3
FOLDING_TIMEOUT = 7200              # 2 hours max wait

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
SPECIFICITY_IDENTITY_THRESHOLD = 0.7    # >70% identity -> non-specific residue

# Whole-patch evaluation (NEW)
MAX_NONSPECIFIC_PERCENT = 15.0          # Max % non-specific residues to accept patch

# DEPRECATED: Legacy sliding window threshold (not used in whole-patch mode)
MAX_NONSPECIFIC_PER_600A2 = 2           # Kept for backward compatibility only

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
LOCAL_BLAST_DB_PATH = PIPELINE_ROOT / "blast_db" / "swissprot" / "swissprot_human"

# ---------------------------------------------------------------------------
# Scoring weights (composite epitope score)
# ---------------------------------------------------------------------------
SCORE_WEIGHTS = {
    "area":          0.25,
    "distance":      0.20,
    "conservation":  0.25,
    "specificity":   0.20,
    "accessibility": 0.10,
}

# ---------------------------------------------------------------------------
# Bispecific / dual-targeting thresholds
# ---------------------------------------------------------------------------
DISTAL_MIN_DISTANCE_A = 60.0         # Distal mode: ECD residues >= 60A
PROXIMAL_MAX_DISTANCE_A = 40.0       # Proximal mode: ECD residues <= 40A
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
