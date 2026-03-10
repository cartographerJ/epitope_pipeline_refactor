#!/bin/bash
set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Local BLAST Database Setup for Epitope Pipeline${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Get script directory and project root
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
DB_DIR="$PROJECT_ROOT/blast_db/swissprot"

# 1. Check for BLAST+ installation
echo -e "${YELLOW}[1/6]${NC} Checking for BLAST+ installation..."
if ! command -v blastp &> /dev/null || ! command -v makeblastdb &> /dev/null; then
    echo -e "${RED}ERROR: BLAST+ tools not found.${NC}"
    echo ""
    echo "Please install NCBI BLAST+ tools:"
    echo ""
    echo "  macOS:"
    echo "    brew install blast"
    echo ""
    echo "  Ubuntu/Debian:"
    echo "    sudo apt-get install ncbi-blast+"
    echo ""
    exit 1
fi
echo -e "${GREEN}✓ BLAST+ tools found${NC}"
blastp -version | head -n 1
echo ""

# 2. Create database directory
echo -e "${YELLOW}[2/6]${NC} Creating database directory..."
mkdir -p "$DB_DIR"
cd "$DB_DIR"
echo -e "${GREEN}✓ Directory created: $DB_DIR${NC}"
echo ""

# 3. Download SwissProt FASTA
echo -e "${YELLOW}[3/6]${NC} Downloading SwissProt database..."
SWISSPROT_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

if [ -f "uniprot_sprot.fasta" ]; then
    echo -e "${YELLOW}SwissProt FASTA already exists. Skipping download.${NC}"
    echo -e "${YELLOW}To re-download, delete: $DB_DIR/uniprot_sprot.fasta${NC}"
else
    echo "Downloading from UniProt (this may take 5-10 minutes)..."
    if command -v wget &> /dev/null; then
        wget -O uniprot_sprot.fasta.gz "$SWISSPROT_URL"
    elif command -v curl &> /dev/null; then
        curl -o uniprot_sprot.fasta.gz "$SWISSPROT_URL"
    else
        echo -e "${RED}ERROR: Neither wget nor curl found. Please install one of them.${NC}"
        exit 1
    fi

    echo "Decompressing..."
    gunzip uniprot_sprot.fasta.gz
    echo -e "${GREEN}✓ SwissProt downloaded${NC}"
fi
echo ""

# 4. Filter to human proteins using Python
echo -e "${YELLOW}[4/6]${NC} Filtering to human proteins..."
python3 - <<'PYTHON_SCRIPT'
import sys
from pathlib import Path

# Try to import BioPython
try:
    from Bio import SeqIO
except ImportError:
    print("\033[0;31mERROR: BioPython not found.\033[0m", file=sys.stderr)
    print("\nPlease install BioPython:", file=sys.stderr)
    print("  pip install biopython", file=sys.stderr)
    sys.exit(1)

input_file = "uniprot_sprot.fasta"
output_file = "swissprot_human.fasta"

print(f"Parsing {input_file}...")
human_count = 0
total_count = 0

with open(output_file, 'w') as out_f:
    for record in SeqIO.parse(input_file, "fasta"):
        total_count += 1

        # Check if organism is Homo sapiens
        # UniProt FASTA headers contain OS=Organism Name
        if "OS=Homo sapiens" in record.description or "OX=9606" in record.description:
            SeqIO.write(record, out_f, "fasta")
            human_count += 1

        # Progress indicator
        if total_count % 10000 == 0:
            print(f"  Processed {total_count:,} sequences, found {human_count:,} human proteins...")

print(f"\n\033[0;32m✓ Filtered {human_count:,} human proteins from {total_count:,} total sequences\033[0m")
print(f"Output: {output_file}")
PYTHON_SCRIPT

if [ $? -ne 0 ]; then
    echo -e "${RED}ERROR: Python filtering failed.${NC}"
    exit 1
fi
echo ""

# 5. Build BLAST database
echo -e "${YELLOW}[5/6]${NC} Building BLAST database..."
makeblastdb \
    -in swissprot_human.fasta \
    -dbtype prot \
    -out swissprot_human \
    -title "SwissProt Human Proteins" \
    -parse_seqids

if [ $? -ne 0 ]; then
    echo -e "${RED}ERROR: makeblastdb failed.${NC}"
    exit 1
fi
echo -e "${GREEN}✓ BLAST database built${NC}"
echo ""

# 6. Verify database
echo -e "${YELLOW}[6/6]${NC} Verifying database..."
TEST_QUERY=">test_sequence
MQRGQALWDFPGKLTDSTRKKTGKTERLQVE"

echo "$TEST_QUERY" | blastp -db swissprot_human -outfmt 6 -max_target_seqs 1 > /dev/null 2>&1

if [ $? -ne 0 ]; then
    echo -e "${RED}ERROR: Database verification failed.${NC}"
    exit 1
fi
echo -e "${GREEN}✓ Database verified and functional${NC}"
echo ""

# Summary
echo -e "${BLUE}========================================${NC}"
echo -e "${GREEN}Setup Complete!${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo "Database location: $DB_DIR/swissprot_human"
echo "Database files:"
ls -lh swissprot_human.* | awk '{print "  " $9 " (" $5 ")"}'
echo ""
echo "To use local BLAST in the pipeline, ensure config.py has:"
echo "  USE_LOCAL_BLAST = True"
echo "  LOCAL_BLAST_DB_PATH = \"$DB_DIR/swissprot_human\""
echo ""
echo -e "${GREEN}You can now run the epitope pipeline with local BLAST!${NC}"
echo ""
