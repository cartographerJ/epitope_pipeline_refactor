#!/usr/bin/env python3
"""PCDHB2 per-domain epitope analysis.

Extracts each cadherin domain, runs SASA + BLAST + patch clustering
using the pipeline's internal functions, and generates per-domain results.
"""

import sys, os, tempfile, subprocess, csv, json
import numpy as np
from pathlib import Path

# Add parent to path for pipeline imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from epitope_pipeline.compute.surface import _calculate_sasa, cluster_surface_patches

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BLAST_DB = "/tmp/blast_pcdhb2_db/swissprot_human"
RUN_DIR = Path(__file__).resolve().parent.parent / "runs" / "260407_1010_pcdhb2"
STRUCTURE = list((RUN_DIR / "Structures").glob("*.pdb"))[0]
RES_TABLE = RUN_DIR / "Supplementary Files" / "Annotated Sequences" / "pcdhb2_residue_table.csv"
OUT_DIR = RUN_DIR / "domain_analysis"
OUT_DIR.mkdir(exist_ok=True)

# Max residue SASA for relative SASA calculation (Tien et al. 2013)
MAX_ASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167,
    "E": 223, "Q": 225, "G": 104, "H": 224, "I": 197,
    "L": 201, "K": 236, "M": 224, "F": 240, "P": 159,
    "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}

# PCDHB2 cadherin domain boundaries (UniProt Q9Y5E7)
DOMAINS = [
    ("EC1",  1,   108),
    ("EC2",  109, 212),
    ("EC3",  213, 325),
    ("EC4",  326, 432),
    ("EC5",  433, 543),
    ("EC6",  544, 660),
]

MIN_PATCH_AREA = 600.0
SASA_REL_THRESHOLD = 0.25


def read_sequence():
    seq = {}
    with open(RES_TABLE) as f:
        reader = csv.DictReader(f)
        for row in reader:
            seq[int(row["residue_num"])] = row["aa"]
    return seq


def extract_domain_pdb(pdb_path, start, end, out_path):
    """Extract residues start-end from PDB chain A."""
    with open(pdb_path) as f_in, open(out_path, "w") as f_out:
        for line in f_in:
            if line.startswith(("ATOM", "HETATM")):
                resnum = int(line[22:26].strip())
                if start <= resnum <= end:
                    f_out.write(line)
        f_out.write("END\n")


def get_ca_coords(pdb_path, start, end):
    """Extract Calpha coordinates from PDB."""
    coords = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                resnum = int(line[22:26].strip())
                if start <= resnum <= end:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords[resnum] = np.array([x, y, z])
    return coords


def compute_patch_area(pdb_path, patch_residues):
    """Compute SASA of a patch (sum of per-residue SASA for patch members)."""
    sasa = _calculate_sasa(pdb_path, "A")
    return sum(sasa.get(r, 0) for r in patch_residues)


def blast_domain(name, sequence):
    """BLAST a domain sequence, return list of (hit_name, pident)."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        tmp.write(">{}\n{}\n".format(name, sequence))
        tmp_path = tmp.name

    cmd = [
        "blastp", "-db", BLAST_DB, "-query", tmp_path,
        "-outfmt", "6 sseqid stitle pident length evalue",
        "-evalue", "1e-5", "-max_target_seqs", "30", "-num_threads", "2",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    os.unlink(tmp_path)

    hits = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        desc = parts[1].split(" OS=")[0]
        pident = float(parts[2])
        if "Protocadherin beta-2" in desc:
            continue
        hits.append((desc, pident))
    return hits


def run_domain(name, start, end, full_seq):
    """Full analysis for one domain."""
    print("\n" + "=" * 60)
    print("  {} (residues {}-{})".format(name, start, end))
    print("=" * 60)

    # Extract domain PDB
    domain_pdb = OUT_DIR / "{}.pdb".format(name.lower())
    extract_domain_pdb(str(STRUCTURE), start, end, str(domain_pdb))

    # Domain sequence
    subseq = "".join(full_seq[r] for r in range(start, end + 1) if r in full_seq)
    print("  Sequence length: {} aa".format(len(subseq)))

    # SASA on domain PDB
    sasa = _calculate_sasa(str(domain_pdb), "A")
    exposed = []
    for r in range(start, end + 1):
        aa = full_seq.get(r, "A")
        max_asa = MAX_ASA.get(aa, 200)
        abs_sasa = sasa.get(r, 0)
        rel_sasa = abs_sasa / max_asa if max_asa > 0 else 0
        if rel_sasa > SASA_REL_THRESHOLD:
            exposed.append(r)
    print("  Surface-exposed residues (>25% rSASA): {}/{}".format(len(exposed), end - start + 1))

    # Calpha coords for clustering
    ca_coords = get_ca_coords(str(STRUCTURE), start, end)

    # Patch clustering
    if len(exposed) >= 2:
        clusters = cluster_surface_patches(exposed, ca_coords)
    else:
        clusters = [exposed] if exposed else []

    # Filter by area and compute patch areas
    patches = []
    for cluster in clusters:
        area = compute_patch_area(str(domain_pdb), cluster)
        if area >= MIN_PATCH_AREA:
            patches.append({"residues": sorted(cluster), "area": area})

    print("  Patches >= {} A2: {}".format(int(MIN_PATCH_AREA), len(patches)))
    for i, p in enumerate(patches):
        rng = p["residues"]
        print("    P{}: {} residues, {:.0f} A2, res {}-{}".format(
            i, len(rng), p["area"], min(rng), max(rng)))

    # BLAST full domain
    hits = blast_domain(name, subseq)
    hits_above_70 = [(h, p) for h, p in hits if p >= 70]
    top = hits[0] if hits else ("---", 0)
    print("  BLAST (full domain): {} hits >= 70% identity".format(len(hits_above_70)))
    print("    Top hit: {} ({:.1f}%)".format(top[0], top[1]))

    # Per-patch specificity BLAST
    for i, p in enumerate(patches):
        patch_res = p["residues"]
        patch_seq = "".join(full_seq.get(r, "X") for r in patch_res)
        patch_hits = blast_domain("{}_P{}".format(name, i), patch_seq)
        patch_top = patch_hits[0] if patch_hits else ("---", 0)
        patch_above_70 = len([h for h, pid in patch_hits if pid >= 70])
        print("    P{} BLAST: top {:.1f}% ({}), {} hits >=70%".format(
            i, patch_top[1], patch_top[0][:35], patch_above_70))

    return {
        "domain": name,
        "start": start,
        "end": end,
        "n_residues": end - start + 1,
        "n_exposed": len(exposed),
        "n_patches": len(patches),
        "patches": [{"residues": p["residues"], "area": p["area"]} for p in patches],
        "top_blast_hit": top[0],
        "top_blast_pident": top[1],
        "n_hits_above_70": len(hits_above_70),
    }


def main():
    # Ensure symlink for BLAST DB
    if not os.path.exists("/tmp/blast_pcdhb2_db"):
        os.symlink(
            str(Path(__file__).resolve().parent.parent / "blast_db" / "swissprot"),
            "/tmp/blast_pcdhb2_db",
        )

    full_seq = read_sequence()
    print("PCDHB2 Per-Domain Epitope Analysis")
    print("Structure: {}".format(STRUCTURE.name))
    print("Total residues: {}".format(len(full_seq)))

    results = []
    for name, start, end in DOMAINS:
        r = run_domain(name, start, end, full_seq)
        results.append(r)

    # Summary
    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print("{:<8} {:<10} {:<10} {:<10} {:<12} {:<10} {}".format(
        "Domain", "Residues", "Exposed", "Patches", "Top Hit %", "#>=70%", "Verdict"))
    print("-" * 80)
    for r in results:
        if r["n_hits_above_70"] == 0 and r["n_patches"] > 0:
            verdict = "TARGETABLE"
        elif r["n_hits_above_70"] == 0:
            verdict = "SPECIFIC (no patches)"
        elif r["n_hits_above_70"] <= 2:
            verdict = "MARGINAL"
        else:
            verdict = "NON-SPECIFIC"
        print("{:<8} {}-{:<6} {:<10} {:<10} {:<12.1f} {:<10} {}".format(
            r["domain"], r["start"], r["end"], r["n_exposed"],
            r["n_patches"], r["top_blast_pident"], r["n_hits_above_70"], verdict))

    # Save JSON
    summary_path = OUT_DIR / "domain_summary.json"
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print("\nResults saved to: {}".format(OUT_DIR))


if __name__ == "__main__":
    main()
