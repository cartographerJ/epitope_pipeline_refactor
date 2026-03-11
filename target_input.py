"""
Step 1: Target Resolution — resolve protein identifiers to full UniProt entries.

Accepts UniProt accession IDs (e.g. P04626) or gene names (e.g. ERBB2),
auto-detects format, and fetches full metadata including sequence, topology
features, PDB cross-references, and cynomolgus monkey ortholog sequence.

All API responses are cached to cache/uniprot/ for efficient re-runs.
"""

import json
import logging
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

import requests

from epitope_pipeline.config import (
    CACHE_DIR,
    CYNO_TAXID,
    HUMAN_TAXID,
    RETRY_UNIPROT,
    UNIPROT_API,
)

logger = logging.getLogger(__name__)

# UniProt accession regex patterns (6 or 10 char, with optional isoform suffix)
_ACCESSION_RE = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9](-\d+)?$"           # 6-char (e.g. P04626 or P56856-2)
    r"|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}(-\d+)?$"  # 10-char with optional isoform
)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class TargetInfo:
    """Resolved target protein metadata from UniProt."""
    uniprot_id: str                         # e.g. "P04626"
    gene_name: str                          # e.g. "ERBB2"
    protein_name: str                       # e.g. "Receptor tyrosine-protein kinase erbB-2"
    organism: str                           # e.g. "Homo sapiens"
    sequence: str                           # Full amino acid sequence
    sequence_length: int = 0
    features: list = field(default_factory=list)  # Raw UniProt features (TM, topological, etc.)
    pdb_ids: list = field(default_factory=list)   # Known PDB cross-references
    cyno_uniprot_id: Optional[str] = None   # Cyno ortholog accession
    cyno_sequence: Optional[str] = None     # Cyno ortholog sequence
    ectodomain_ranges: list = field(default_factory=list)  # [(start, end), ...] extracellular regions
    ectodomain_length: int = 0              # Total extracellular residue count

    def __post_init__(self):
        if not self.sequence_length:
            self.sequence_length = len(self.sequence)


class TargetResolutionError(Exception):
    """Could not resolve a target identifier to a UniProt entry."""
    pass


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def resolve_targets(identifiers):
    """
    Resolve a list of protein identifiers to TargetInfo objects.

    Each identifier is auto-detected as either a UniProt accession ID or a
    gene name. Gene names are searched against reviewed human Swiss-Prot entries.

    Args:
        identifiers: List of UniProt IDs or gene names (can be mixed).

    Returns:
        List of TargetInfo objects, one per successfully resolved target.

    Raises:
        TargetResolutionError: If any identifier cannot be resolved.
    """
    targets = []
    for ident in identifiers:
        ident = ident.strip()
        logger.info("Resolving target: %s", ident)

        if _is_accession(ident):
            uniprot_id = ident.upper()
        else:
            # Check for isoform notation (e.g., CLDN18.2 or CLDN18-2)
            isoform_suffix = ""
            if "." in ident and ident.split(".")[-1].isdigit():
                # GENE.X format (e.g., CLDN18.2)
                gene_base, isoform_num = ident.rsplit(".", 1)
                isoform_suffix = f"-{isoform_num}"
                ident = gene_base
                logger.info("  Detected isoform notation: %s → %s (isoform %s)",
                           f"{gene_base}.{isoform_num}", gene_base, isoform_num)

            base_id = _search_uniprot_by_gene(ident)
            uniprot_id = base_id + isoform_suffix

        # Fetch full UniProt entry
        entry = _fetch_uniprot_entry(uniprot_id)
        target = _parse_uniprot_entry(entry)

        # Resolve cynomolgus monkey ortholog
        cyno_id, cyno_seq = _resolve_cyno_ortholog(target.gene_name)
        target.cyno_uniprot_id = cyno_id
        target.cyno_sequence = cyno_seq
        if cyno_id:
            logger.info("  Cyno ortholog: %s (%d aa)", cyno_id, len(cyno_seq))
        else:
            logger.warning("  No cyno ortholog found for %s", target.gene_name)

        logger.info(
            "  Resolved: %s (%s) — %d aa, %d PDB xrefs, ectodomain %d aa (%s)",
            target.gene_name, target.uniprot_id,
            target.sequence_length, len(target.pdb_ids),
            target.ectodomain_length,
            ", ".join("{}-{}".format(s, e) for s, e in target.ectodomain_ranges),
        )
        targets.append(target)

    return targets


# ---------------------------------------------------------------------------
# Identifier detection
# ---------------------------------------------------------------------------

def _is_accession(identifier):
    """Check if the identifier matches UniProt accession format."""
    return bool(_ACCESSION_RE.match(identifier.strip().upper()))


# ---------------------------------------------------------------------------
# UniProt API calls (with caching and retry)
# ---------------------------------------------------------------------------

def _fetch_uniprot_entry(uniprot_id):
    """
    Fetch a full UniProt entry as JSON.

    Caches responses to cache/uniprot/{accession}.json so re-runs skip
    the API call.

    Args:
        uniprot_id: UniProt accession (e.g. "P04626").

    Returns:
        Parsed JSON dict of the UniProt entry.

    Raises:
        TargetResolutionError: If the entry cannot be fetched.
    """
    cache_path = CACHE_DIR / "uniprot" / "{}.json".format(uniprot_id)
    if cache_path.exists():
        logger.debug("  Cache hit: %s", cache_path)
        with open(cache_path) as f:
            return json.load(f)

    url = "{}/uniprotkb/{}.json".format(UNIPROT_API, uniprot_id)
    data = _api_get(url, "UniProt entry {}".format(uniprot_id))

    # Cache the response
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(data, f, indent=2)

    return data


def _search_uniprot_by_gene(gene_name):
    """
    Search UniProt for a human reviewed entry matching the gene name.

    Args:
        gene_name: Gene symbol (e.g. "ERBB2").

    Returns:
        UniProt accession string.

    Raises:
        TargetResolutionError: If no matching entry found.
    """
    # Search for exact gene name in reviewed human entries
    query = "gene:{} AND organism_id:{} AND reviewed:true".format(
        gene_name, HUMAN_TAXID
    )
    url = "{}/uniprotkb/search".format(UNIPROT_API)
    params = {
        "query": query,
        "format": "json",
        "size": 5,
        "fields": "accession,gene_names,protein_name,organism_name",
    }

    data = _api_get(url, "UniProt search for gene '{}'".format(gene_name), params=params)
    results = data.get("results", [])

    if not results:
        raise TargetResolutionError(
            "No reviewed human UniProt entry found for gene name '{}'. "
            "Try using a UniProt accession ID directly.".format(gene_name)
        )

    # Find result where primary gene name matches (not just an alias)
    exact_match = None
    for result in results:
        gene_names = result.get("genes", [])
        if gene_names:
            primary_gene = gene_names[0].get("geneName", {}).get("value", "")
            if primary_gene.upper() == gene_name.upper():
                exact_match = result
                break

    # Use exact match if found, otherwise fall back to first result with warning
    if exact_match:
        result = exact_match
    else:
        result = results[0]
        gene_names = result.get("genes", [])
        if gene_names:
            primary_gene = gene_names[0].get("geneName", {}).get("value", "")
            logger.warning(
                "  No exact primary gene name match for '%s'. Using UniProt %s (primary name: '%s') "
                "where '%s' is listed as an alias. Consider using '%s' or UniProt ID directly.",
                gene_name, result["primaryAccession"], primary_gene, gene_name, primary_gene
            )

    accession = result["primaryAccession"]
    logger.info("  Gene '%s' resolved to UniProt %s", gene_name, accession)
    return accession


def _resolve_cyno_ortholog(gene_name):
    """
    Find the cynomolgus macaque ortholog for a given gene.

    Searches UniProt for the gene in Macaca fascicularis (taxid 9541).
    Prefers reviewed entries but falls back to unreviewed (TrEMBL).

    Args:
        gene_name: Human gene symbol (e.g. "ERBB2").

    Returns:
        Tuple of (accession, sequence) or (None, None) if not found.
    """
    cache_path = CACHE_DIR / "cyno_sequences" / "{}.json".format(gene_name)
    if cache_path.exists():
        with open(cache_path) as f:
            cached = json.load(f)
            return cached.get("accession"), cached.get("sequence")

    # Try cynomolgus (fascicularis) first, then rhesus (mulatta) as fallback.
    # Rhesus and cyno are >99% identical genome-wide and routinely
    # substituted for cross-reactive antibody assessments.
    from epitope_pipeline.config import RHESUS_TAXID
    species_ids = [
        (CYNO_TAXID, "macaca_fascicularis"),
        (RHESUS_TAXID, "macaca_mulatta"),
    ]

    for taxid, species_name in species_ids:
        # Try reviewed first, then unreviewed
        for reviewed in [True, False]:
            query = "gene:{} AND organism_id:{}".format(gene_name, taxid)
            if reviewed:
                query += " AND reviewed:true"
            url = "{}/uniprotkb/search".format(UNIPROT_API)
            params = {
                "query": query,
                "format": "json",
                "size": 3,
                "fields": "accession,gene_names,sequence",
            }
            try:
                data = _api_get(
                    url,
                    "Cyno ortholog search for '{}'".format(gene_name),
                    params=params,
                )
            except TargetResolutionError:
                continue

            results = data.get("results", [])
            if results:
                entry = results[0]
                accession = entry["primaryAccession"]
                sequence = entry.get("sequence", {}).get("value", "")

                if taxid != CYNO_TAXID:
                    logger.info(
                        "  Using %s ortholog %s as cyno proxy (>99%% genome-wide identity)",
                        species_name, accession,
                    )

                # Cache
                cache_path.parent.mkdir(parents=True, exist_ok=True)
                with open(cache_path, "w") as f:
                    json.dump({"accession": accession, "sequence": sequence}, f)

                return accession, sequence

    # Fallback: Ensembl ortholog lookup (try both species)
    logger.info("  UniProt cyno search failed, trying Ensembl ortholog lookup...")
    for taxid, species_name in species_ids:
        ensembl_result = _resolve_cyno_via_ensembl(gene_name, target_species=species_name)
        if ensembl_result[0] is not None:
            accession, sequence = ensembl_result
            if species_name != "macaca_fascicularis":
                logger.info("  Using %s Ensembl ortholog as cyno proxy", species_name)
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_path, "w") as f:
                json.dump({"accession": accession, "sequence": sequence}, f)
            return accession, sequence

    # Nothing found
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump({"accession": None, "sequence": None}, f)
    return None, None


def _resolve_cyno_via_ensembl(gene_name, target_species="macaca_fascicularis"):
    """
    Resolve macaque ortholog via Ensembl REST API.

    Uses the homology endpoint to find orthologs, then fetches the
    protein sequence for the best match.

    Args:
        gene_name: Human gene symbol (e.g. "CEACAM5").
        target_species: Ensembl species name (default: macaca_fascicularis).

    Returns:
        Tuple of (ensembl_gene_id, sequence) or (None, None).
    """
    ENSEMBL_REST = "https://rest.ensembl.org"

    try:
        # Step 1: Find orthologs
        url = "{}/homology/symbol/homo_sapiens/{}".format(ENSEMBL_REST, gene_name)
        params = {
            "target_species": target_species,
            "type": "orthologues",
            "content-type": "application/json",
        }
        resp = requests.get(url, params=params, timeout=30)
        if resp.status_code != 200:
            logger.debug("  Ensembl ortholog lookup returned %d", resp.status_code)
            return None, None

        data = resp.json()
        homologies = data.get("data", [{}])[0].get("homologies", [])
        if not homologies:
            logger.debug("  Ensembl: no cyno orthologs found for %s", gene_name)
            return None, None

        # Take the best match (highest identity)
        best = max(homologies, key=lambda h: h.get("target", {}).get("perc_id", 0))
        target_info = best.get("target", {})
        protein_id = target_info.get("protein_id", "")
        gene_id = target_info.get("id", "")
        identity = target_info.get("perc_id", 0)
        homology_type = best.get("type", "")

        if not protein_id:
            logger.debug("  Ensembl: no protein_id for cyno ortholog")
            return None, None

        logger.info(
            "  Ensembl cyno ortholog: %s (%s, %.1f%% identity, %s)",
            gene_id, protein_id, identity, homology_type,
        )

        # Step 2: Fetch protein sequence
        seq_url = "{}/sequence/id/{}".format(ENSEMBL_REST, protein_id)
        seq_resp = requests.get(
            seq_url,
            params={"content-type": "application/json"},
            timeout=30,
        )
        if seq_resp.status_code != 200:
            logger.warning("  Ensembl: failed to fetch sequence for %s", protein_id)
            return None, None

        sequence = seq_resp.json().get("seq", "")
        if not sequence:
            return None, None

        # Use Ensembl gene ID as the accession identifier
        return gene_id, sequence

    except requests.RequestException as e:
        logger.warning("  Ensembl ortholog lookup failed: %s", e)
        return None, None


# ---------------------------------------------------------------------------
# Entry parsing
# ---------------------------------------------------------------------------

def _parse_uniprot_entry(entry):
    """
    Parse a UniProt JSON entry into a TargetInfo dataclass.

    Extracts: accession, gene name, protein name, sequence, features
    (transmembrane, topological domain, etc.), and PDB cross-references.
    """
    accession = entry.get("primaryAccession", "")

    # Gene name: take the first from the first gene entry
    genes = entry.get("genes", [])
    gene_name = ""
    if genes:
        gene_name = genes[0].get("geneName", {}).get("value", "")

    # Protein name: recommended name > submitted name
    protein_desc = entry.get("proteinDescription", {})
    rec_name = protein_desc.get("recommendedName", {})
    protein_name = rec_name.get("fullName", {}).get("value", "")
    if not protein_name:
        sub_names = protein_desc.get("submissionNames", [])
        if sub_names:
            protein_name = sub_names[0].get("fullName", {}).get("value", "")

    # Organism
    organism = entry.get("organism", {}).get("scientificName", "")

    # Sequence
    seq_info = entry.get("sequence", {})
    sequence = seq_info.get("value", "")

    # Features: extract topology-relevant annotations
    features = []
    for feat in entry.get("features", []):
        feat_type = feat.get("type", "")
        if feat_type in (
            "Transmembrane", "Topological domain", "Intramembrane",
            "Signal peptide", "Signal", "Transit peptide", "Chain",
            "Domain", "Region", "Lipidation",
        ):
            location = feat.get("location", {})
            start = location.get("start", {}).get("value")
            end = location.get("end", {}).get("value")
            description = feat.get("description", "")
            features.append({
                "type": feat_type,
                "start": start,
                "end": end,
                "description": description,
            })

    # PDB cross-references
    pdb_ids = _extract_pdb_xrefs(entry)

    # InterPro domain annotations (ECD subdomains etc.)
    interpro_domains = _fetch_interpro_domains(accession)
    features.extend(interpro_domains)

    # Compute ectodomain ranges from topology features
    ec_ranges = _extract_ectodomain_ranges(features, len(sequence))
    ec_length = sum(end - start + 1 for start, end in ec_ranges)

    return TargetInfo(
        uniprot_id=accession,
        gene_name=gene_name,
        protein_name=protein_name,
        organism=organism,
        sequence=sequence,
        features=features,
        pdb_ids=pdb_ids,
        ectodomain_ranges=ec_ranges,
        ectodomain_length=ec_length,
    )


def _extract_ectodomain_ranges(features, sequence_length):
    """
    Extract extracellular residue ranges from UniProt topology features.

    For TM proteins: collects "Topological domain" features with
    "Extracellular" in their description.

    For GPI-anchored proteins (no Topological domain annotations):
    infers ectodomain as signal peptide end+1 through GPI anchor.

    Returns list of (start, end) tuples (1-based, inclusive).
    """
    ranges = []
    for feat in features:
        if feat.get("type") == "Topological domain":
            desc = feat.get("description", "")
            if "extracellular" in desc.lower():
                start = feat.get("start")
                end = feat.get("end")
                if start is not None and end is not None:
                    ranges.append((int(start), int(end)))

    if ranges:
        ranges.sort(key=lambda x: x[0])
        return ranges

    # Fallback for GPI-anchored proteins: no Topological domain features.
    # Ectodomain = after signal peptide to GPI anchor (or sequence end).
    sp_end = 0
    gpi_residue = None
    for feat in features:
        ftype = feat.get("type", "")
        if ftype in ("Signal peptide", "Signal") and feat.get("end"):
            sp_end = int(feat["end"])
        if ftype == "Lipidation" and "gpi" in feat.get("description", "").lower():
            if feat.get("start"):
                gpi_residue = int(feat["start"])

    ec_start = sp_end + 1 if sp_end > 0 else 1
    ec_end = gpi_residue if gpi_residue else sequence_length
    ranges.append((ec_start, ec_end))
    return ranges


def _extract_pdb_xrefs(entry):
    """
    Extract PDB IDs from UniProt cross-references.

    Returns a list of PDB ID strings (e.g. ["1N8Z", "3BE1", ...]).
    """
    pdb_ids = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "PDB":
            pdb_id = xref.get("id", "")
            if pdb_id:
                pdb_ids.append(pdb_id)
    return pdb_ids


# ---------------------------------------------------------------------------
# InterPro domain annotations
# ---------------------------------------------------------------------------

def _fetch_interpro_domains(uniprot_id):
    """
    Fetch domain-type InterPro annotations for a protein.

    Returns a list of feature dicts with type="InterPro" that can be
    appended to the target's feature list.  Only "domain" type entries
    are included (not family, superfamily, or site entries).

    Results are cached to cache/interpro/{accession}.json.
    """
    cache_path = CACHE_DIR / "interpro" / "{}.json".format(uniprot_id)
    if cache_path.exists():
        logger.debug("  InterPro cache hit: %s", cache_path)
        with open(cache_path) as f:
            return json.load(f)

    url = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{}/" \
          "?page_size=200&format=json".format(uniprot_id)

    features = []
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code != 200:
            logger.warning("  InterPro API returned %d for %s — skipping",
                           resp.status_code, uniprot_id)
            return features

        data = resp.json()
        # Deduplicate by (start, end) to avoid overlapping InterPro entries
        # for the same region (e.g. Protein kinase from both UniProt and InterPro)
        seen_ranges = set()

        # First pass: collect domain entries
        family_entries = []
        for result in data.get("results", []):
            metadata = result.get("metadata", {})
            entry_type = metadata.get("type", "")

            if entry_type == "domain":
                name = metadata.get("name", "")
                proteins = result.get("proteins", [])
                for prot in proteins:
                    for loc_group in prot.get("entry_protein_locations", []):
                        for fragment in loc_group.get("fragments", []):
                            start = fragment.get("start")
                            end = fragment.get("end")
                            if start is None or end is None:
                                continue
                            key = (start, end)
                            if key in seen_ranges:
                                continue
                            seen_ranges.add(key)
                            features.append({
                                "type": "InterPro",
                                "start": start,
                                "end": end,
                                "description": name,
                            })
            elif entry_type == "family":
                # Stash family entries as fallback
                family_entries.append(result)

        # Fallback: if no domain entries found, use family entries
        # (e.g. MSLN only has InterPro family annotations)
        if not features and family_entries:
            for result in family_entries:
                metadata = result.get("metadata", {})
                name = metadata.get("name", "")
                proteins = result.get("proteins", [])
                for prot in proteins:
                    for loc_group in prot.get("entry_protein_locations", []):
                        for fragment in loc_group.get("fragments", []):
                            start = fragment.get("start")
                            end = fragment.get("end")
                            if start is None or end is None:
                                continue
                            key = (start, end)
                            if key in seen_ranges:
                                continue
                            seen_ranges.add(key)
                            features.append({
                                "type": "InterPro",
                                "start": start,
                                "end": end,
                                "description": name,
                            })
            logger.info("  InterPro: %d family annotations (no domains) for %s",
                        len(features), uniprot_id)
        else:
            logger.info("  InterPro: %d domain annotations for %s", len(features), uniprot_id)
    except requests.RequestException as e:
        logger.warning("  InterPro fetch failed for %s: %s", uniprot_id, e)

    # Cache
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(features, f, indent=2)

    return features


# ---------------------------------------------------------------------------
# HTTP helper with retry and backoff
# ---------------------------------------------------------------------------

def _api_get(url, description, params=None):
    """
    GET request with retry logic for UniProt API.

    Args:
        url: Full URL to fetch.
        description: Human-readable description for error messages.
        params: Optional query parameters dict.

    Returns:
        Parsed JSON response.

    Raises:
        TargetResolutionError: After all retries exhausted.
    """
    max_attempts = RETRY_UNIPROT["attempts"]
    backoff = RETRY_UNIPROT["backoff_seconds"]

    for attempt in range(1, max_attempts + 1):
        try:
            resp = requests.get(url, params=params, timeout=30)
            if resp.status_code == 200:
                return resp.json()
            elif resp.status_code in (429, 503):
                # Rate limited or service unavailable — retry
                wait = backoff * attempt
                logger.warning(
                    "  %s returned %d, retrying in %ds (attempt %d/%d)",
                    description, resp.status_code, wait, attempt, max_attempts,
                )
                time.sleep(wait)
            else:
                raise TargetResolutionError(
                    "{} failed: HTTP {} — {}".format(
                        description, resp.status_code, resp.text[:200]
                    )
                )
        except requests.RequestException as e:
            if attempt < max_attempts:
                wait = backoff * attempt
                logger.warning(
                    "  %s request error: %s, retrying in %ds",
                    description, e, wait,
                )
                time.sleep(wait)
            else:
                raise TargetResolutionError(
                    "{} failed after {} attempts: {}".format(
                        description, max_attempts, e
                    )
                )

    raise TargetResolutionError(
        "{} failed after {} retries (rate limited)".format(description, max_attempts)
    )
