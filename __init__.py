"""
Epitope Pipeline — Structural bioinformatics pipeline for identifying
druggable VHH epitope space on human membrane protein targets.

Identifies ectodomain surface patches that are:
  (a) >=150 Angstroms from the membrane plane
  (b) >98% conserved with cynomolgus monkey
  (c) Unique to the target within the human proteome (<95% identity)
"""

__version__ = "0.1.0"
