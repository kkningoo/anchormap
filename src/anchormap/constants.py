"""Shared constants, with the provenance of each.

Keeping these in one module rather than duplicated per script is what stopped
three accessibility tables in this project from using three different maxima.
"""

from __future__ import annotations

# Tien et al. (2013) THEORETICAL maxima, Gly-X-Gly context, A^2.
# Use these everywhere. An earlier script here used MET 203, SER 130, THR 142
# -- values from a different table -- which inflated methionine RSA by ~10%
# and threonine by ~21% and made two tables in the same document
# incomparable.
TIEN_THEORETICAL = {
    'ALA': 129, 'ARG': 274, 'ASN': 195, 'ASP': 193, 'CYS': 167,
    'GLU': 223, 'GLN': 225, 'GLY': 104, 'HIS': 224, 'ILE': 197,
    'LEU': 201, 'LYS': 236, 'MET': 224, 'PHE': 240, 'PRO': 159,
    'SER': 155, 'THR': 172, 'TRP': 285, 'TYR': 263, 'VAL': 174,
}

# Availability threshold for a coupling site.
AVAILABILITY_RSA = 0.30

# Positions within this of the threshold are chain-dependent on a large
# assembly and should be reported as a range. Lys86 of B. subtilis GS spans
# 0.173-0.310 across twelve chains and sits below 0.30 in ten of them.
THRESHOLD_MARGIN = 0.05

# Crosslinker reach, Angstrom. The LOWER bound is NOT a property of the
# reagent: a bifunctional crosslinker is limited by how far it can span, not
# by how close two amines may be, and genipin self-polymerises into short
# dimers and trimers so its span is a distribution with short members. A 10 A
# floor applied here once excluded four of six real cross-ring contacts at
# 8.4-8.5 A. The floor is steric, not chemical.
CROSSLINKER_SPAN = {
    'genipin': (5.0, 16.0),
    'glutaraldehyde': (5.0, 12.0),
}

# Interfaces below this are usually crystallographic rather than biological.
BIOLOGICAL_INTERFACE_MIN_AREA = 800.0

# A post-relax total score this far above its siblings means the run has not
# converged, not that the interface is unfavourable.
CONVERGENCE_SCORE_MARGIN = 50.0

# Chemical liability motifs. The RSA reported for a two-residue motif is that
# of the FIRST residue, which is the one that reacts.
DEAMIDATION_FOLLOWERS = 'GSTNH'   # Asn followed by these
ISOMERISATION_FOLLOWERS = 'GST'   # Asp followed by these
OXIDISABLE = 'MCW'

LIABILITY_RSA = 0.20

# Crystallisation additives that contaminate ligand-proximity selections.
# Sulfates pulled two lysines into a GS active-site definition in 2 of 12
# chains here.
ADDITIVES = {'SO4', 'PO4', 'GOL', 'EDO', 'PEG', 'MPD', 'ACT', 'CL', 'NA', 'K',
             'CA', 'ZN', 'IOD', 'BR', 'FMT', 'DMS', 'TRS', 'EPE'}
