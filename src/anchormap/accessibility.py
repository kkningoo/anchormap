"""Solvent accessibility, normalised, compared between oligomeric states.

The two-state comparison is the point. A map computed on a monomer lists sites
that do not exist on the assembled particle: five of twenty-one exposed GS
lysines cross below threshold on assembly, and one falls from 0.71 to 0.28.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from statistics import mean

import freesasa

from .constants import AVAILABILITY_RSA, THRESHOLD_MARGIN, TIEN_THEORETICAL
from .guards import require_minimised


@dataclass(frozen=True)
class ResidueRSA:
    chain: str
    number: str
    name: str
    rsa_subunit: float
    rsa_assembly: float

    @property
    def delta(self) -> float:
        return self.rsa_assembly - self.rsa_subunit

    @property
    def available(self) -> bool:
        return self.rsa_assembly > AVAILABILITY_RSA

    @property
    def lost_to_assembly(self) -> bool:
        return self.rsa_subunit > AVAILABILITY_RSA and not self.available

    @property
    def near_threshold(self) -> bool:
        """Chain-dependent on a large assembly; report as a range."""
        return abs(self.rsa_assembly - AVAILABILITY_RSA) < THRESHOLD_MARGIN


def extract_chain(src: str, dst: str, chain: str = 'A') -> str:
    """Polymer atoms of one chain, so FreeSASA sees an isolated subunit."""
    with open(dst, 'w') as out:
        for line in open(src):
            if line.startswith('ATOM') and line[21] == chain:
                out.write(line)
        out.write('END\n')
    return dst


def rsa_by_chain(pdb: str, chain: str = 'A') -> dict[str, tuple[str, float]]:
    """{resnum: (resname, rsa)} for one chain, Tien theoretical maxima.

    Values above 1.0 occur at free termini, where the Gly-X-Gly assumption
    does not hold. Reported as computed rather than clipped.
    """
    areas = freesasa.calc(freesasa.Structure(pdb)).residueAreas()
    out = {}
    for num, a in areas[chain].items():
        m = TIEN_THEORETICAL.get(a.residueType)
        if m:
            out[num] = (a.residueType, a.total / m)
    return out


def compare_states(assembly_pdb: str, chain: str = 'A',
                   subunit_pdb: str | None = None,
                   resname_filter: str | None = 'LYS',
                   allow_unminimised: bool = False) -> list[ResidueRSA]:
    """Shrake-Rupley on both oligomeric states, matched by residue number.

    If subunit_pdb is not given, the chain is extracted from the assembly.
    Filter to lysines by default, since that is the coupling chemistry; pass
    None for every residue.
    """
    require_minimised(assembly_pdb, allow=allow_unminimised)

    if subunit_pdb is None:
        subunit_pdb = f'_chain{chain}_only.pdb'
        if not os.path.exists(subunit_pdb):
            extract_chain(assembly_pdb, subunit_pdb, chain)

    asm = rsa_by_chain(assembly_pdb, chain)
    sub = rsa_by_chain(subunit_pdb, chain)

    out = []
    for num in sorted(set(asm) & set(sub), key=lambda n: int(n)):
        name, a = asm[num]
        _, s = sub[num]
        if resname_filter and name != resname_filter:
            continue
        out.append(ResidueRSA(chain, num, name, s, a))
    return out


def summarise(records: list[ResidueRSA]) -> dict:
    lost = [r for r in records if r.lost_to_assembly]
    return {
        'total': len(records),
        'available_subunit': sum(r.rsa_subunit > AVAILABILITY_RSA for r in records),
        'available_assembly': sum(r.available for r in records),
        'lost_to_assembly': lost,
        'near_threshold': [r for r in records if r.near_threshold],
    }


def chain_variation(assembly_pdb: str, resnum: str) -> dict:
    """One residue across every chain.

    Run this for any position within THRESHOLD_MARGIN of the cutoff before
    stating a single value. It also validates the assembly: a symmetric
    biological unit gives identical values for symmetry-related chains, and a
    merge or chain-renaming error breaks that pattern.
    """
    areas = freesasa.calc(freesasa.Structure(assembly_pdb)).residueAreas()
    vals = {}
    for ch in sorted(areas):
        a = areas[ch].get(resnum)
        if a:
            m = TIEN_THEORETICAL.get(a.residueType)
            if m:
                vals[ch] = a.total / m
    if not vals:
        return {'residue': resnum, 'found': False}
    v = list(vals.values())
    return {
        'residue': resnum,
        'found': True,
        'per_chain': vals,
        'min': min(v),
        'max': max(v),
        'spread': max(v) - min(v),
        'below_threshold': sum(1 for x in v if x <= AVAILABILITY_RSA),
        'n_chains': len(v),
        'chain_dependent': min(v) <= AVAILABILITY_RSA < max(v),
    }


def terminal_accessibility(assembly_pdb: str, chain: str,
                           n_term: list[str], c_term: list[str],
                           subunit_pdb: str | None = None) -> dict:
    """Terminal exposure in both states, for tag and fusion-point placement.

    Report the per-residue table, not a single mean. A mean over a terminal
    segment depends on how many residues you include and which maxima you
    normalise by, and two defensible choices here differ by twenty points.
    """
    recs = compare_states(assembly_pdb, chain, subunit_pdb, resname_filter=None)
    by_num = {r.number: r for r in recs}

    def block(nums):
        rows = [by_num[n] for n in nums if n in by_num]
        if not rows:
            return {}
        return {
            'residues': rows,
            'mean_subunit': mean(r.rsa_subunit for r in rows),
            'mean_assembly': mean(r.rsa_assembly for r in rows),
            'range_subunit': (min(r.rsa_subunit for r in rows),
                              max(r.rsa_subunit for r in rows)),
            'range_assembly': (min(r.rsa_assembly for r in rows),
                               max(r.rsa_assembly for r in rows)),
        }

    return {'n_terminus': block(n_term), 'c_terminus': block(c_term)}
