"""Chemical degradation liabilities on the assembled state.

A statement about position, not about consequence. No rates are estimated, so
these cannot be ranked against the physical inactivation routes -- subunit
dissociation, leaching, unfolding -- that the interface analysis addresses.
"""

from __future__ import annotations

from dataclasses import dataclass

from .constants import (DEAMIDATION_FOLLOWERS, ISOMERISATION_FOLLOWERS,
                        LIABILITY_RSA, OXIDISABLE)
from .accessibility import rsa_by_chain

THREE2ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E',
    'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V',
}


@dataclass(frozen=True)
class Liability:
    kind: str
    motif: str
    resnum: str
    rsa: float
    in_critical_set: bool


def chain_sequence(pdb: str, chain: str) -> tuple[str, list[str]]:
    seq, nums, seen = [], [], set()
    for l in open(pdb):
        if l.startswith('ATOM') and l[21] == chain and l[12:16].strip() == 'CA':
            r = l[22:26].strip()
            if r in seen:
                continue
            seen.add(r)
            seq.append(THREE2ONE.get(l[17:20], 'X'))
            nums.append(r)
    return ''.join(seq), nums


def scan(pdb: str, chain: str = 'A', critical: set[str] = frozenset(),
         rsa_cut: float = LIABILITY_RSA) -> list[Liability]:
    """Exposed deamidation, isomerisation and oxidation sites.

    The RSA reported for a two-residue motif is that of the FIRST residue,
    which is the one that reacts. An Asn-Thr motif is therefore normalised by
    the asparagine maximum, not the threonine one -- getting this wrong is an
    easy way to shift a value by twenty percent.
    """
    seq, nums = chain_sequence(pdb, chain)
    rsa = rsa_by_chain(pdb, chain)
    found = []

    for i in range(len(seq) - 1):
        if seq[i] == 'N' and seq[i + 1] in DEAMIDATION_FOLLOWERS:
            found.append(('deamidation',
                          f'N{nums[i]}{seq[i+1]}{nums[i+1]}', nums[i]))
        if seq[i] == 'D' and seq[i + 1] in ISOMERISATION_FOLLOWERS:
            found.append(('isomerisation',
                          f'D{nums[i]}{seq[i+1]}{nums[i+1]}', nums[i]))
    for i, aa in enumerate(seq):
        if aa in OXIDISABLE:
            found.append(('oxidation', f'{aa}{nums[i]}', nums[i]))

    out = []
    for kind, motif, num in found:
        _, r = rsa.get(num, ('?', 0.0))
        if r >= rsa_cut:
            out.append(Liability(kind, motif, num, r, num in critical))
    return sorted(out, key=lambda x: -x.rsa)


def protection_by_assembly(assembly_pdb: str, chain: str,
                           subunit_pdb: str | None = None,
                           critical: set[str] = frozenset()) -> list[dict]:
    """Does assembly bury any liability?

    Usually not, for the fast motifs. An Asn-Gly that stays fully exposed on
    the particle is the residue most likely to limit shelf life, and assembly
    offers it nothing.
    """
    from .accessibility import compare_states
    recs = {r.number: r for r in
            compare_states(assembly_pdb, chain, subunit_pdb, resname_filter=None)}
    out = []
    for lia in scan(assembly_pdb, chain, critical):
        r = recs.get(lia.resnum)
        out.append({
            'liability': lia,
            'rsa_subunit': r.rsa_subunit if r else None,
            'rsa_assembly': lia.rsa,
            'delta': r.delta if r else None,
            'protected': bool(r and r.delta < -0.05),
        })
    return out
