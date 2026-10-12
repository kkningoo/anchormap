"""Candidate coupling sites, exclusions, and multipoint reachability."""

from __future__ import annotations

import itertools
import math
from collections import defaultdict
from dataclasses import dataclass

import freesasa

from .constants import (ADDITIVES, AVAILABILITY_RSA, CROSSLINKER_SPAN,
                        TIEN_THEORETICAL)
from .guards import check_selection_size, warn_crystallisation_additives


@dataclass(frozen=True)
class Site:
    chain: str
    number: str
    xyz: tuple[float, float, float]

    def distance_to(self, other: 'Site') -> float:
        return math.sqrt(sum((a - b) ** 2 for a, b in zip(self.xyz, other.xyz)))


def load_sites(pdb: str, resname: str = 'LYS', atom: str = 'NZ') -> dict:
    """Coupling-atom coordinates, keyed by (chain, resnum)."""
    out = {}
    for line in open(pdb):
        if (line.startswith('ATOM') and line[17:20] == resname
                and line[12:16].strip() == atom):
            key = (line[21], line[22:26].strip())
            out[key] = Site(key[0], key[1],
                            (float(line[30:38]), float(line[38:46]),
                             float(line[46:54])))
    return out


def available_sites(pdb: str, resname: str = 'LYS',
                    exclude: set[str] = frozenset()) -> dict:
    """Sites above the RSA threshold and not on the exclusion list."""
    sites = load_sites(pdb, resname)
    areas = freesasa.calc(freesasa.Structure(pdb)).residueAreas()
    maxacc = TIEN_THEORETICAL[resname]
    rsa = {(ch, n): a.total / maxacc
           for ch in areas for n, a in areas[ch].items()
           if a.residueType == resname}
    return {k: v for k, v in sites.items()
            if rsa.get(k, 0) > AVAILABILITY_RSA and k[1] not in exclude}


def reachable_pairs(sites: dict, crosslinker: str | None = None,
                    span: tuple[float, float] | None = None) -> list:
    """Pairs within the crosslinker span.

    One of `crosslinker` or `span` must be given; there is no default. The
    span changes the answer materially -- a 10 A floor rather than 5 A
    excluded four of six real cross-ring contacts on a GS dodecamer -- so the
    window is always stated by the caller and recorded with the result.
    """
    if (crosslinker is None) == (span is None):
        raise ValueError(
            'give exactly one of crosslinker= or span=. The window is not '
            'optional: it changes the counts and must be reported with them.'
        )
    lo, hi = span if span else CROSSLINKER_SPAN[crosslinker]
    out = []
    for a, b in itertools.combinations(sorted(sites), 2):
        d = sites[a].distance_to(sites[b])
        if lo <= d <= hi:
            out.append((a, b, d, a[0] != b[0]))
    return out


def coverage(sites: dict, pairs: list) -> dict:
    """How many sites have at least one partner, and the partner histogram.

    Multipoint attachment, not any point mutation, is what extends immobilised
    lifetime, by restricting the unfolding pathway. A site with no partner
    contributes attachment but not stabilisation.
    """
    nbrs = defaultdict(set)
    for a, b, _, _ in pairs:
        nbrs[a].add(b)
        nbrs[b].add(a)

    hist = defaultdict(int)
    per_chain = defaultdict(int)
    for k in sites:
        n = len(nbrs[k])
        hist[n] += 1
        if n:
            per_chain[k[0]] += 1

    return {
        'sites': len(sites),
        'with_partner': sum(1 for k in sites if nbrs[k]),
        'isolated': [k for k in sites if not nbrs[k]],
        'histogram': dict(sorted(hist.items())),
        'max_partners': max(hist) if hist else 0,
        'per_chain': dict(sorted(per_chain.items())),
        'cross_chain_pairs': sum(1 for p in pairs if p[3]),
        'neighbours': nbrs,
    }


def clusters(sites: dict, pairs: list) -> list[list]:
    """Connected components of the reachability graph.

    A single large component means the crosslinker can network the whole
    assembly as one covalent unit, which is a stronger statement than a count
    of pairs.
    """
    adj = defaultdict(set)
    for a, b, _, _ in pairs:
        adj[a].add(b)
        adj[b].add(a)

    seen, comps = set(), []
    for k in sites:
        if k in seen:
            continue
        stack, comp = [k], []
        while stack:
            n = stack.pop()
            if n in seen:
                continue
            seen.add(n)
            comp.append(n)
            stack.extend(adj[n] - seen)
        comps.append(sorted(comp))
    return sorted(comps, key=len, reverse=True)


def stratify_by_group(pairs: list, group_of: dict[str, int]) -> dict:
    """Split cross-chain pairs by whether they cross a structural boundary.

    For a stacked-ring assembly, pass {chain: ring_number}. Returns within-
    group and across-group counts, which is what tells you whether a
    crosslinker reaches the interface that matters.
    """
    within = [p for p in pairs if p[3] and group_of[p[0][0]] == group_of[p[1][0]]]
    across = [p for p in pairs if p[3] and group_of[p[0][0]] != group_of[p[1][0]]]
    return {'within_group': within, 'across_group': across,
            'n_within': len(within), 'n_across': len(across)}


def span_sweep(sites: dict, group_of: dict[str, int],
               floors: tuple[float, ...] = (5.0, 10.0),
               uppers: tuple[float, ...] = (16, 18, 20, 22, 25, 30)) -> list[dict]:
    """Counts against both floor and upper bound.

    Reporting the sweep rather than one window is what makes the window
    dependence visible. A count that changes threefold with the floor is not
    a property of the structure.
    """
    rows = []
    for up in uppers:
        row = {'upper': up}
        for lo in floors:
            p = reachable_pairs(sites, span=(lo, up))
            s = stratify_by_group(p, group_of)
            row[f'within_{lo:g}'] = s['n_within']
            row[f'across_{lo:g}'] = s['n_across']
        rows.append(row)
    return rows


def active_site_residues(pdb: str, ligand_resnames: set[str],
                         chain: str | None = None, cutoff: float = 5.0,
                         max_expected_atoms: int = 500) -> dict:
    """Residues within `cutoff` of any ligand atom, by chain.

    Iterates over all atoms and deduplicates by residue. Selecting on CA alone
    misses partner residues whose side chains reach the ligand but whose CA
    does not, which nearly produced a false negative on an inter-subunit site
    here -- a 39:4 split reported as "chain A only".
    """
    warn_crystallisation_additives(ligand_resnames)
    clean = ligand_resnames - ADDITIVES

    lig, prot = [], []
    for line in open(pdb):
        if line[:6] not in ('ATOM  ', 'HETATM'):
            continue
        rn, ch, ri = line[17:20].strip(), line[21], line[22:26].strip()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if rn in clean and (chain is None or ch == chain):
            lig.append(xyz)
        elif line.startswith('ATOM'):
            prot.append((ch, ri, rn, xyz))

    check_selection_size(len(lig), max_expected_atoms, 'ligand selection')

    c2 = cutoff ** 2
    hits = {}
    for ch, ri, rn, p in prot:
        for q in lig:
            if sum((a - b) ** 2 for a, b in zip(p, q)) <= c2:
                hits[(ch, ri)] = rn
                break

    by_chain = defaultdict(list)
    for (ch, ri), rn in hits.items():
        by_chain[ch].append((ri, rn))
    return {ch: sorted(v, key=lambda x: int(x[0])) for ch, v in by_chain.items()}


def partner_chain_contribution(pdb: str, host_chain: str,
                               ligand_resnames: set[str],
                               cutoff: float = 5.0) -> dict:
    """How many residues each chain contributes to one chain's ligand set.

    The ratio decides whether monomer calculations are admissible. A 39:4
    split means they are not; an 88:6 atom split with one partner residue at
    the edge means they nearly are.
    """
    hits = active_site_residues(pdb, ligand_resnames, chain=host_chain,
                                cutoff=cutoff)
    counts = {ch: len(v) for ch, v in hits.items()}
    partners = {ch: v for ch, v in hits.items() if ch != host_chain}
    return {
        'counts': counts,
        'host': counts.get(host_chain, 0),
        'partners': partners,
        'inter_subunit': bool(partners),
        'verdict': ('inter-subunit: monomer calculations are not admissible'
                    if partners else
                    'self-contained: monomer calculations are admissible'),
    }
