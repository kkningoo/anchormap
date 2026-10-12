"""Interface area and energetics, with convergence reported rather than hidden.

The question this answers for immobilisation is narrow: which interface opens
first, so crosslinker reach can be checked against it. Stitch count is not a
proxy for interface strength -- assuming it was inverted a reinforcement
strategy here.
"""

from __future__ import annotations

import itertools
import os
from dataclasses import dataclass
from statistics import median

import freesasa

from .constants import BIOLOGICAL_INTERFACE_MIN_AREA, CONVERGENCE_SCORE_MARGIN
from .requirements import require


@dataclass(frozen=True)
class InterfaceResult:
    pair: tuple[str, str]
    buried_area: float
    dG: float | None = None
    packstat: float | None = None
    hbond_energy: float | None = None
    n_residues: int | None = None
    score_pre_relax: float | None = None
    score_post_relax: float | None = None
    converged: bool = True
    backbone_free: bool = False

    @property
    def dG_per_1000A2(self) -> float | None:
        if self.dG is None or not self.buried_area:
            return None
        return 1000.0 * self.dG / self.buried_area

    @property
    def likely_biological(self) -> bool:
        return self.buried_area >= BIOLOGICAL_INTERFACE_MIN_AREA


def split_chains(pdb: str, chains: str, outdir: str = '.') -> dict[str, str]:
    """One file per chain, polymer only. Needed before pairwise areas."""
    paths = {}
    handles = {c: open(os.path.join(outdir, f'chain_{c}.pdb'), 'w') for c in chains}
    for line in open(pdb):
        if line.startswith('ATOM') and line[21] in handles:
            handles[line[21]].write(line)
    for c, h in handles.items():
        h.write('END\n')
        h.close()
        paths[c] = os.path.join(outdir, f'chain_{c}.pdb')
    return paths


def buried_areas(pdb: str, chains: str, outdir: str = '.') -> dict:
    """Buried surface per side for every chain pair, by FreeSASA.

    (SASA_a + SASA_b - SASA_ab) / 2. Pairs below the biological-relevance
    threshold are flagged rather than dropped, because a small interface is a
    result: a 296 A^2 reading is what identified a deposited "dimer" as the
    outlier in a family of tetramers.
    """
    require('freesasa')
    paths = split_chains(pdb, chains, outdir)
    total = lambda p: freesasa.calc(freesasa.Structure(p)).totalArea()
    single = {c: total(p) for c, p in paths.items()}

    out = {}
    for c1, c2 in itertools.combinations(chains, 2):
        pair_path = os.path.join(outdir, f'pair_{c1}{c2}.pdb')
        with open(pair_path, 'w') as f:
            for c in (c1, c2):
                f.write(open(paths[c]).read())
        bsa = (single[c1] + single[c2] - total(pair_path)) / 2
        out[(c1, c2)] = InterfaceResult(pair=(c1, c2), buried_area=bsa)
    return out


def analyse(pair_pdb: str, pair: tuple[str, str], relax_cycles: int = 5,
            backbone_free: bool = False) -> InterfaceResult:
    """Rosetta InterfaceAnalyzer. Relaxation is mandatory.

    Unrelaxed crystal poses returned roughly +550 REU here, dominated by clash
    repulsion. Both scores are returned so a reader can see the run converged.

    backbone_free changes the protocol and the numbers. Use it only to rescue
    a pair that will not converge with the backbone fixed, and never compare
    its output against fixed-backbone results.
    """
    require('pyrosetta')
    if relax_cycles < 1:
        raise ValueError(
            'relax_cycles must be at least 1. Unrelaxed poses return positive '
            'dG dominated by clash repulsion and are not interpretable.'
        )

    from pyrosetta import init, pose_from_pdb
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.core.scoring import get_score_function
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    from pyrosetta.rosetta.protocols.relax import FastRelax

    init('-mute all -ignore_unrecognized_res -ignore_zero_occupancy false')
    sf = get_score_function()
    pose = pose_from_pdb(pair_pdb)

    pre = sf(pose)
    mm = MoveMap()
    mm.set_bb(backbone_free)
    mm.set_chi(True)
    fr = FastRelax(sf, relax_cycles)
    fr.set_movemap(mm)
    fr.apply(pose)
    post = sf(pose)

    ia = InterfaceAnalyzerMover(1, False, sf, True, False, True, False)
    ia.apply(pose)
    d = ia.get_all_data()

    return InterfaceResult(
        pair=pair,
        buried_area=ia.get_interface_delta_sasa(),
        dG=ia.get_interface_dG(),
        packstat=d.packstat,
        hbond_energy=d.total_hb_E,
        n_residues=ia.get_num_interface_residues(),
        score_pre_relax=pre,
        score_post_relax=post,
        backbone_free=backbone_free,
    )


def flag_nonconverged(results: list[InterfaceResult]) -> list[InterfaceResult]:
    """Runs whose post-relax score sits far above their siblings.

    Two interfaces in this project did not converge under the fixed-backbone
    protocol and are reported separately rather than ranked with the rest.
    """
    scores = [r.score_post_relax for r in results if r.score_post_relax is not None]
    if len(scores) < 3:
        return []
    med = median(scores)
    return [r for r in results
            if r.score_post_relax is not None
            and r.score_post_relax - med > CONVERGENCE_SCORE_MARGIN]


def rank(results: list[InterfaceResult]) -> list[InterfaceResult]:
    """Weakest first, by energy density.

    Per-area normalisation carries real weight when interfaces differ in size,
    so report buried area alongside any ranking that uses it. Mixing protocols
    is refused: a backbone-free result is not comparable to a fixed-backbone
    one.
    """
    usable = [r for r in results if r.converged and r.dG is not None]
    protocols = {r.backbone_free for r in usable}
    if len(protocols) > 1:
        raise ValueError(
            'results mix fixed-backbone and backbone-free runs, which are not '
            'comparable. Rank them separately and say which protocol each used.'
        )
    return sorted(usable, key=lambda r: r.dG_per_1000A2, reverse=True)


def crosslinker_covers_weakest(ranked: list[InterfaceResult],
                               pairs_by_interface: dict) -> dict:
    """Does the crosslinker reach the interfaces most likely to open?

    This is the question the interface analysis exists to answer. Here the
    answer was yes, by coincidence rather than design: genipin's within-ring
    stitches sat at the weaker joints.
    """
    return {r.pair: len(pairs_by_interface.get(r.pair, [])) for r in ranked}
