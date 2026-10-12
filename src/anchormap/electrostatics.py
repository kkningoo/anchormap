"""PROPKA output: what is usable and what is not.

The isoelectric point and individual shifted pKa values are usable. The
free-energy-of-folding curve is an empirical electrostatic sum over every
titratable residue, not a thermodynamic quantity, which is why values run into
the hundreds and thousands. Only its shape is interpretable, and it predicts
electrostatic stability rather than catalytic activity.
"""

from __future__ import annotations

from dataclasses import dataclass

MODEL_PKA = {'ASP': 3.80, 'GLU': 4.50, 'HIS': 6.50, 'CYS': 9.00,
             'TYR': 10.00, 'LYS': 10.50, 'ARG': 12.50,
             'N+': 8.00, 'C-': 3.20}

PHYSICAL_RANGE = (0.0, 14.0)


@dataclass(frozen=True)
class Titratable:
    restype: str
    number: str
    chain: str
    pka: float

    @property
    def shift(self) -> float:
        return self.pka - MODEL_PKA[self.restype]

    @property
    def physical(self) -> bool:
        return PHYSICAL_RANGE[0] <= self.pka <= PHYSICAL_RANGE[1]


def parse(path: str) -> dict:
    rows, curve = [], []
    pi_folded = pi_unfolded = opt_ph = None
    mode = None
    for line in open(path):
        s = line.strip()
        if s.startswith('SUMMARY OF THIS PREDICTION'):
            mode = 'sum'
            continue
        if s.startswith('Free energy of'):
            mode = 'curve'
            continue
        if s.startswith('The pH of optimum stability'):
            mode = None
            try:
                opt_ph = float(s.split('is')[1].split('for')[0])
            except (IndexError, ValueError):
                pass
            continue
        if s.startswith('The pI'):
            parts = s.replace('(folded)', '').replace('(unfolded)', '').split()
            nums = [p for p in parts if p.replace('.', '', 1).isdigit()]
            if len(nums) >= 2:
                pi_folded, pi_unfolded = float(nums[0]), float(nums[1])
            continue
        if mode == 'sum':
            p = s.split()
            if len(p) >= 4 and p[0] in MODEL_PKA:
                try:
                    rows.append(Titratable(p[0], p[1], p[2], float(p[3])))
                except ValueError:
                    pass
        elif mode == 'curve':
            p = s.split()
            if len(p) == 2:
                try:
                    curve.append((float(p[0]), float(p[1])))
                except ValueError:
                    pass

    return {'titratable': rows, 'n_titratable': len(rows), 'curve': curve,
            'pI_folded': pi_folded, 'pI_unfolded': pi_unfolded,
            'pH_optimum': opt_ph}


def shifted(result: dict, min_shift: float = 1.5) -> list[Titratable]:
    """Residues whose pKa departs from the model value, physical ones only.

    A depressed lysine pKa is a signature of a positively charged pocket; one
    depressed by 3.45 units independently supported a polyanion channel
    assignment here.
    """
    out = [t for t in result['titratable']
           if t.physical and abs(t.shift) >= min_shift]
    return sorted(out, key=lambda t: -abs(t.shift))


def unphysical(result: dict) -> list[Titratable]:
    """Values outside 0-14, which must be discarded rather than reported.

    PROPKA extrapolates without bound for buried charges lacking a
    counter-charge. A metal-stripped run here returned +22.0; the same
    structure with metals returned -5.62. Both are meaningless, and the
    presence of either is a signal that the run needs metals in place.
    """
    return [t for t in result['titratable'] if not t.physical]


def normalised_curve(result: dict) -> list[tuple[float, float]]:
    """Curve relative to the protein's own optimum, per titratable residue.

    This normalisation is a construction, not a PROPKA output. It exists only
    to make proteins of very different size comparable, and the shape is the
    only interpretable part. Do not present the absolute values.
    """
    curve, n = result['curve'], result['n_titratable']
    if not curve or not n:
        return []
    lo = min(dg for _, dg in curve)
    return [(ph, (dg - lo) / n) for ph, dg in curve]
