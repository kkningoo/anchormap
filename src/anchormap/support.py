"""Support compatibility: will the enzymes load, and will substrate reach them?

Charge and size, not kinetics. These results reach the device models directly
without passing through any rate law.
"""

from __future__ import annotations

from dataclasses import dataclass

VDW_ALLOWANCE = 3.0  # Angstrom, added to the limiting extent for pore sizing


@dataclass(frozen=True)
class AdsorptionPrediction:
    enzyme: str
    pI: float
    operating_pH: float
    support_charge: str

    @property
    def net_charge(self) -> str:
        return "negative" if self.operating_pH > self.pI else "positive"

    @property
    def adsorbs_readily(self) -> bool:
        return (
            (self.support_charge == "cationic" and self.net_charge == "negative")
            or (self.support_charge == "anionic" and self.net_charge == "positive")
        )


def differential_adsorption(predictions: list[AdsorptionPrediction]) -> dict:
    """Will a mixed-solution load give the intended ratio?

    If the enzymes carry opposite net charge against a charged support, the
    ratio in the bead is not the ratio in the coupling buffer. Direction is
    predictable; magnitude is not, and needs a loading experiment.
    """
    readily = [p.enzyme for p in predictions if p.adsorbs_readily]
    repelled = [p.enzyme for p in predictions if not p.adsorbs_readily]
    return {
        "favoured": readily,
        "disfavoured": repelled,
        "ratio_biased": bool(readily) and bool(repelled),
        "note": "Direction known, magnitude unknown. Measure the achieved "
                "ratio by densitometry or paired activity assay.",
    }


def minimum_pore_size(extents: tuple[float, float, float]) -> float:
    """Pore diameter required for internal loading of an assembled particle.

    A cylinder enters a pore end-on more easily than side-on, so the limiting
    dimension is the smallest extent, not the largest.
    """
    return min(extents) + VDW_ALLOWANCE


def polyelectrolyte_barrier_risk(substrate_charge: str, support_charge: str) -> str | None:
    """Flag a self-reinforcing transport barrier.

    A polyanionic substrate meeting a polycationic support ionically
    crosslinks the outer shell as it arrives, tightening it against the next
    chain. The barrier strengthens in proportion to substrate supply, so the
    enzyme starves amid substrate. A Thiele analysis treating the substrate as
    a neutral diffusing species will not show this.
    """
    opposite = {("polyanionic", "cationic"), ("polycationic", "anionic")}
    if (substrate_charge, support_charge) in opposite:
        return (
            "Self-reinforcing barrier likely: substrate gels the support on "
            "contact. Measure uptake against cumulative exposure, not at a "
            "single time point."
        )
    return None
