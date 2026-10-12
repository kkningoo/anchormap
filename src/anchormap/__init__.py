"""Where can an assembled multimeric enzyme be attached to a solid support?

anchormap takes coordinates of an oligomeric protein and answers the questions
that decide an immobilisation protocol: which surface sites exist on the
assembly rather than the subunit, which would destroy activity if coupled,
whether multipoint attachment is geometrically possible, which subunit
interface opens first, what degrades first, and whether the support chemistry
is compatible with the enzyme and its substrates.

It does not predict structures, dock ligands, or run dynamics. OpenMM,
FreeSASA, Rosetta, Reduce and MolProbity are wrapped, not reimplemented.

Three of the external tools are not pip-installable. Run
``python -m anchormap.requirements`` to see what is present; every wrapper
checks before it runs and names the tool it needs.

Order of use::

    prep.load_assembly -> prep.add_hydrogens_with_flips -> prep.minimise
                                   -> prep.validate
                                      |
          +-----------------+---------+---------+-----------------+
          |                 |                   |                 |
    accessibility      coupling            interfaces          support
    .compare_states    .active_site_       .buried_areas       .differential_
    .chain_variation    residues           .analyse             adsorption
    .terminal_         .available_sites    .rank               .minimum_pore_
     accessibility     .reachable_pairs    .flag_nonconverged   size
                       .coverage                               .polyelectrolyte
    liabilities        .clusters                                _barrier_risk
    .scan              .span_sweep
    .protection_by_    .partner_chain_     electrostatics
     assembly           contribution       .parse / .shifted

The preconditions in ``guards`` each correspond to a mistake that produced a
plausible wrong answer during development. They raise rather than warn.

Two interfaces are deliberately narrow for the same reason.
``coupling.reachable_pairs`` has no default span, because the window changes
the counts materially and must be reported with them. ``interfaces.rank``
refuses to mix fixed-backbone and backbone-free runs.
"""

from . import (accessibility, constants, coupling, electrostatics, guards,
               interfaces, liabilities, prep, requirements, support)

__all__ = [
    "accessibility", "constants", "coupling", "electrostatics", "guards",
    "interfaces", "liabilities", "prep", "requirements", "support",
]
__version__ = "0.2.0"
