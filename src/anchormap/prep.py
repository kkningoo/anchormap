"""Structure preparation. Everything else depends on it.

Refinement is a gate, not housekeeping: an unminimised build here produced a
coupling-site map that was wrong in both directions.
"""

from __future__ import annotations

from dataclasses import dataclass

from .guards import (check_biological_assembly, check_metals_present,
                     check_numbering, mark_minimised, verify_mutation)
from .requirements import require

RESTRAINT_K = 1000.0          # kJ mol^-1 nm^-2, heavy atoms
MAX_MIN_ITER = 2000
FORCEFIELD = 'amber14-all.xml'
SOLVENT = 'implicit/obc2.xml'


@dataclass
class PreparedStructure:
    path: str
    minimised: bool = False
    n_chains: int = 0
    n_metals: int = 0
    first_residue: int = 1
    energy_before: float | None = None
    energy_after: float | None = None

    def first_residue_number(self) -> int:
        return self.first_residue


def inspect(path: str) -> dict:
    """Chain count, metal count and first residue, without loading a toolkit."""
    chains, metals, first = set(), 0, None
    known_metals = {'MG', 'MN', 'ZN', 'CA', 'FE', 'NI', 'CO', 'CU'}
    for line in open(path):
        if line.startswith('ATOM'):
            chains.add(line[21])
            if first is None:
                first = int(line[22:26])
        elif line.startswith('HETATM') and line[17:20].strip() in known_metals:
            metals += 1
    return {'n_chains': len(chains), 'chains': sorted(chains),
            'n_metals': metals, 'first_residue': first or 1}


def load_assembly(path: str, expected_chains: int, expected_metals: int,
                  source_first_residue: int) -> PreparedStructure:
    """Load a biological assembly, with the three checks that caught errors.

    Use the .pdb1 assembly file, not the asymmetric unit, and rename chains on
    merge. Keep metals: stripping them changes electrostatic output by a
    factor of two and makes pKa values unusable.
    """
    info = inspect(path)
    s = PreparedStructure(path=path, n_chains=info['n_chains'],
                          n_metals=info['n_metals'],
                          first_residue=info['first_residue'])
    check_biological_assembly(s.n_chains, expected_chains)
    check_metals_present(s.n_metals, expected_metals)
    check_numbering(s, source_first_residue)
    return s


def minimise(structure: PreparedStructure, restraint_k: float = RESTRAINT_K,
             out_path: str | None = None, ph: float = 7.0) -> PreparedStructure:
    """Restrained minimisation in OpenMM.

    Restrained rather than free: free minimisation relaxes the structure
    toward the force field's preference and away from the coordinates being
    refined, which defeats the purpose for surface calculations.

    Expect small regressions in bond angles, rotamers and Cbeta deviation --
    the cost of the restraints. An unblemished improvement means they were
    not binding.

    Two warnings from experience. Ligands and metals must usually be stripped
    first, because the force field has no parameters for them; coordination
    distances then lengthen, so a wild-type-versus-mutant comparison stays
    valid while absolute geometry does not. And a full large assembly may not
    finish: a 5,316-residue dodecamer in implicit solvent with no cutoff ran
    nine hours here without converging and was replaced by a dimer carrying
    one complete inter-subunit site.
    """
    require('openmm', 'pdbfixer')
    from openmm import CustomExternalForce, LangevinIntegrator, unit
    from openmm.app import ForceField, Modeller, NoCutoff, PDBFile, Simulation

    pdb = PDBFile(structure.path)
    ff = ForceField(FORCEFIELD, SOLVENT)
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(ff, pH=ph)

    system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff,
                             constraints=None)

    force = CustomExternalForce('k*periodicdistance(x,y,z,x0,y0,z0)^2')
    force.addGlobalParameter(
        'k', restraint_k * unit.kilojoules_per_mole / unit.nanometer ** 2)
    for p in ('x0', 'y0', 'z0'):
        force.addPerParticleParameter(p)
    for atom in modeller.topology.atoms():
        if atom.element.symbol != 'H':
            pos = modeller.positions[atom.index]
            force.addParticle(atom.index, [pos.x, pos.y, pos.z])
    system.addForce(force)

    integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                    0.002 * unit.picoseconds)
    sim = Simulation(modeller.topology, system, integrator)
    sim.context.setPositions(modeller.positions)

    before = sim.context.getState(getEnergy=True).getPotentialEnergy()
    sim.minimizeEnergy(maxIterations=MAX_MIN_ITER)
    after = sim.context.getState(getEnergy=True).getPotentialEnergy()

    out_path = out_path or structure.path.replace('.pdb', '_min.pdb')
    state = sim.context.getState(getPositions=True)
    with open(out_path, 'w') as f:
        PDBFile.writeFile(modeller.topology, state.getPositions(), f)

    mark_minimised(out_path)
    return PreparedStructure(
        path=out_path, minimised=True, n_chains=structure.n_chains,
        n_metals=structure.n_metals, first_residue=structure.first_residue,
        energy_before=before.value_in_unit(unit.kilojoule_per_mole),
        energy_after=after.value_in_unit(unit.kilojoule_per_mole),
    )


def add_hydrogens_with_flips(path: str, out_path: str | None = None) -> str:
    """Reduce, with Asn/Gln/His flip correction.

    OpenMM's addHydrogens does not do flips. Run this first if flip correction
    matters: four His in one tetramer here needed flipping, and hydrogens plus
    flips alone improved the clashscore by 0.71 before minimisation
    contributed anything.
    """
    require('reduce')
    import subprocess
    out_path = out_path or path.replace('.pdb', '_H.pdb')
    with open(out_path, 'w') as f:
        subprocess.run(['reduce', '-FLIP', path], stdout=f, check=True)
    return out_path


def validate(path: str) -> dict:
    """MolProbity, before and after. Report both.

    The B-factor filter is inverted on predicted models: the B column holds
    pLDDT, where high is good, so a low-B selection retains the worst
    residues. A clashscore of zero from such a filter is an empty selection.
    """
    require('molprobity')
    raise NotImplementedError(
        'Port your MolProbity invocation here. Return a dict with clashscore, '
        'molprobity_score, rama_favoured, rama_outliers, poor_rotamers, '
        'cbeta_deviations, bad_bonds and bad_angles, and report the same keys '
        'before and after refinement.'
    )
