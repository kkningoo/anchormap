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


def parse_molprobity(path: str) -> dict:
    """Read a saved MolProbity summary into a dict.

    MolProbity has no stable programmatic interface, so the run itself is
    done by hand at molprobity.biochem.duke.edu and the summary saved. This
    reads that file, whether saved as text or as the HTML summary page, so
    the numbers reach a report without being retyped.

    Retyping is the failure this avoids: a validation table is long, every
    value has a percentage beside it, and a transposed digit in a clashscore
    is invisible to a reader.

    Returns the metrics it finds; keys absent from the file are absent from
    the dict rather than guessed.
    """
    import re

    text = open(path, encoding="utf-8", errors="replace").read()
    # Strip HTML tags only. A bare "<" or ">" is ordinary MolProbity text --
    # goals are written "Goal: <0.3%" and "Cbeta deviations >0.25A" -- and a
    # greedy <...> pattern deletes everything between them, silently removing
    # the Ramachandran and MolProbity score rows.
    text = re.sub(r"</?[A-Za-z!][^>]*>", " ", text)
    text = re.sub(r"&[a-z]+;", " ", text)
    text = re.sub(r"[ \t]+", " ", text)

    num = r"(-?\d+(?:\.\d+)?)"
    fields = {
        "clashscore":          rf"Clashscore,?\s*all atoms:?\s*{num}",
        "molprobity_score":    rf"MolProbity score[^0-9\-]*{num}",
        "poor_rotamers":       rf"Poor rotamers\s*{num}",
        "favoured_rotamers":   rf"Favou?red rotamers\s*{num}",
        "rama_outliers":       rf"Ramachandran outliers\s*{num}",
        "rama_favoured":       rf"Ramachandran favou?red\s*{num}",
        "rama_z":              rf"Rama distribution Z-score[^0-9\-]*{num}",
        "cbeta_deviations":    rf"C\S* deviations[^0-9]*0\.25[^0-9]*{num}",
        "bad_bonds":           rf"Bad bonds:?\s*{num}",
        "bad_angles":          rf"Bad angles:?\s*{num}",
        "cablam_outliers":     rf"CaBLAM outliers\s*{num}",
        "ca_geometry_outliers": rf"CA Geometry outliers\s*{num}",
        "cis_prolines":        rf"Cis Prolines:?\s*{num}",
        "chiral_outliers":     rf"Chiral volume outliers\s*{num}",
    }

    out = {}
    for key, pattern in fields.items():
        m = re.search(pattern, text, re.I)
        if m:
            out[key] = float(m.group(1))

    if not out:
        raise ValueError(
            f"no MolProbity metrics found in {path}. Save the summary page "
            "or its text export, not a screenshot."
        )
    return out


def compare_validation(before: str, after: str) -> dict:
    """Before and after refinement, side by side.

    Report both. A single post-refinement table says the structure is good;
    the pair says the refinement did something, and shows what it cost. Minor
    regressions in bond angles, rotamers and C-beta deviation are expected
    under positional restraint and an unblemished improvement usually means
    the restraints were not binding.
    """
    b, a = parse_molprobity(before), parse_molprobity(after)
    rows = {}
    for k in sorted(set(b) | set(a)):
        rows[k] = {
            "before": b.get(k),
            "after": a.get(k),
            "change": (a[k] - b[k]) if k in a and k in b else None,
        }
    return rows


def validate(path: str) -> dict:
    """MolProbity cannot be driven from here; run it and parse the result.

    Submit the structure at molprobity.biochem.duke.edu, save the summary,
    then use parse_molprobity() on the saved file, or compare_validation()
    on the pre- and post-refinement pair.
    """
    raise NotImplementedError(
        "MolProbity has no stable programmatic interface and is not wrapped. "
        "Run it at molprobity.biochem.duke.edu, save the summary, then call "
        "parse_molprobity(saved_file) or compare_validation(before, after)."
    )
