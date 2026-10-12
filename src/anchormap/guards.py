"""Preconditions that caught real errors in this project.

Every guard here corresponds to a mistake we made. They raise by default and
can be bypassed with allow=True, which is deliberate: the override has to be
written down in the calling code, where a reader can see it.
"""

from __future__ import annotations

import warnings


class PreconditionError(RuntimeError):
    """A structure or selection failed a check that previously caught an error."""


def mark_minimised(path: str) -> None:
    """Record that a file came out of a refinement step."""
    import os
    open(os.fspath(path) + '.minimised', 'w').close()


def require_minimised(structure, allow: bool = False) -> None:
    """Accessibility and pKa calculations on unrefined coordinates are unreliable.

    An unminimised predicted tetramer in this project carried 0.6 A and 2.5 A
    contacts and produced a lysine exclusion list (Lys8, Lys12, Lys73) that did
    not survive rebuilding. None of those three is actually buried on assembly.
    """
    if allow:
        return
    if getattr(structure, "minimised", False):
        return
    import os
    if isinstance(structure, (str, os.PathLike)):
        if os.path.exists(os.fspath(structure) + '.minimised'):
            return
    raise PreconditionError(
        "Structure is not marked as minimised. Unrefined coordinates produce "
        "false burial calls. Run prep.minimise() first, or pass allow=True and "
        "say why in a comment."
    )


def check_numbering(structure, source_first_residue: int) -> None:
    """PDBFixer silently renumbers from 1.

    4LNK starts at residue 2, so any structure passed through PDBFixer is
    offset by -1 against the literature numbering. Three mutation-position
    errors in this project traced to this.
    """
    first = structure.first_residue_number()
    if first != source_first_residue:
        raise PreconditionError(
            f"Residue numbering starts at {first}, source starts at "
            f"{source_first_residue}. Offset {first - source_first_residue}. "
            "Renumber before mapping any literature residue position."
        )


def check_selection_size(n_atoms: int, expected_max: int, label: str) -> None:
    """A ligand selection by residue name matches every copy in the protein.

    `resn GLU` returned 27,550 atoms instead of an active site here. An active
    site is tens of atoms, not thousands.
    """
    if n_atoms > expected_max:
        raise PreconditionError(
            f"{label}: selection returned {n_atoms} atoms, expected at most "
            f"{expected_max}. Selecting by residue name alone matches every "
            "copy in the structure; constrain by chain and residue number."
        )


def warn_crystallisation_additives(resnames: set[str]) -> None:
    """Additives contaminate active-site definitions.

    Sulfates pulled Lys86 and Lys88 into the GS active-site selection in 2 of
    12 chains. Excluding additives removed them.
    """
    additives = {"SO4", "PO4", "GOL", "EDO", "PEG", "MPD", "ACT", "CL", "NA", "K"}
    found = resnames & additives
    if found:
        warnings.warn(
            f"Crystallisation additives present in selection: {sorted(found)}. "
            "Exclude them before defining a binding site.",
            stacklevel=2,
        )


def check_metals_present(n_metals: int, expected: int) -> None:
    """Protonation prediction is unreliable at metal sites in a stripped structure.

    Stripping the 24 Mg from GS roughly doubled the apparent high-pH penalty
    and returned pKa values from -5.6 to +22.
    """
    if n_metals < expected:
        raise PreconditionError(
            f"{n_metals} metal ions present, {expected} expected. Run with "
            "metals in place; pKa and electrostatic output is not interpretable "
            "without them."
        )


def check_biological_assembly(n_chains: int, expected: int) -> None:
    """The asymmetric unit is not the biological assembly.

    4LNK's ASU is one hexameric ring; the dodecamer requires the .pdb1 file
    and chain renaming on merge.
    """
    if n_chains != expected:
        raise PreconditionError(
            f"{n_chains} chains present, biological assembly has {expected}. "
            "Load the .pdb1 assembly file, not the asymmetric unit."
        )


def verify_mutation(pdb: str, resnum: int, expected_resname: str,
                    expected_count: int) -> None:
    """Every copy of a mutated residue carries the substitution.

    A PyMOL mutagenesis wizard silently applied a substitution to one chain of
    a twelve-chain assembly here. Nine hours of minimisation later, the
    structure was still mostly wild type. One grep would have caught it.

        grep " CA " f.pdb | awk 'substr($0,23,4)+0==195 {print substr($0,18,3)}'
    """
    found = []
    for line in open(pdb):
        if (line.startswith('ATOM') and line[12:16].strip() == 'CA'
                and line[22:26].strip() == str(resnum)):
            found.append((line[21], line[17:20].strip()))

    wrong = [(ch, rn) for ch, rn in found if rn != expected_resname]
    if len(found) != expected_count or wrong:
        raise PreconditionError(
            f"Position {resnum}: expected {expected_resname} in "
            f"{expected_count} chains, found {len(found)} copies"
            + (f", wrong in {wrong}" if wrong else "")
            + ". A mutagenesis step applied to one chain and not the rest is "
              "the usual cause. Check before spending compute on it."
        )


def check_chain_symmetry(per_chain: dict[str, float], tol: float = 0.005) -> dict:
    """Symmetry-related chains should agree. Disagreement means a build error.

    A correctly merged stacked-ring assembly gives identical values for each
    chain and its symmetry mate. Twelve chains collapsing into six matched
    pairs to three decimal places is direct evidence the biological assembly
    was built correctly; a merge or renaming error breaks the pattern.
    """
    vals = sorted(per_chain.items())
    n = len(vals)
    if n % 2:
        return {'symmetric': None, 'reason': 'odd chain count'}
    half = n // 2
    pairs, ok = [], True
    for i in range(half):
        a, b = vals[i], vals[i + half]
        match = abs(a[1] - b[1]) <= tol
        ok &= match
        pairs.append({'chains': (a[0], b[0]), 'values': (a[1], b[1]),
                      'match': match})
    return {'symmetric': ok, 'pairs': pairs,
            'reason': None if ok else 'symmetry-related chains disagree; '
                                      'check the merge and chain renaming'}


def require_metals(path: str, expected: int, context: str = "") -> None:
    """Fail loudly when a file that should carry metals does not.

    Minimisation usually strips metals, because the force field has no
    parameters for them, and the stripped file is often the one left in the
    working directory under the obvious name. Coordination distances, metal
    site validation and active site definitions computed on that file return
    nothing -- or worse, return a plausible number from a site that no longer
    has a metal in it.

    This cost real time here: a magnesium site superposition ran against a
    minimised structure containing no magnesium, and the failure surfaced only
    as an empty PyMOL selection.

    Call before any analysis that depends on a metal being present, and point
    it at the file that kept its metals rather than the minimised one.
    """
    known = {'MG', 'MN', 'ZN', 'CA', 'FE', 'NI', 'CO', 'CU'}
    found = sum(1 for line in open(path)
                if line.startswith('HETATM') and line[17:20].strip() in known)
    if found != expected:
        raise PreconditionError(
            f"{path} contains {found} metal atoms, expected {expected}"
            + (f" ({context})" if context else "")
            + ". Minimisation strips metals because the force field has no "
              "parameters for them; the file that kept them is usually the "
              "pre-minimisation structure, or one where they were reattached "
              "afterwards. Coordination geometry from a stripped file is "
              "meaningless, and an empty selection is the mildest way it can "
              "fail."
        )


def which_file_has_metals(*paths: str) -> dict:
    """Metal count per file, to find the right input without guessing.

    Run this once on every candidate before starting a coordination analysis.
    """
    known = {'MG', 'MN', 'ZN', 'CA', 'FE', 'NI', 'CO', 'CU'}
    out = {}
    for p in paths:
        try:
            out[p] = sum(1 for line in open(p)
                         if line.startswith('HETATM')
                         and line[17:20].strip() in known)
        except OSError as exc:
            out[p] = f"unreadable: {exc}"
    return out
