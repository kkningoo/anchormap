"""Smoke tests for the parts of anchormap that are fully implemented.

This is not a test suite. It checks that the guards fire, that the
pure-computation helpers return the values the package was developed against,
and that the package imports cleanly. Functions wrapping external tools are
not covered.
"""

import pytest

from anchormap import coupling, guards, interfaces, support


class _Unminimised:
    minimised = False

    def first_residue_number(self):
        return 1


def test_minimisation_guard_raises():
    with pytest.raises(guards.PreconditionError):
        guards.require_minimised(_Unminimised())


def test_minimisation_guard_is_overridable():
    guards.require_minimised(_Unminimised(), allow=True)


def test_numbering_guard_catches_offset():
    with pytest.raises(guards.PreconditionError):
        guards.check_numbering(_Unminimised(), source_first_residue=2)


def test_selection_size_guard():
    with pytest.raises(guards.PreconditionError):
        guards.check_selection_size(27550, expected_max=500, label="ligand")


def test_assembly_guard_catches_asymmetric_unit():
    with pytest.raises(guards.PreconditionError):
        guards.check_biological_assembly(n_chains=6, expected=12)


def test_metals_guard():
    with pytest.raises(guards.PreconditionError):
        guards.check_metals_present(n_metals=0, expected=24)


def test_differential_adsorption_direction():
    """Opposite net charge against a cationic support biases the loaded ratio."""
    acidic = support.AdsorptionPrediction("acidic", 5.40, 7.0, "cationic")
    basic = support.AdsorptionPrediction("basic", 9.66, 7.0, "cationic")
    result = support.differential_adsorption([acidic, basic])
    assert result["favoured"] == ["acidic"]
    assert result["disfavoured"] == ["basic"]
    assert result["ratio_biased"] is True


def test_no_bias_when_both_enzymes_share_sign():
    a = support.AdsorptionPrediction("a", 5.0, 7.0, "cationic")
    b = support.AdsorptionPrediction("b", 5.5, 7.0, "cationic")
    assert support.differential_adsorption([a, b])["ratio_biased"] is False


def test_pore_size_uses_smallest_extent():
    """A cylinder enters a pore end-on, so the limiting dimension is smallest."""
    assert support.minimum_pore_size((146.5, 140.9, 114.1)) == pytest.approx(117.1)


def test_barrier_flagged_only_for_opposite_polyelectrolyte():
    assert support.polyelectrolyte_barrier_risk("polyanionic", "cationic") is not None
    assert support.polyelectrolyte_barrier_risk("polyanionic", "anionic") is None


def test_interface_ranking_is_weakest_first():
    def mk(pair, area, dG):
        return interfaces.InterfaceResult(pair, area, dG, None, None, 550.0, -10.0, True)

    results = [
        mk(("A", "H"), 4250.8, -90.13),
        mk(("A", "F"), 3459.8, -44.44),
        mk(("A", "B"), 3477.3, -35.75),
    ]
    assert [r.pair for r in interfaces.rank(results)] == [("A", "B"), ("A", "F"), ("A", "H")]


def test_nonconverged_runs_are_excluded_from_ranking():
    def mk(pair, dG, converged):
        return interfaces.InterfaceResult(pair=pair, buried_area=3000.0, dG=dG,
                                          converged=converged)

    results = [mk(("A", "B"), -30.0, True), mk(("C", "D"), -25.0, False)]
    assert [r.pair for r in interfaces.rank(results)] == [("A", "B")]


def test_relax_cycles_cannot_be_zero():
    """Without PyRosetta the tool check fires first; either refusal is correct."""
    from anchormap.requirements import MissingToolError
    with pytest.raises((ValueError, MissingToolError)):
        interfaces.analyse("nonexistent.pdb", ("A", "B"), relax_cycles=0)


def test_crosslinker_spans_are_ranges_not_lengths():
    lo, hi = coupling.CROSSLINKER_SPAN["genipin"]
    assert lo < hi


# --- added in 0.2.0, covering the ported implementations -------------------

from anchormap import constants, electrostatics, liabilities


def test_tien_maxima_are_theoretical_not_empirical():
    """MET/SER/THR are the three that differ between the two Tien tables."""
    t = constants.TIEN_THEORETICAL
    assert (t['MET'], t['SER'], t['THR']) == (224, 155, 172)


def test_crosslinker_floor_is_not_ten():
    """A 10 A floor excluded four of six real cross-ring contacts."""
    lo, _ = constants.CROSSLINKER_SPAN['genipin']
    assert lo == 5.0


def test_reachable_pairs_requires_an_explicit_window():
    with pytest.raises(ValueError):
        coupling.reachable_pairs({})
    with pytest.raises(ValueError):
        coupling.reachable_pairs({}, crosslinker='genipin', span=(5.0, 16.0))


def test_rank_refuses_to_mix_protocols():
    def mk(pair, bb):
        return interfaces.InterfaceResult(pair, 3000.0, -30.0, backbone_free=bb)
    with pytest.raises(ValueError):
        interfaces.rank([mk(('A', 'B'), False), mk(('C', 'D'), True)])


def test_mutation_verification_catches_partial_application():
    import tempfile, os
    lines = []
    for i, ch in enumerate('ABCDEFGHIJKL'):
        rn = 'TYR' if ch != 'F' else 'HIS'          # one chain not mutated
        lines.append(f"ATOM  {i:5d}  CA  {rn} {ch} 195      0.000   0.000   0.000\n")
    fd, path = tempfile.mkstemp(suffix='.pdb')
    with os.fdopen(fd, 'w') as f:
        f.writelines(lines)
    try:
        with pytest.raises(guards.PreconditionError):
            guards.verify_mutation(path, 195, 'TYR', 12)
    finally:
        os.unlink(path)


def test_chain_symmetry_detects_a_broken_merge():
    good = {'A': 0.277, 'B': 0.271, 'C': 0.310,
            'G': 0.277, 'H': 0.271, 'I': 0.310}
    assert guards.check_chain_symmetry(good)['symmetric'] is True
    bad = dict(good, I=0.190)
    assert guards.check_chain_symmetry(bad)['symmetric'] is False


def test_unphysical_pka_values_are_flagged():
    r = {'titratable': [
        electrostatics.Titratable('GLU', '196', 'A', 22.0),
        electrostatics.Titratable('GLU', '132', 'A', -5.62),
        electrostatics.Titratable('LYS', '43', 'A', 7.05),
    ]}
    bad = electrostatics.unphysical(r)
    assert len(bad) == 2
    assert all(not t.physical for t in bad)


def test_depressed_lysine_is_reported_as_shifted():
    r = {'titratable': [electrostatics.Titratable('LYS', '43', 'A', 7.05)]}
    s = electrostatics.shifted(r)
    assert len(s) == 1 and s[0].shift == pytest.approx(-3.45)


def test_liability_motifs_recognised():
    assert 'G' in constants.DEAMIDATION_FOLLOWERS
    assert 'G' in constants.ISOMERISATION_FOLLOWERS
    assert set('MCW') == set(constants.OXIDISABLE)


def test_metal_guard_catches_a_stripped_file():
    """The failure that cost an afternoon: minimised file, no metals."""
    import tempfile, os
    fd, path = tempfile.mkstemp(suffix='.pdb')
    with os.fdopen(fd, 'w') as f:
        f.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")
    try:
        with pytest.raises(guards.PreconditionError) as e:
            guards.require_metals(path, 4, 'PPK2 tetramer')
        assert 'strips metals' in str(e.value)
        assert guards.which_file_has_metals(path) == {path: 0}
    finally:
        os.unlink(path)
