# anchormap

Surface, interface and support-compatibility analysis for immobilising multimeric enzymes.

Attaching an oligomeric enzyme to a solid support raises five questions that have to be answered in order. `anchormap` answers them from coordinates.

| Module | Question |
| --- | --- |
| `prep` | Are the coordinates fit to compute on? |
| `accessibility` | Which surface sites exist on the assembly, not the subunit? |
| `coupling` | Which are usable, and can they be crosslinked in pairs? |
| `interfaces` | Which subunit interface opens first? |
| `liabilities` | What degrades first, and does assembly protect it? |
| `electrostatics` | What is usable in a PROPKA run, and what is not? |
| `support` | Will the enzymes load, and will substrate reach them? |
| `constants` | Every shared threshold and maximum, with its provenance |

## Install

Install directly from the repository:

```bash
pip install git+https://github.com/kkningoo/anchormap.git
```

That pulls in numpy, biopython and freesasa. To work on the code instead, clone it and install in editable mode:

```bash
git clone https://github.com/kkningoo/anchormap.git
cd anchormap
pip install -e .            # anchormap, numpy, biopython, freesasa
pip install -e ".[refine]"  # adds OpenMM and PDBFixer
pip install -e ".[dev]"     # adds pytest
```

The Python package and three of its dependencies come from pip. **Three tools do not**, and must be installed separately.

For OpenMM and PDBFixer, conda-forge is the better-supported route:

```bash
conda install -c conda-forge openmm pdbfixer
```

### External tools you must obtain yourself

| Tool | Needed for | How to get it |
| --- | --- | --- |
| PyRosetta | `interfaces.analyse` | Not on PyPI. Free academic licence from [RosettaCommons](https://rosettacommons.org/software/licensing), then their installer. iGEM teams qualify. |
| MolProbity | `prep.parse_molprobity` | Not on PyPI and has no stable programmatic interface. Run it at the [web server](http://molprobity.biochem.duke.edu), save the summary, then parse the saved file with `prep.parse_molprobity()` or `prep.compare_validation()`. |
| Reduce | hydrogen addition with flips, inside `prep.minimise` | Not on PyPI. Ships with MolProbity, or build from [rlabduke/reduce](https://github.com/rlabduke/reduce). |

Everything that does not need these still works. The whole of `support`, the ranking and convergence logic in `interfaces`, the reachability and coverage functions in `coupling`, the MolProbity summary parser in `prep`, and all of `guards` run on a bare install.

### Check what you have

```bash
anchormap-check
```

Prints each tool, what it is used for, whether it is present, and the install route for anything missing. Every wrapper calls the same check before it runs, so a missing tool surfaces as `MissingToolError` naming the tool and its source rather than an `ImportError` from inside a dependency.

## Use

```python
from anchormap import prep, accessibility, coupling, interfaces, support

s = prep.load_assembly("4LNK.pdb1", expected_chains=12,
                       expected_metals=24, source_first_residue=2)
s = prep.minimise(s)

rsa = accessibility.compare_states(monomer, s, resname_filter="LYS")
summary = accessibility.summarise(rsa)          # how many sites assembly removes

excluded = coupling.active_site_residues(s, ligand_resnames={"GLU", "ADP"})
pairs = coupling.reachable_pairs(sites, crosslinker="genipin")
coverage = coupling.multipoint_coverage(sites, pairs)

results = [interfaces.analyse(s, p) for p in subunit_pairs]
weakest = interfaces.rank(results)              # weakest first, by energy density
interfaces.flag_nonconverged(results)           # reported, never silently ranked
```

The support module needs no structure and runs standalone:

```python
gs   = support.AdsorptionPrediction("GS",   pI=5.40, operating_pH=7.0, support_charge="cationic")
ppk2 = support.AdsorptionPrediction("PPK2", pI=9.66, operating_pH=7.0, support_charge="cationic")

support.differential_adsorption([gs, ppk2])
# {'favoured': ['GS'], 'disfavoured': ['PPK2'], 'ratio_biased': True, ...}

support.minimum_pore_size((146.5, 140.9, 114.1))        # 117.1 A
support.polyelectrolyte_barrier_risk("polyanionic", "cationic")
```

## Why the guards raise

`guards.py` encodes six preconditions. Each corresponds to an error made during development that produced a plausible wrong answer and survived review:

| Guard | What it caught |
| --- | --- |
| `require_minimised` | Accessibility on unrefined coordinates gave three false lysine exclusions |
| `check_numbering` | PDBFixer renumbering from 1 caused three mutation-position errors |
| `check_selection_size` | Selecting a ligand by residue name returned 27,550 atoms |
| `warn_crystallisation_additives` | Sulfates added two residues to an active-site definition |
| `check_metals_present` | Stripping metals doubled an electrostatic penalty and gave pKa values from −5.6 to +22 |
| `check_biological_assembly` | Loading an asymmetric unit instead of a biological assembly halved a particle |
| `verify_mutation` | A PyMOL mutagenesis applied to one chain of twelve, caught only after nine hours of minimisation |
| `check_chain_symmetry` | Symmetry-related chains disagreeing, which means a merge or renaming error |

They raise rather than warn, and each can be overridden — but the override has to be written in the calling code, where a reviewer can see it.

Four interfaces are deliberately narrow for the same reason.

`accessibility.compare_states` requires both oligomeric states and has no single-state path, because a single-state result is what produced the wrong coupling map.

`coupling.reachable_pairs` has no default span and raises if you give neither a named crosslinker nor an explicit window. A 10 Å floor rather than 5 Å excluded four of six real cross-ring contacts here, and the counts are meaningless without the window beside them.

`interfaces.analyse` rejects `relax_cycles=0`, because unrelaxed crystal poses return positive ΔG dominated by clash repulsion.

`interfaces.rank` refuses to rank fixed-backbone and backbone-free runs together, because backbone-free relaxation gives systematically more favourable ΔG and the two are not comparable.

## Scope

Written for one project and exercised on two enzymes: a twelve-subunit glutamine synthetase built from a crystal structure, and a four-subunit polyphosphate kinase from a prediction. It has not been run on anything else.

Implemented and exercised: solvent accessibility in both oligomeric states, per-chain variation and assembly symmetry checking, terminal accessibility, coupling-site enumeration with active-site and additive exclusion, reachability and clustering against a stated crosslinker span, span sweeps stratified by structural group, inter-subunit contact analysis, buried surface by FreeSASA, interface energetics through Rosetta with convergence flagging, chemical liability scanning, and PROPKA parsing.

`prep.validate` still raises `NotImplementedError`: MolProbity has no stable programmatic interface here and the project drove it by hand.

It is not a packaged tool with a stable API, and it does not bundle or install the structural-biology software it drives. It is a readable pipeline with its failure modes written down in a form that stops them recurring.

## Tests

```bash
pytest
```

Twenty-three smoke tests covering the guards, the shared constants and the pure-computation helpers. Not a suite; external-tool wrappers are not covered.

## Licence

MIT. See `LICENSE`.

---

Developed by iGEM VIT Vellore 2026 (Team Nishvasa) at the Vellore Institute of Technology, for the Nishvasa project. MIT licensed.
