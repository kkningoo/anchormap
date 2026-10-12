"""Which external tools are present, and what to do when one is not.

Three of the tools this package drives cannot be installed with pip. Rather
than letting a missing one surface as an ImportError from somewhere deep in a
wrapper, every wrapper calls `require()` first and fails with the tool's name
and where to get it.

Run `python -m anchormap.requirements` to see what is available.
"""

from __future__ import annotations

import importlib.util
import shutil
from dataclasses import dataclass


class MissingToolError(RuntimeError):
    """An external tool this function drives is not installed."""


@dataclass(frozen=True)
class Tool:
    key: str
    name: str
    used_for: str
    kind: str  # "python" (import) or "binary" (on PATH)
    probe: str
    install: str
    pip_installable: bool


TOOLS = {
    "numpy": Tool(
        "numpy", "NumPy", "arithmetic throughout",
        "python", "numpy", "pip install numpy", True,
    ),
    "biopython": Tool(
        "biopython", "Biopython", "structure parsing, Shrake-Rupley accessibility",
        "python", "Bio", "pip install biopython", True,
    ),
    "freesasa": Tool(
        "freesasa", "FreeSASA", "buried surface area at subunit interfaces",
        "python", "freesasa", "pip install freesasa", True,
    ),
    "openmm": Tool(
        "openmm", "OpenMM", "restrained energy minimisation (prep.minimise)",
        "python", "openmm",
        "conda install -c conda-forge openmm  (recommended)  |  pip install openmm",
        True,
    ),
    "pdbfixer": Tool(
        "pdbfixer", "PDBFixer", "missing atoms and hydrogens before minimisation",
        "python", "pdbfixer",
        "conda install -c conda-forge pdbfixer  (recommended)  |  pip install pdbfixer",
        True,
    ),
    "pyrosetta": Tool(
        "pyrosetta", "PyRosetta", "interface energetics (interfaces.analyse)",
        "python", "pyrosetta",
        "Not on PyPI. Free academic licence from RosettaCommons at "
        "rosettacommons.org/software/licensing, then follow their install "
        "instructions. iGEM teams qualify for the academic licence.",
        False,
    ),
    "molprobity": Tool(
        "molprobity", "MolProbity", "all-atom validation (prep.validate)",
        "binary", "molprobity.molprobity",
        "Not on PyPI. Either use the web server at molprobity.biochem.duke.edu "
        "and record the output manually, or build locally from "
        "github.com/rlabduke/MolProbity.",
        False,
    ),
    "reduce": Tool(
        "reduce", "Reduce", "hydrogen addition with Asn/Gln/His flip correction",
        "binary", "reduce",
        "Not on PyPI. Ships with MolProbity, or build from "
        "github.com/rlabduke/reduce.",
        False,
    ),
}


def available(key: str) -> bool:
    tool = TOOLS[key]
    if tool.kind == "python":
        return importlib.util.find_spec(tool.probe) is not None
    return shutil.which(tool.probe) is not None


def require(*keys: str) -> None:
    """Called at the top of every wrapper. Raises naming the tool and its source."""
    missing = [TOOLS[k] for k in keys if not available(k)]
    if not missing:
        return
    lines = [f"{t.name} is required for this step ({t.used_for}).\n  Install: {t.install}"
             for t in missing]
    raise MissingToolError("\n\n".join(lines))


def report() -> dict[str, bool]:
    return {k: available(k) for k in TOOLS}


def main() -> None:
    status = report()
    width = max(len(t.name) for t in TOOLS.values())
    print("anchormap external tool check\n")
    for key, ok in status.items():
        t = TOOLS[key]
        mark = "yes" if ok else "NO "
        pip = "" if t.pip_installable else "   [not pip-installable]"
        print(f"  {mark}  {t.name:<{width}}  {t.used_for}{pip}")
    missing = [TOOLS[k] for k, ok in status.items() if not ok]
    if missing:
        print("\nMissing:\n")
        for t in missing:
            print(f"  {t.name}\n    {t.install}\n")
    else:
        print("\nAll tools present.")


# Also installed as the console command `anchormap-check`.
if __name__ == "__main__":
    main()
