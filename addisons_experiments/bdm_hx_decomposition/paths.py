"""Output locations for this experiment.

FIGURES holds plots and metrics that were explicitly asked for; it is tracked in
git and is the only place a reader should need to look.  AGENT_ARTIFACTS holds
everything an agent generated for its own verification -- resolution sweeps,
convergence data, raw DOF vectors -- and is gitignored.

Default new output to AGENT_ARTIFACTS.  Write to FIGURES only for a deliverable
that was requested by name.
"""

from pathlib import Path

_ROOT = Path(__file__).parent

FIGURES = _ROOT / "figures"
AGENT_ARTIFACTS = _ROOT / "agent_artifacts"

FIGURES.mkdir(parents=True, exist_ok=True)
AGENT_ARTIFACTS.mkdir(parents=True, exist_ok=True)
