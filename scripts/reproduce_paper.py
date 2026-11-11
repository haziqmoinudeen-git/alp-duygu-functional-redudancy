"""Run or resume all paired trajectories, validate, and regenerate figures."""

import subprocess
import sys

for args in [
    ["-m", "enzyme_specialisation", "simulate", "--config", "configs/paper.yaml"],
    ["scripts/validate_outputs.py"],
    ["-m", "enzyme_specialisation", "plot"],
]:
    subprocess.run([sys.executable, *args], check=True)
