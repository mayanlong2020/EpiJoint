import os
from utils.logger import console
from scripts.phase4_qa import run_phase4

runs_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "runs")
latest_run = sorted([os.path.join(runs_dir, d) for d in os.listdir(runs_dir) if d.startswith("run_")])[-1]
run_phase4(latest_run)
