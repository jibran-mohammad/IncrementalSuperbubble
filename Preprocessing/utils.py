import time
import subprocess

def run_command(cmd, label=None):
    """
    Run a shell command with:
      - live stdout/stderr
      - timing
      - consistent logging
      - proper failure detection
    """

    tag = f" [{label}]" if label else ""

    print(f"Running command{tag}:\n{cmd}\n")

    start_time = time.time()
    proc = subprocess.run(cmd, shell=True)
    elapsed = time.time() - start_time

    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed{tag} (exit code {proc.returncode})"
        )

    print(f"--- Time : {elapsed:.2f} seconds ---{tag}")
