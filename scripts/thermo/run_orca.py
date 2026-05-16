# ~/orca_6_0_1/orca co2_280.inp --use-hwthread-cpus > ./output/co2_280.out 
# # scp -P 2019 co2_280.out pedro@146.164.73.120:~/Desktop/Artigos/ayahuasca/compounds/acetaldehyde/quantum/co2_280.out

from datetime import datetime
from pathlib import Path
import time
from contextlib import contextmanager

from fabric import Connection, Config
from paramiko import SSHException
from paramiko.ssh_exception import NoValidConnectionsError
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.traceback import install

from data import COMPOUNDS_PATH, compounds

install()
console = Console()

# =========================
# SSH / ORCA CONFIG
# =========================

REMOTE_ORCA_PATH = "~/orca_6_0_1/orca"

SSH_HOST = "146.164.73.124"
SSH_USER = "pedro"
SSH_PORT = 22
SSH_PASSWORD = "bi0inf0"

MAX_RETRIES = 3
RETRY_DELAY = 5


# =========================
# SSH CONNECTION
# =========================

@contextmanager
def ssh_connection():
    config = Config(overrides={"connect_kwargs": {"password": SSH_PASSWORD}})
    conn = Connection(
        host=SSH_HOST,
        user=SSH_USER,
        port=SSH_PORT,
        config=config,
    )
    try:
        yield conn
    finally:
        conn.close()


# =========================
# ORCA JOB
# =========================

def run_orca_job(conn: Connection, cpd: str, cpd_id: str):
    """Run a single ORCA job using an existing SSH connection."""

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    input_file = f"{cpd}_{cpd_id}.inp"

    # ---- Local paths
    local_quantum_dir = Path(COMPOUNDS_PATH) / cpd / "quantum"
    # local_quantum_dir.mkdir(parents=True, exist_ok=True)

    local_input = local_quantum_dir / input_file
    local_output = local_quantum_dir / f"{cpd}_{cpd_id}.out"

    # ---- Remote paths
    remote_base = f"/home/pedro/desktop/projects/ayahuasca/analysis/quantum_chemistry/{cpd}"
    remote_job_dir = f"{remote_base}/job_{timestamp}"
    remote_input = f"{remote_job_dir}/{input_file}"
    remote_output = f"{remote_job_dir}/{cpd}_{cpd_id}.out"

    console.log(f"[bold blue]Starting ORCA job for {cpd} ({cpd_id})[/bold blue]")

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            # Prepare remote job dir
            conn.run(f"mkdir -p {remote_job_dir}")

            # Upload input
            conn.put(str(local_input), remote=remote_job_dir)

            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}")
            ) as progress:
                task = progress.add_task("[green]Running ORCA...", start=False)
                progress.start_task(task)

                conn.run(
                    f"bash -l -c '{REMOTE_ORCA_PATH} {remote_input} --use-hwthread-cpus > {remote_output}'",
                    hide=True, pty=True
                )

                progress.update(task, description="[green]ORCA finished")

            # Retrieve output
            conn.get(remote_output, str(local_output))

            console.log(
                f"[bold green]Finished {cpd}. Output saved to {local_output}[/bold green]"
            )
            return

        except (SSHException, NoValidConnectionsError) as e:
            console.log(f"[bold red]Attempt {attempt} failed: {e}[/bold red]")
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_DELAY)
            else:
                raise RuntimeError(f"ORCA job failed for {cpd}") from e


# =========================
# MAIN LOOP
# =========================

if __name__ == "__main__":
    with ssh_connection() as conn:
        for cpd, info in compounds.items():
            # Check if there already is a g_minus_eel_kcalpmol value on an output file
            output_path = Path(COMPOUNDS_PATH) / cpd / "thermo" / f"{cpd}_{info['pubchem_id']}.json"
            if output_path.exists():
                with open(output_path, "r") as f:
                    data = f.read()
                    if "g_minus_eel_kcalpmol" in data:
                        console.log(f"[yellow]Skipping {cpd}, already has g_minus_eel_kcalpmol[/yellow]")
                        continue
            run_orca_job(conn, cpd, info["pubchem_id"])
