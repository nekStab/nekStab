import os, sys, glob, subprocess, shutil

try:
    sys.path.insert(0, os.path.join(os.environ["NEKSTAB_SOURCE_ROOT"], "bin"))
except KeyError:
    raise EnvironmentError("NEKSTAB_SOURCE_ROOT environment variable is not set.")
from nekStab_tools import *

# ---- SIMULATION TARGET SETTINGS ----
target_end_time = 200.0  # Set your desired end time here
single_job_time = 50.0  # Set your desired field write interval here

pbs_script = "jz.pbs"
local_sh_script = "run_nohup.sh"

# ---- CASE-SPECIFIC SETTINGS ----
case_name = get_case_name_from_usr()
case_initial = case_name[0]
history_file = f"{case_name}.his"
parameter_glob = "*.par"


# NOTE: For test runs only. Not used for production cases.
def prepare_par_for_restart(parameter_glob, single_job_time, case_preserving_adjust_par_file):
    """
    Update the .par file to set startFrom = rst_1cyl0.f00001 and endTime = single_job_time.
    This is intended for starting from scratch in test runs only.
    """
    import glob

    par_files = glob.glob(parameter_glob)
    if not par_files:
        raise FileNotFoundError(f"No .par file matching {parameter_glob} found.")
    par_file = par_files[0]
    case_preserving_adjust_par_file(par_file, {("GENERAL", "startFrom"): "rst_1cyl0.f00001", ("GENERAL", "endTime"): str(single_job_time)})


log_file = "logfile"
log_error_file = "logerror"
log_restart_file = "logrestart"

log_files = [log_file, log_error_file, log_restart_file]

lift_drag_file = "lift_drag.dat"
total_energy_file = "total_energy.dat"
total_enstrophy_file = "total_enstrophy.dat"

consolidate_file_patterns = [
    ("lift_drag.dat_*", "lift_drag.all"),
    ("total_energy.dat_*", "total_energy.all"),
    ("total_enstrophy.dat_*", "total_enstrophy.all"),
    (f"{case_name}.his_*_*", f"{case_name}.all"),
]
default_job_name = lambda cinit, re, t: f"{cinit}{re}_{t}"[:8]

# ---- MAIN WORKFLOW ----
if __name__ == "__main__":
    parameter_files = glob.glob(parameter_glob)
    if not parameter_files:
        print(f"No parameter file found in current directory matching {parameter_glob}.")
        sys.exit(1)
    parameter_file = parameter_files[0]

    if "--consolidate" in sys.argv:
        consolidate_files(consolidate_file_patterns)
        sys.exit(0)

    working_directory = os.getcwd()
    case_reynolds = os.path.basename(working_directory)

    print("\n==== [SIMULATION INFO] ====")
    print(f"[INFO] Target final time: {target_end_time}")
    print(f"[INFO] Case name: {case_name}")
    print(f"[INFO] Parameter file: {parameter_file}")
    print("==========================\n")

    final_time_in_his_str, final_time_in_his = get_time_str_and_float(history_file, which_time="final")
    initial_time_in_his_str, initial_time_in_his = get_time_str_and_float(history_file, which_time="initial")
    latest_restart_file = find_latest_final_dns_file(case_name, verbose=True)
    velocity_file_time = extract_time_from_binary(latest_restart_file)
    if np.isnan(velocity_file_time):
        print(f"The file {latest_restart_file} does not exist.")
        sys.exit()
    print("---- [TIME CHECK] ----")
    print(f"[INFO] Initial time in {history_file}: {initial_time_in_his}")
    print(f"[INFO] Final time in {history_file}: {final_time_in_his}")
    print(f"[INFO] Time in velocity file: {velocity_file_time}")
    print("----------------------\n")

    # Generate the expected job name for this case
    if velocity_file_time < target_end_time:
        next_time = velocity_file_time + single_job_time
        if velocity_file_time < single_job_time:
            next_time = single_job_time
        expected_job_name = default_job_name(case_initial, case_reynolds, next_time)
    else:
        expected_job_name = None

    # Check if any job is running (PBS/SLURM or local)
    job_running = False
    
    # Check for PBS/SLURM jobs first (check for exact job name)
    if expected_job_name:
        try:
            if shutil.which("qstat") is not None:  # PBS
                result = subprocess.run(['qstat', '-u', os.environ.get('USER', 'unknown')], capture_output=True, text=True)
                if result.returncode == 0:
                    for line in result.stdout.strip().split('\n')[2:]:  # Skip header lines
                        if line.strip() and expected_job_name in line:
                            print(f"[INFO] Found specific PBS job '{expected_job_name}' running. Skipping all checks.. only consolidating.")
                            job_running = True
                            break
            elif shutil.which("squeue") is not None:  # SLURM
                result = subprocess.run(['squeue', '-u', os.environ.get('USER', 'unknown')], capture_output=True, text=True)
                if result.returncode == 0:
                    for line in result.stdout.strip().split('\n')[1:]:  # Skip header line
                        if line.strip() and expected_job_name in line:
                            print(f"[INFO] Found specific SLURM job '{expected_job_name}' running. Skipping all checks.. only consolidating.")
                            job_running = True
                            break
        except Exception:
            pass
    
    # Check for local nek5000 processes (simpler check for local case)
    if not job_running:
        try:
            result = subprocess.run(['pgrep', '-f', 'nek5000'], capture_output=True, text=True)
            if result.returncode == 0:
                pids = result.stdout.strip().split('\n')
                pid_count = len(pids)
                print(f"[INFO] Found {pid_count} local nek5000 job{'s' if pid_count > 1 else ''} running. Skipping all checks.. only consolidating.")
                job_running = True
        except Exception:
            pass
    
    if velocity_file_time < target_end_time:
        next_time = velocity_file_time + single_job_time
        if velocity_file_time < single_job_time:
            next_time = single_job_time

        if not job_running:
            print("==== [PARAMETER UPDATE & JOB SUBMISSION] ====")
            print(f"[INFO] Updating parameter file '{parameter_file}' for next run...")
            case_preserving_adjust_par_file(
                parameter_file,
                {
                    ("GENERAL", "startFrom"): str(latest_restart_file),
                    ("GENERAL", "endTime"): str(next_time),
                    ("VELOCITY", "viscosity"): f"-{float(case_reynolds)}",
                },
            )
            job_name = default_job_name(case_initial, case_reynolds, next_time)
            print(f"[INFO] Submitting job '{job_name}' using script '{local_sh_script}'...")
            submit_job(job_name, local_sh_script, working_directory, pbs_script)
            print("============================================\n")
        else:
            print("==== [JOB STATUS] ====")
            print("[INFO] Job is currently running. Skipping parameter update and job submission.")
            print("======================\n")
    else:
        print("==== [JOB STATUS] ====")
        print(f"[INFO] Current time: {final_time_in_his} >= Target time: {target_end_time}")
        print(f"[INFO] Simulation has reached target time. No further action needed.")
        print("======================\n")

    # Only do backups when NO job is running
    if not job_running and initial_time_in_his >= 0.0:
        if initial_time_in_his == 0 and final_time_in_his == 0:
            backup_suffix = f"_{initial_time_in_his_str}"
        else:
            backup_suffix = f"_{initial_time_in_his_str}_{final_time_in_his_str}"

        print(f"==== [BACKUP & CLEANUP] ====")
        print(f"[INFO] Backup suffix: {backup_suffix}")

        # Check if backup files already exist and skip backing up if they do
        files_to_backup = [latest_restart_file, history_file, lift_drag_file, total_energy_file, total_enstrophy_file] + log_files
        files_needing_backup = []
        for file in files_to_backup:
            if os.path.exists(file):
                # Skip binary files that already have backup suffixes (they are already backups)
                if "_" in file and (file.endswith(".f00001") or ".f" in file):
                    print(f"[INFO] File '{file}' is already a backup binary file. Skipping backup.")
                    continue
                    
                backup_file = f"{file}{backup_suffix}"
                if os.path.exists(backup_file):
                    print(f"[INFO] Backup file '{backup_file}' already exists. Skipping backup of '{file}'.")
                else:
                    print(f"[INFO] Need to backup '{file}' to '{backup_file}'")
                    files_needing_backup.append(file)
            else:
                print(f"[INFO] Original file '{file}' does not exist. Skipping backup.")
        
        if files_needing_backup:
            backup_and_cleanup_files(files_needing_backup, backup_suffix, None, lambda x: None)
        else:
            print("[INFO] All files already backed up. No backup needed.")
        print("========================\n")

    print("==== [CONSOLIDATION] ====")
    consolidate_files(consolidate_file_patterns)
    print("========================\n")


# NOTE: For test runs only. Not used for production cases.
def prepare_par_for_restart(parameter_glob, single_job_time, case_preserving_adjust_par_file):
    """
    Update the .par file to set startFrom = rst_1cyl0.f00001 and endTime = single_job_time.
    This is intended for starting from scratch in test runs only.
    """
    import glob

    par_files = glob.glob(parameter_glob)
    if not par_files:
        raise FileNotFoundError(f"No .par file matching {parameter_glob} found.")
    par_file = par_files[0]
    case_preserving_adjust_par_file(par_file, {("GENERAL", "startFrom"): "rst_1cyl0.f00001", ("GENERAL", "endTime"): str(single_job_time)})
