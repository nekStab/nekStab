import os, sys

try:
    sys.path.insert(0, os.path.join(os.environ["NEKSTAB_SOURCE_ROOT"], "bin"))
except KeyError:
    raise EnvironmentError("NEKSTAB_SOURCE_ROOT environment variable is not set.")
from nekStab_tools import *

# ---- SIMULATION TARGET SETTINGS ----
target_end_time = 100.0  # Set your desired end time here
single_job_time = 10.0  # Set your desired field write interval here

pbs_script = "jz.pbs"
local_sh_script = "run_nohup.sh"

# ---- CASE-SPECIFIC SETTINGS ----
case_name = get_case_name_from_usr()
case_initial = case_name[0]
history_file = f"{case_name}.his"
parameter_glob = "*.par"

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

    final_time_in_his = check_time(history_file, which_time="final")
    initial_time_in_his = check_time(history_file, which_time="initial")
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

    if initial_time_in_his >= 0.0:
        if initial_time_in_his == 0 and final_time_in_his == 0:
            backup_suffix = f"_{initial_time_in_his}"
        else:
            backup_suffix = f"_{initial_time_in_his}_{final_time_in_his}"

        print(f"==== [BACKUP & CLEANUP] ====")
        print(f"[INFO] Backup suffix: {backup_suffix}")

        # Backup all files (restart, history, lift_drag, energy, enstrophy, logs) using the utility from nekStab_tools
        files_to_backup = [latest_restart_file, history_file, lift_drag_file, total_energy_file, total_enstrophy_file] + log_files
        backup_and_cleanup_files(files_to_backup, backup_suffix, None, lambda x: None)

    if velocity_file_time < target_end_time:
        next_time = velocity_file_time + single_job_time
        if velocity_file_time < single_job_time:
            next_time = single_job_time

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
        print(f"[INFO] Current time: {final_time_in_his} == Final time: {target_end_time}")
        print(f"[INFO] Current time is less than final time. No action taken.")
        print("======================\n")

    print("==== [CONSOLIDATION] ====")
    consolidate_files(consolidate_file_patterns)
    print("========================\n")
