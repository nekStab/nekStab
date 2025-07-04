"""
-------------------------------------------------------------------------------
 nekStab_tools.py - General Automation and Utility Routines for NekStab Workflows
-------------------------------------------------------------------------------

 This module contains general-purpose functions for file management, job submission,
 backup, consolidation, and other utilities used in NekStab auto-restart and
 simulation consolidation workflows. All case-specific logic should be kept in the
 main driver scripts (e.g., restart_manager.py). This module is intended to be imported
 and reused across different cases and automation tasks.

 Author: Ricardo Frantz
 Last updated: 2025-07-04
-------------------------------------------------------------------------------
"""
import os, re, shutil, subprocess, glob, sys, datetime
import numpy as np


def dist(rec, N, s):
    y = np.cos(np.pi * np.arange(N) / (N - 1))
    x = rec - s * y / np.sqrt(1 + ((s / (N - 1)) ** 2) - y**2)
    x = np.round(x, 0)
    x = list(map(str, x.astype(int)))
    print("range distribution around ", rec, x)
    return x


def build_nek(folder_path, usr_file, verbose=False):
    from subprocess import Popen, PIPE, STDOUT
    from pathlib import Path

    source_root = os.environ["NEKSTAB_SOURCE_ROOT"]
    cwd = os.getcwd() + "/" + folder_path
    print("Compiling nek5000...")
    print(f'    Using working directory "{cwd}"')
    print(f'    Using .usr file "{usr_file}"')
    my_env = os.environ.copy()
    makenek_in = Path(source_root) / "bin" / "mks"
    logfile = Path(cwd) / "build.log"
    proc = Popen([makenek_in, "clean"], cwd=cwd, env=my_env, stdin=PIPE, text=True)
    proc.communicate(input="Y\n")
    proc.wait()
    proc = Popen([makenek_in, usr_file], cwd=cwd, env=my_env, stdin=PIPE, stderr=STDOUT)
    proc.wait()

    if proc.returncode != 0:
        with open(logfile, "r") as file:
            text = file.read()
        print(text)
        exit(-1)


def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.exists(d):
            try:
                shutil.rmtree(d)
            except Exception as e:
                print(e)
                os.unlink(d)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def cSZ(infile, outfile, params):
    with open(infile, "r") as f:
        lines = f.readlines()
    # Substitute all the variables
    for key, value in params.items():
        print("Modf.", key, "to", value, "in", infile)
        if value:
            lines = [
                re.sub(
                    r"(.*\bparameter\b.*\b{0} *= *)\S+?( *[),])".format(key),
                    r"\g<1>{0}\g<2>".format(value),
                    l,
                    flags=re.I,
                )
                for l in lines
            ]
    with open(outfile, "w") as f:
        f.writelines(lines)


def check_keyword_in_file(folder_path, file_name, keyword):
    file_path = os.path.join(folder_path, file_name)
    if os.path.exists(file_path):
        with open(file_path, "r") as file:
            file_contents = file.read()
            if keyword in file_contents:
                print(f"The keyword '{keyword}' was found in the file '{file_name}'.")
                return True
            else:
                print(f"The keyword '{keyword}' was not found in the file '{file_name}'.")
                return False
    else:
        print(f"The file '{file_name}' does not exist in the folder '{folder_path}'.")
        return False


def write_to_file(filename, content, target_line):
    with open(filename, "r") as file:
        lines = file.readlines()
    # content_str = ' '.join(map(str, content))
    lines[target_line - 1] = content + "\n"
    with open(filename, "w") as file:
        file.writelines(lines)


def compile_neks(casename):
    command = ["makeneks", casename]
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    output, error = process.communicate(input=b"y\n")  # Send 'y' and newline ('\n')

    if process.returncode != 0:  # If the command failed...
        print(f"Compilation for {casename} failed.")
        if error:
            print("Error message:", error.decode("utf-8"))
    else:
        print("Command output:", output.decode("utf-8"))


def c_pf(infile, outfile, opts):
    import configparser

    parfile = configparser.ConfigParser()
    parfile.read(infile)
    for section, name_vals in opts.items():
        for name, val in name_vals.items():
            print("In", section, ":", name, "set to", val, "in", infile)
            parfile.set(section, name, val)
    with open(outfile, "w") as f:
        parfile.write(f)


def check_last_value(filename, tolerance):
    with open(filename, "r") as file:
        last_line = file.readlines()[-1]  # Read the last line of the file
    last_value = float(last_line.split()[-1])  # Split the line into parts and convert the last part to float
    return last_value < tolerance


def delete_files(pattern):
    [os.remove(file) for file in glob.glob(pattern)]


def copy_bf(file, residu_file, oldbfs, tolerance):
    if os.path.isfile(file) and os.path.isfile(residu_file):
        print(f"Found '{file}' and '{residu_file}'")
        with open(residu_file, "r") as f:
            lines = f.readlines()
            last_value = float(lines[-1].split()[-1])
            if last_value < tolerance:
                print(f" Last value of '{residu_file}' is {last_value} < {tolerance}")
                curr_folder = os.path.basename(os.getcwd())
                dest_path = os.path.join(oldbfs, curr_folder)
                print(f" Copying '{os.getcwd() + '/' + file}' to '{dest_path}'")
                shutil.copy(file, dest_path)


def check_job_status(job_name):
    result = subprocess.run(["squeue", "-u", "rvpo014"], capture_output=True, text=True)
    job_status = job_name in result.stdout
    if job_status:
        print(f"Job {job_name} is running.")
        print("")
        print("")
        sys.exit()
    else:
        print(f"Job {job_name} is not running.")
    return job_status


def check_job_exist(job_name):
    result = subprocess.run(["squeue", "-u", "rvpo014"], capture_output=True, text=True)
    # print(result)
    job_status = job_name in result.stdout
    if job_status:
        print(f"Job {job_name} exist.")
    else:
        print(f"Job {job_name} does not exist.")
    return job_status


def submit_job(job_name, job_script, working_directory, batch_file):
    """Submit a job: use SLURM if available, otherwise run locally with a script."""
    if shutil.which("squeue") is not None:
        if not check_job_exist(job_name):
            resubmit_job(batch_file, working_directory, job_name)
    else:
        print(f"[INFO] squeue not found, running locally with {job_script}")
        subprocess.call(["bash", job_script])


def resubmit_job(pbs_file, folder_path, job_name):
    print(f" Job name: {job_name}")
    pbs_file_path = os.path.join(folder_path, pbs_file)
    print(f" Adjusting file '{pbs_file_path}'")
    with open(pbs_file_path, "r") as file:
        filedata = file.read()
    match = re.search(r"(#SBATCH -J )\S+", filedata)
    if match:
        old_line = match.group(0)  # Get the entire matched line
        new_line = f"#SBATCH -J {job_name}"
        filedata = filedata.replace(old_line, new_line)
        print(f"Replacing '{old_line}' with '{new_line}'")
    with open(pbs_file_path, "w") as file:
        file.write(filedata)
    current_dir = os.getcwd()
    os.chdir(folder_path)
    try:
        result = subprocess.run(["sbatch", pbs_file], check=True, stdout=subprocess.PIPE)
        print(f"Job {job_name} submitted successfully.")
        print("Command output:", result.stdout.decode("utf-8"))  # Print the command output
    except subprocess.CalledProcessError:
        print(f"Job submission for {job_name} failed.")
    os.chdir(current_dir)


def get_closest_filename(target_dir, reference):
    file_names = os.listdir(target_dir)
    reference = int(reference)  # assuming reference is a numerical value

    print(f" Looking for closest match to {reference} in {target_dir}")
    print(f" Found {file_names} files in {target_dir}")

    closest_diff = float("inf")
    closest_filename = None

    for file_name in file_names:
        # assuming filenames are numbers
        file_number = int(file_name)
        diff = abs(reference - file_number)

        if diff < closest_diff:
            closest_diff = diff
            closest_filename = file_name

    print(f" Closest match to {reference} is {closest_filename} with a difference of {closest_diff}")
    return closest_filename


def adjust_and_submit_job(folder_path, job_name):
    """Adjust the PBS file and submit the job."""
    pbs_file_path = os.path.join(folder_path, pbs_file)
    with open(pbs_file_path, "r") as file:
        filedata = file.read()

    filedata = filedata.replace("JOBNAME", job_name)
    filedata = filedata.replace("CASENAME", f'"{cn}"')

    with open(pbs_file_path, "w") as file:
        file.write(filedata)

    try:
        result = subprocess.run(["sbatch", pbs_file], check=True, stdout=subprocess.PIPE)
        print(f"Job {job_name} submitted successfully.")
        print("Command output:", result.stdout.decode("utf-8"))
    except subprocess.CalledProcessError:
        print(f"Job submission for {job_name} failed.")


def check_time(filename, which_time="final"):
    print(f"Checking {which_time} time in {filename}")
    with open(filename, "r") as file:
        lines = file.readlines()
        # print(lines)
        if not lines:  # if the file is empty
            raise ValueError(f"The file '{filename}' is empty.")
        num = int(lines[0].split()[0]) + 1  # get the number on the first line
        print(f"  Number of lines: {num}")
        print(f"  Number of lines in file: {len(lines)}")
        if len(lines) > num:
            if which_time == "final":
                line = lines[-1]
            elif which_time == "initial":
                line = lines[num]  # skip the first line + 'num' lines
            value = float(line.split()[0])
        else:
            print(f"  File '{filename}' does not contain enough lines.")
            value = 0.0
    print(f"  Time: {value}")
    return value


def reset_his_file(filename):
    temp_filename = filename + ".tmp"
    with open(filename, "r") as file, open(temp_filename, "w") as temp_file:
        num = int(next(file).split()[0])  # read the first line
        temp_file.write(str(num) + "\n")  # write the first line to the temp file
        for _ in range(num):  # copy 'num' lines
            temp_file.write(next(file))
    os.remove(filename)  # delete the original file
    shutil.move(temp_filename, filename)  # move the temp file to the original file's location


def extract_time_from_binary(filename):
    # Check if the file exists
    # if not os.path.isfile(filename):
    #   raise FileNotFoundError(f"The file {filename} does not exist.")
    #   sys.exit()
    try:
        result = subprocess.run(["head", "-1", filename], stdout=subprocess.PIPE)
        first_line = result.stdout.decode("utf-8", errors="ignore")
        words = first_line.split()

        # Check if words[7] is empty
        # if len(words) < 8 or not words[7]:
        #   raise ValueError(f"The file {filename} does not contain a value at index 7.")

        # try:
        time = float(words[7])  # adjust this index based on where the value is located
    except:
        time = np.NaN
    return time


def append_files(pattern, output_file):
    print(f"[INFO] Consolidating files matching '{pattern}' into '{output_file}':")
    files = sorted(glob.glob(pattern), key=lambda f: [float(num) for num in re.findall(r"\d+\.\d+", f)])
    if not files:
        print(f"  [WARN] No files matching pattern '{pattern}' found. Skipping consolidation for '{output_file}'.")
        return
    with open(output_file, "w") as f:
        print(f"  [INFO] Created empty output file: {output_file}")
        pass  # create the file
    with open(output_file, "w") as outfile:
        for i, file in enumerate(files):
            print(f"    [INFO] File {i}: {file}")
            if ".his_" in file:
                with open(file, "r") as infile:
                    lines = infile.readlines()  # Read all lines into a list
                    if i == 0:
                        num = int(lines[0].split()[0])
                        outfile.writelines(lines[:-num])
                        print(f"      [INFO] Found '.his_', skipping {num} lines from the end")
                    elif i > 0 and i < len(files) - 1:
                        outfile.writelines(lines[num + 1 : -num])
                    else:
                        outfile.writelines(lines[num + 1 :])
            else:
                with open(file, "r") as infile:
                    shutil.copyfileobj(infile, outfile)



def find_latest_final_dns_file(cn, verbose=False):
    pattern = f"{cn}0.f0*"
    files = glob.glob(pattern)
    if not files:
        if verbose:
            print(f"No files found matching pattern: {pattern}")
        return None
    # Sort by modification time (most recent last)
    files_sorted = sorted(files, key=os.path.getmtime)
    if verbose:
        print(f"Found {len(files_sorted)} files matching pattern '{pattern}':")
        for f in files_sorted:
            mtime = os.path.getmtime(f)
            print(f"  {f}  (modified: {datetime.datetime.fromtimestamp(mtime)})")
        print(f"Selected latest file: {files_sorted[-1]}")
    return files_sorted[-1]


def get_case_name_from_usr():
    usr_files = glob.glob("*.usr")
    if not usr_files:
        raise FileNotFoundError("No .usr file found in the current directory.")
    # If multiple .usr files exist, pick the first or implement your own logic
    return os.path.splitext(os.path.basename(usr_files[0]))[0]


def adjust_par_file(pf, param_dict):
    """Adjust parameters in the .par file for the simulation using a dictionary of section/key: value."""
    for (section, key), value in param_dict.items():
        c_pf(pf, pf, {section: {key: value}})


def case_preserving_adjust_par_file(filename, param_dict):
    """
    Update parameters in a .par file, preserving case, spaces, and comments.
    param_dict: { (section, key): value }
    """
    import re
    with open(filename, "r") as f:
        lines = f.readlines()

    current_section = None
    section_re = re.compile(r'^\s*\[(.+?)\]\s*$')
    key_re_template = r'^(\s*{key}\s*=\s*)([^#;\n]*)(.*)$'

    # Build a lookup for faster access
    param_lookup = {(sect, key): val for (sect, key), val in param_dict.items()}
    updated = set()

    for i, line in enumerate(lines):
        section_match = section_re.match(line)
        if section_match:
            current_section = section_match.group(1)
            continue

        if current_section is None:
            continue

        # For each key in this section, try to match and replace
        for (sect, key), value in param_lookup.items():
            if sect == current_section:
                key_re = re.compile(key_re_template.format(key=re.escape(key)))
                m = key_re.match(line)
                if m:
                    # Reconstruct line: keep original spaces and comments
                    lines[i] = f"{m.group(1)}{value}{m.group(3)}\n"
                    updated.add((sect, key))
                    break  # Only one key per line

    # Optionally: warn if any keys were not found
    not_found = set(param_lookup.keys()) - updated
    if not_found:
        print(f"Warning: The following parameters were not found and not updated: {not_found}")

    with open(filename, "w") as f:
        f.writelines(lines)

# Usage:
# case_preserving_adjust_par_file('1cyl.par', {('GENERAL', 'startFrom'): 'rst_1cyl0.f00001'})


def check_required_files_exist(file_list):
    """Check that all files in file_list exist. Exit if any are missing."""
    missing = [f for f in file_list if not os.path.exists(f)]
    if missing:
        print(f"Missing required files: {', '.join(missing)}")
        sys.exit(1)


def backup_and_cleanup_files(files, backup_suffix, history_file, reset_his_file_func):
    """
    Binary file: copy for backup.
    .his file: copy for backup, then reset using provided function.
    All other files: MOVE (not copy) to backup name.
    files: list of files to backup (binary, his, then rest)
    backup_suffix: string to append to backup files
    history_file: the main .his file to reset
    reset_his_file_func: function to reset the history file
    """
    import shutil, os
    if not files:
        print("[WARN] No files provided for backup.")
        return
    # Binary file (first in list): copy
    binary_file = files[0] if len(files) > 0 else None
    his_file = files[1] if len(files) > 1 else None
    rest_files = files[2:] if len(files) > 2 else []

    # Backup binary file
    if binary_file and os.path.exists(binary_file):
        backup_name = binary_file + backup_suffix
        print(f"  [INFO] Binary file: {binary_file} → {backup_name}")
        shutil.copy(binary_file, backup_name)
    elif binary_file:
        print(f"  [WARN] Binary file '{binary_file}' not found, skipping backup.")

    # Backup and reset .his file
    if his_file and os.path.exists(his_file):
        backup_name = his_file + backup_suffix
        print(f"  [INFO] .his file: {his_file} → {backup_name} (reset after backup)")
        shutil.copy(his_file, backup_name)
        reset_his_file_func(his_file)
    elif his_file:
        print(f"  [WARN] .his file '{his_file}' not found, skipping backup and reset.")

    # Move all other files
    for f in rest_files:
        if os.path.exists(f):
            backup_name = f + backup_suffix
            print(f"  [INFO] Moving: {f} → {backup_name}")
            shutil.move(f, backup_name)
        else:
            print(f"  [WARN] File '{f}' not found, skipping move.")


def consolidate_files(pattern_output_list):
    """
    Consolidate files matching patterns into single output files.
    pattern_output_list: list of (pattern, output_file) tuples
    """
    for pattern, output_file in pattern_output_list:
        append_files(pattern, output_file)
