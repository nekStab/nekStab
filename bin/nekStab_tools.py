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
    """Generates N integer values distributed around a reference point using cosine spacing.
    Args: rec (float), N (int), s (float)
    Returns: list of string integers
    """
    y = np.cos(np.pi * np.arange(N) / (N - 1))
    x = rec - s * y / np.sqrt(1 + ((s / (N - 1)) ** 2) - y**2)
    x = np.round(x, 0)
    x = list(map(str, x.astype(int)))
    print("range distribution around ", rec, x)
    return x


def build_nek(folder_path, usr_file, verbose=False):
    """Cleans and compiles Nek5000 from .usr file, exits on failure.
    Args: folder_path (str), usr_file (str), verbose (bool)
    Returns: None
    """
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
    """Copies directory tree from src to dst, overwriting existing files.
    Args: src (str), dst (str), symlinks (bool), ignore (callable)
    Returns: None
    """
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
    """Modifies Fortran parameter statements using regex substitution.
    Args: infile (str), outfile (str), params (dict)
    Returns: None
    """
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
    """Checks if keyword exists in file, returns True if found.
    Args: folder_path (str), file_name (str), keyword (str)
    Returns: bool
    """
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
    """Replaces specific line (1-indexed) in file with new content.
    Args: filename (str), content (str), target_line (int)
    Returns: None
    """
    with open(filename, "r") as file:
        lines = file.readlines()
    # content_str = ' '.join(map(str, content))
    lines[target_line - 1] = content + "\n"
    with open(filename, "w") as file:
        file.writelines(lines)


def compile_neks(casename):
    """Compiles Nek5000 case using makeneks script, automatically answers 'y' to prompts.
    Args: casename (str)
    Returns: None
    """
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
    """Modifies parameters in .par file using ConfigParser.
    Args: infile (str), outfile (str), opts (dict of {section: {param: value}})
    Returns: None
    """
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
    """Checks if last numerical value in file is below tolerance.
    Args: filename (str), tolerance (float)
    Returns: bool (True if last value < tolerance)
    """
    with open(filename, "r") as file:
        last_line = file.readlines()[-1]  # Read the last line of the file
    last_value = float(last_line.split()[-1])  # Split the line into parts and convert the last part to float
    return last_value < tolerance


def delete_files(pattern):
    """Deletes all files matching glob pattern.
    Args: pattern (str) - glob pattern like '*.log' or 'temp_*'
    Returns: None
    """
    [os.remove(file) for file in glob.glob(pattern)]


def copy_bf(file, residu_file, oldbfs, tolerance):
    """Copies baseflow file if residual converged below tolerance.
    Args: file (str), residu_file (str), oldbfs (str), tolerance (float)
    Returns: None
    """
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
    """Checks if SLURM job is running and exits if found.
    Args: job_name (str)
    Returns: bool (True if running, exits program if found)
    """
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
    """Checks if SLURM job exists in queue.
    Args: job_name (str)
    Returns: bool (True if job exists)
    """
    result = subprocess.run(["squeue", "-u", "rvpo014"], capture_output=True, text=True)
    # print(result)
    job_status = job_name in result.stdout
    if job_status:
        print(f"Job {job_name} exist.")
    else:
        print(f"Job {job_name} does not exist.")
    return job_status


def submit_job(job_name, job_script, working_directory, batch_file):
    """Submits job using SLURM if available, otherwise runs locally with script.
    Args: job_name (str), job_script (str), working_directory (str), batch_file (str)
    Returns: None
    """
    if shutil.which("squeue") is not None:
        if not check_job_exist(job_name):
            resubmit_job(batch_file, working_directory, job_name)
    else:
        print(f"[INFO] squeue not found, running locally with {job_script}")
        subprocess.call(["bash", job_script])


def resubmit_job(pbs_file, folder_path, job_name):
    """Resubmits SLURM job by adjusting job name in PBS file.
    Args: pbs_file (str), folder_path (str), job_name (str)
    Returns: None
    """
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
    """Finds filename closest to reference number in target directory.
    Args: target_dir (str), reference (str or int)
    Returns: str (closest filename)
    """
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
    """Adjusts PBS file and submits job with updated parameters.
    Args: folder_path (str), job_name (str)
    Returns: None
    """
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
    """Checks time value from history file, either initial or final.
    Args: filename (str), which_time (str) - 'final' or 'initial'
    Returns: float (time value)
    """
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
    """Resets history file by keeping only header and coordinate lines.
    Args: filename (str)
    Returns: None
    """
    temp_filename = filename + ".tmp"
    with open(filename, "r") as file, open(temp_filename, "w") as temp_file:
        num = int(next(file).split()[0])  # read the first line
        temp_file.write(str(num) + "\n")  # write the first line to the temp file
        for _ in range(num):  # copy 'num' lines
            temp_file.write(next(file))
    os.remove(filename)  # delete the original file
    shutil.move(temp_filename, filename)  # move the temp file to the original file's location


def extract_time_from_binary(filename):
    """Extracts time value from binary file header.
    Args: filename (str)
    Returns: float (time value or NaN if extraction fails)
    """
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
    """Consolidates files matching pattern into single output file.
    Args: pattern (str), output_file (str)
    Returns: None
    """
    print(f"[INFO] Consolidating files matching '{pattern}' into '{output_file}':")
    
    # Sort files by time range suffix for chronological order
    files = sorted(glob.glob(pattern), key=lambda f: float(f.split('_')[-1]))
    
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
                        # First file: include header (2 lines) + all data (no footer to skip)
                        outfile.writelines(lines)
                        print(f"      [INFO] Found '.his_', including all data from first file")
                    else:
                        # Subsequent files: skip header (2 lines), include all data
                        outfile.writelines(lines[2 :])
            else:
                with open(file, "r") as infile:
                    shutil.copyfileobj(infile, outfile)



def find_latest_final_dns_file(cn, verbose=False):
    """Finds latest DNS final file, prioritizing original files over backup files.
    Args: cn (str), verbose (bool)
    Returns: str (latest filename or None if not found)
    """
    pattern = f"{cn}0.f0*"
    files = glob.glob(pattern)
    if not files:
        if verbose:
            print(f"No files found matching pattern: {pattern}")
        return None
    
    # Separate original files from backup files
    original_files = []
    backup_files = []
    
    for file in files:
        # Original files have pattern: case0.f00001, case0.f00002, etc.
        # Backup files have pattern: case0.f00001_time1_time2, etc.
        if '_' in file:
            backup_files.append(file)
        else:
            original_files.append(file)
    
    if verbose:
        print(f"Found {len(files)} files matching pattern '{pattern}':")
        for f in sorted(files, key=os.path.getmtime):
            mtime = os.path.getmtime(f)
            file_type = "[ORIGINAL]" if '_' not in f else "[BACKUP]"
            print(f"  {f}  (modified: {datetime.datetime.fromtimestamp(mtime)}) {file_type}")
    
    # Prioritize original files
    if original_files:
        # Sort original files by modification time and select the latest
        original_files_sorted = sorted(original_files, key=os.path.getmtime)
        selected_file = original_files_sorted[-1]
        if verbose:
            print(f"Selected latest original file: {selected_file}")
        return selected_file
    elif backup_files:
        # If no original files, fall back to backup files
        backup_files_sorted = sorted(backup_files, key=os.path.getmtime)
        selected_file = backup_files_sorted[-1]
        if verbose:
            print(f"No original files found, selected latest backup file: {selected_file}")
        return selected_file
    else:
        if verbose:
            print("No files found")
        return None


def get_case_name_from_usr():
    """Gets case name from .usr file in current directory.
    Args: None
    Returns: str (case name without extension)
    """
    usr_files = glob.glob("*.usr")
    if not usr_files:
        raise FileNotFoundError("No .usr file found in the current directory.")
    # If multiple .usr files exist, pick the first or implement your own logic
    return os.path.splitext(os.path.basename(usr_files[0]))[0]


def adjust_par_file(pf, param_dict):
    """Adjusts parameters in .par file using dictionary of section/key: value.
    Args: pf (str), param_dict (dict)
    Returns: None
    """
    for (section, key), value in param_dict.items():
        c_pf(pf, pf, {section: {key: value}})


def case_preserving_adjust_par_file(filename, param_dict):
    """Updates parameters in .par file, preserving case, spaces, and comments.
    Args: filename (str), param_dict (dict of {(section, key): value})
    Returns: None
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
                    # Reconstruct line: ensure exactly one space before comment if present
                    comment = m.group(3)
                    if comment:
                        comment_stripped = comment.lstrip()
                        if comment_stripped.startswith('#'):
                            lines[i] = f"{m.group(1)}{value} {comment_stripped}\n"
                        else:
                            lines[i] = f"{m.group(1)}{value}{comment}\n"
                    else:
                        lines[i] = f"{m.group(1)}{value}\n"
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
    """Checks that all files exist, exits if any are missing.
    Args: file_list (list)
    Returns: None
    """
    missing = [f for f in file_list if not os.path.exists(f)]
    if missing:
        print(f"Missing required files: {', '.join(missing)}")
        sys.exit(1)


def backup_and_cleanup_files(files, backup_suffix, history_file, reset_his_file_func):
    """Backs up and cleans files with different strategies per file type.
    Args: files (list), backup_suffix (str), history_file (str), reset_his_file_func (callable)
    Returns: None
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
    """Consolidates files matching patterns into single output files.
    Args: pattern_output_list (list of tuples)
    Returns: None
    """
    for pattern, output_file in pattern_output_list:
        append_files(pattern, output_file)

# Read initial and final time as both float (for logic) and string (for suffix)
from decimal import Decimal, localcontext

def format_time_str(time_str):
    """Formats time string for filenames.
    - If integer, returns only the integer digits.
    - If decimal and has only one nonzero digit after the decimal (e.g., 0.1, 0.2), use the shortest representation.
    - Otherwise, use exactly 9 digits after the decimal (zero-padded).
    Args: time_str (str)
    Returns: str (formatted time)
    """
    try:
        with localcontext() as ctx:
            ctx.prec = 16
            val = Decimal(time_str)
            if val == val.to_integral():
                return str(val.to_integral())
            # Check if can be represented as a short decimal (e.g., 0.1, 0.2, 0.5)
            s = format(val.normalize(), 'f')
            if '.' in s and len(s.rstrip('0').split('.')[-1]) == 1:
                # Only one nonzero digit after decimal
                s = s.rstrip('0').rstrip('.')
                return s
            # Otherwise, use 9 digits after decimal, zero-padded
            return f"{val:.9f}"
    except Exception:
        return time_str

def get_time_str_and_float(filename, which_time):
    """Gets time value from file as both formatted string and float.
    Args: filename (str), which_time (str)
    Returns: tuple (formatted time string, time float)
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
        if not lines:
            raise ValueError(f"The file '{filename}' is empty.")
        num = int(lines[0].split()[0]) + 1
        if len(lines) > num:
            if which_time == "final":
                line = lines[-1]
            elif which_time == "initial":
                line = lines[num]
            time_str = line.split()[0]
            time_float = float(time_str)
        else:
            time_str = "0"
            time_float = 0.0
    return format_time_str(time_str), time_float

def periodogram_rfft(x, fs, scaling="spectrum"):
    """Computes power spectrum of real-valued signal using periodogram.
    Args: x (array), fs (float), scaling (str)
    Returns: tuple (frequencies, power spectrum)
    """
    from scipy import signal

    # Efficient real FFT-based periodogram for real-valued signals
    freqs, psd = signal.periodogram(x, fs, scaling=scaling)
    return freqs, psd


def find_peaks(st, psd, threshold=0.01):
    """Finds peaks in power spectrum above threshold fraction of maximum.
    Args: st (array), psd (array), threshold (float)
    Returns: tuple (peak frequencies, peak values)
    """
    from scipy import signal
    import numpy as np

    # Handle empty or all-zero input
    if len(psd) == 0 or not np.any(psd):
        return np.array([]), np.array([])

    # Detect peaks above threshold * max(psd)
    peak_indices = signal.find_peaks(psd, height=max(psd) * threshold)[0]

    return st[peak_indices], psd[peak_indices]

def interpolate_signal(t, y, num_points=None, t_new=None, fft_safe=True):
    """Interpolates signal to new time base with constant time step.
    Args: 
        t (array): Original time array
        y (array): Original signal array  
        num_points (int): Number of points for interpolation (if t_new not provided)
        t_new (array): New time array (if provided, overrides num_points)
        fft_safe (bool): If True, use max dt from original data for FFT-safe interpolation
    Returns: tuple (new time array, interpolated signal)
    """
    from scipy.interpolate import interp1d
    import numpy as np

    if len(t) < 2:
        # Cannot interpolate if there are not enough points
        return t, y

    if t_new is None:
        if fft_safe and len(t) > 2:
            # For FFT analysis with variable time steps: use max dt to avoid artificial high-freq content
            dt_orig = np.diff(t)
            max_dt = np.max(dt_orig)
            min_dt = np.min(dt_orig)
            mean_dt = np.mean(dt_orig)
            
            # Calculate number of points based on max_dt for conservative interpolation
            total_time = t[-1] - t[0]
            num_points_conservative = int(np.ceil(total_time / max_dt)) + 1
            
            print(f"Variable time step detected: min_dt={min_dt:.6f}, max_dt={max_dt:.6f}, mean_dt={mean_dt:.6f}")
            print(f"Time range: {t[0]:.6f} to {t[-1]:.6f} (duration: {total_time:.6f})")
            print(f"Using max_dt={max_dt:.6f} for FFT-safe interpolation ({num_points_conservative} points)")
            
            # Check if time array gets cropped due to max_dt constraint
            t_new_end = t[0] + (num_points_conservative - 1) * max_dt
            if t_new_end < t[-1] - max_dt/10:  # Allow small tolerance
                print(f"WARNING: Time array cropped from {t[-1]:.6f} to {t_new_end:.6f} due to max_dt constraint")
            
            t_new = np.arange(t[0], t[-1] + max_dt/2, max_dt)
        else:
            # Original behavior: uniform spacing with specified number of points
            if num_points is None:
                num_points = len(t)
            t_new = np.linspace(t.min(), t.max(), num_points)

    # Use linear interpolation for robustness, cubic can be unstable with noisy data
    f = interp1d(t, y, kind='linear', bounds_error=False, fill_value="extrapolate")
    y_new = f(t_new)

    return t_new, y_new


class LiftDragLoader(object):
    """Loads and processes lift/drag data from simulation files.
    Args: filename (str), flip (bool)
    Returns: LiftDragLoader object with t, dgx, dgy, dgz attributes
    """
    def __init__(self, filename, flip=False):
        print("Reading " + filename)
        data = []
        with open(filename, "r") as f:
            for line in f:
                parts = line.strip().split()
                # Skip lines that are too short or not data
                if len(parts) < 3 or any(c in line for c in ['C', 's', 't']) or line.startswith('#'):
                    continue
                try:
                    # Extract time from column 1 (index 1) if available, else column 0
                    time = float(parts[1]) if len(parts) > 1 else float(parts[0])
                    # Extract dragx and dragy from available columns
                    dragx = float(parts[2]) if len(parts) > 2 else 0.0
                    dragy = float(parts[5]) if len(parts) > 5 else (float(parts[3]) if len(parts) > 3 else 0.0)
                    # Try to get dragz if present (3D), else None
                    dragz = float(parts[8]) if len(parts) > 8 else None
                    data.append([time, dragx, dragy, dragz])
                except (ValueError, IndexError):
                    continue
        if not data:
            print(f"Warning: No data read from {filename}")
            self.t = self.dgx = self.dgy = self.dgz = np.array([])
            return

        d = np.transpose(data)
        self.t = d[0]
        self.dgx = d[1]
        self.dgy = d[2]
        self.dgz = d[3] if d.shape[0] > 3 and np.all([x is not None for x in d[3]]) else None
        if flip and self.dgz is not None:
            self.dgy, self.dgz = self.dgz, self.dgy

def zero_his(filename):
    """Zeros out a .his file by truncating it after the header section.

    This function reads the number of header lines from the first line of the file,
    then truncates the file to keep only the header and one subsequent line,
    effectively resetting the history data.

    Args:
        filename (str): The path to the .his file.
    """
    print(f"Zeroing history file: {filename}")
    try:
        with open(filename, 'r+') as f:
            lines = f.readlines()
            if not lines:
                print(f"Warning: File is empty, cannot zero: {filename}")
                return

            try:
                # First line should contain the number of header lines
                n = int(lines[0].strip())
            except (ValueError, IndexError):
                print(f"Error: Could not parse header line count from {filename}")
                return

            # Keep header (n lines) + column names (1 line)
            # The original bash script was sed 'n+2,$d', which keeps lines 1 to n+1
            f.seek(0)
            f.writelines(lines[:n + 1])
            f.truncate()
            f.write('\n')  # Add a blank line at the end

        print(f"Time zeroed in file: {filename}")
    except IOError as e:
        print(f"Error processing file {filename}: {e}")

def zero_binary(filename):
    """Zeros out the time in a Nek5000 binary file header (e.g., .fld).

    This function reads the first line of the specified file, finds a
    floating-point time value in scientific notation, and replaces it with zero.

    Args:
        filename (str): The path to the binary-format file.
    """
    print(f"Zeroing binary file: {filename}")
    try:
        with open(filename, 'r+') as f:
            lines = f.readlines()
            if not lines:
                print(f"Warning: File is empty, cannot zero: {filename}")
                return

            first_line = lines[0]
            # Pattern for Nek5000 time format: 0.1234567890123E+00
            pattern = r'0\.\d{13}E[+-]\d{2}'
            replacement = '0.0000000000000E+00'

            # Replace the pattern in the first line
            new_first_line, num_replacements = re.subn(pattern, replacement, first_line, count=1)

            if num_replacements > 0:
                lines[0] = new_first_line
                f.seek(0)
                f.writelines(lines)
                f.truncate()
                print(f"Time zeroed in file: {filename}")
            else:
                print(f"Time pattern not found in {filename}. File not modified.")

    except IOError as e:
        print(f"Error processing file {filename}: {e}")
