# auto_restart: a DNS to a target time through a chain of restarts

A direct numerical simulation (DNS) that needs more time than one job may run (a wall-time limit of the queue, a node that is shared) is cut into jobs of `single_job_time` time units. Each job reads the last field file, runs to its `endTime`, and writes the field file again. At the end of each job `check_restart.py` backs up the outputs of that job, joins the histories of all jobs into one file, writes the next `startFrom` and `endTime` into the `.par` file and starts the next job. The chain stops when the time in the field file reaches `target_end_time`.

The two folders are the same case, a 2D flow around a circular cylinder (`1cyl`, 2128 elements, `lx1 = 6`, uniform inflow of 1, one probe at (5, 0, 0)), at two Reynolds numbers: `50/` is Re = 50 and `150/` is Re = 150. The folder name is the Reynolds number: `check_restart.py` writes `viscosity = -<folder name>` into the `.par` file at every restart (a negative viscosity in Nek5000 is 1/Re). The scripts are identical in the two folders; only `1cyl.par` (residual tolerances, viscosity) and `1cyl.usr` differ.

## Files of a case folder

| File | What it does |
|---|---|
| `1cyl.par`, `1cyl.usr`, `SIZE`, `1cyl.re2`, `1cyl.ma2` | The case, as for any nekStab run. `SIZE` sets `lpmin = 2`: use at least 2 ranks. |
| `1cyl.his` | Probe file. Before the first job it holds only the header (number of probes, then their coordinates). The solver appends the probe values to it (`hpts` in `userchk`, every 3rd step). |
| `rst_1cyl0.f00001` | Initial field of the first job (a field file of this mesh, time 0 here). Not changed by the chain. |
| `check_restart.py` | The logic. Settings at the top: `target_end_time` and `single_job_time`. Reads the time of the last field file, decides, backs up, joins, updates the `.par` file and submits. |
| `run.sh` | One job: writes `SESSION.NAME`, runs `mpiexec -n $NP ./nek5000 > logfile 2> logerror` (`NP` defaults to 5), then writes the exit status of `mpiexec` to `logstatus` and runs `python3 check_restart.py > logrestart`. |
| `run_nohup.sh` | Starts `run.sh` in the background of the shell (used where there is no Slurm). |
| `job.pbs` | Slurm script of one job. It calls `run.sh`. `check_restart.py` rewrites its `#SBATCH -J` line at each submission. |
| `clean_and_start.sh` | Starts a chain from scratch: cancels the Slurm jobs and stops the processes of this folder, deletes the previous outputs, sets `startFrom = rst_1cyl0.f00001` and the first `endTime`, rebuilds with `mks`, submits job 1. |
| `p_signals.py` | Plots the joined probe history (`1cyl.all`) and `lift_drag.all` and prints the peak frequencies: `python3 p_signals.py`. `p_signals_debug.py` is an older, shorter version of it. |

The Python functions are in `bin/nekStab_tools.py` of the nekStab tree.

## Use on your own case

1. Copy the case into a folder named by its Reynolds number, for example `example/mycase/200/`, with `.par`, `.usr`, `SIZE`, mesh, the initial field file and the probe file `<case>.his`. Copy `check_restart.py`, `run.sh`, `run_nohup.sh`, `job.pbs`, `clean_and_start.sh` and `p_signals.py` from here. The case name is the name of the `.usr` file and the scripts take it from there. If the folder name is not a number, `check_restart.py` leaves the viscosity of the `.par` file as it is.
2. In the `.par` file set `startFrom` (the first job starts from `rst_<case>0.f00001`, which `clean_and_start.sh` writes), `stopAt = endTime`, `writeControl = runTime` and a `writeInterval` larger than `single_job_time`, so that each job writes its field file once, at its end. Set `[PRESSURE] writeToFieldFile = yes`: see "Measured result" for the reason.
3. In `check_restart.py` set `target_end_time` (final time of the whole run) and `single_job_time` (time covered by one job; the last job is shortened to reach the target exactly). Choose `single_job_time` so that one job ends well inside the wall-time limit: a job that is killed by the limit writes no field file, and its whole segment is lost.
4. In `job.pbs` set the partition, `--ntasks` (at least `lpmin` of `SIZE`), `--mem` and the time limit; add the lines that load the compiler and MPI. `NEKSTAB_SOURCE_ROOT` must be set in the environment of the job: `check_restart.py` imports `$NEKSTAB_SOURCE_ROOT/bin/nekStab_tools.py`. `python3` needs `numpy`. Without Slurm, set the same variables in your shell and `NP` for the number of ranks.
5. Start: `./clean_and_start.sh`. It needs `NEKSTAB_SOURCE_ROOT`, `NEK_SOURCE_ROOT` and `mks` in the environment. It cancels the pending, running, suspended and completing Slurm jobs of this folder and stops its processes, then deletes the outputs of an earlier chain in this folder (`<case>0.f00*`, `<case>.his_*`, `lift_drag.*`, `total_energy.*`, `total_enstrophy.*`, `log*`), which includes the backup field files, so copy what you want to keep first. With `sbatch` on the path it submits `job.pbs`, otherwise it runs `run_nohup.sh`. For the Slurm route `sbatch`, `squeue` and `scancel` must all be on the path: if `squeue` is there and `scancel` is missing, or `squeue` or `scancel` fails, the script stops with an error before it deletes anything. If a job of the old chain is just ending, wait until the chain has stopped before you run it: a job that ends after the cleanup starts its successor in the cleaned folder. To continue an existing chain, do not use it.
6. Watch: `tail -f logfile` shows the running job; `squeue -u $USER` shows the queue (job names are `<first letter of the case name><Re>_<end time>`, cut at 8 characters); `tail logrestart_*` shows what the last `check_restart.py` decided. The progress of the whole run is the last line of `<case>.all` (time, then the probe values).
7. Change the target of a running or finished chain: edit `target_end_time` in `check_restart.py`. The next `check_restart.py` reads the new value. If the chain has already stopped at the old target, start the next job by running `python3 check_restart.py` by hand in the folder (with `NEKSTAB_SOURCE_ROOT` set). Tested: 20, then 30, then 40 time units.
8. Stop: cancel the queued or running job (`scancel <jobid>`); nothing resubmits it, because the next job is submitted only by the end of a job that has run to the end. Without Slurm, kill the `run.sh` shell first and then `nek5000`/`mpiexec` of this folder, in that order; if you kill only `nek5000`, `run.sh` goes on to `check_restart.py`, which stops with the error below. `clean_and_start.sh` also stops the processes of the folder, but then starts again from scratch.
9. Results: `<case>.all` (probe history of all jobs, one header, times increasing), `lift_drag.all`, the last field file `<case>0.f00001` and, for each job, the backups `<case>0.f00001_<t0>_<t1>`, `<case>.his_<t0>_<t1>`, `lift_drag.dat_<t0>_<t1>`, `logfile_<t0>_<t1>`, `logerror_<t0>_<t1>`, `logrestart_<t0>_<t1>`, with `t0` the time of the first and `t1` the time of the last probe record of that job.

## What check_restart.py does at the end of a job

1. Reads the time `T` from the header of the newest `<case>0.f0*` file (the original names are preferred to the backups) and the first and last time of `<case>.his`.
2. If `T = 0` and `logstatus` does not exist (a first start by hand from a field file at time 0, before any solver has run), it moves old `*.dat*`, `*.his_*`, `log*` and `*.all` into `old_dat/` and resets the history file.
3. Looks for a running job: the Slurm or PBS queue for the job name of the next job, and a `nek5000` process of this user whose working directory is this folder. If one runs, it only joins the histories.
4. Stops with an error if the job that just ended did not reach its `endTime` (see "Failure").
5. If `T < target_end_time`: sets `startFrom` to the newest field file, `endTime` to `T + min(single_job_time, target_end_time - T)` and the viscosity to `-Re`; backs up the field file, history, `lift_drag.dat`, `total_energy.dat`, `total_enstrophy.dat` and the three log files with the suffix `_<t0>_<t1>` (history and logs are reset or moved, so the next job starts with an empty history); joins all `*_<t0>_<t1>` histories into `<case>.all` (the header of the first, the data of all) and the `lift_drag`, `total_energy`, `total_enstrophy` files into `.all` files; submits the next job (`sbatch job.pbs` if `squeue` exists, else `bash run_nohup.sh`; `qstat` is used only to look for a running job).
6. If `T >= target_end_time`: backs up and joins the last job, and does nothing else.

## Failure

If a job dies (crash, MPI error, out of memory, a rank count below `lpmin`), `run.sh` still calls `check_restart.py`, after it has written the exit status of the solver to `logstatus`. The file tells `check_restart.py` that a solver has run: a field file at time 0 together with `logstatus` is a first job that wrote nothing, not a first start, so the failure check below applies to it too (without it the first job was taken for a first start, its logs were moved to `old_dat/` and the same job was submitted again, for ever; reproduced with a stand-in solver that wrote probe records and died: 55 jobs in 8 s). Before the check described above was added, the field file was still the one of the previous job, so the script saw an old time below the target, found nothing new to back up and submitted the same job again, for ever (observed with `NP=1`, which is below `lpmin = 2`: a new job every few seconds). Now `check_restart.py` compares the time in the field file with the `endTime` of the `.par` file, prints

```
[ERROR] Time in 1cyl0.f00001 is 30.0, below the endTime 35.0 of 1cyl.par.
```

and a second line with the exit status of the solver, and stops with exit code 1, before it changes any file. Read `logfile` and `logerror`, remove the cause, and repeat the job with `./run_nohup.sh` or `sbatch job.pbs`: the `.par` file is already set for it. Tested: after the failure above, `./run_nohup.sh` continued the chain to the target with a history that is complete and strictly increasing in time (1335 records from 0.03 to 40.01, spacing 0.02 to 0.03).

## Limits

- The chain needs the field file to carry everything the next job needs. Nek5000 restarts from the velocity (and the pressure and temperature, if written) at one time level. The time-stepping history (previous steps for BDF/extrapolation) is not stored, so each restart begins with the low-order start-up steps: this gives the small step of the next section. The nekStab filters and the forcing of an SFD run, if the `.usr` file uses them, also start again.
- `target_end_time` is reached exactly with the variable time step of the shipped `.par`. With a fixed `dt` each job ends at the first step at or after its `endTime`: a job meant to end at 5.0 ended at 5.01 with `dt = 0.01`, and the next job started from there. Choose `single_job_time` as a multiple of `dt`.
- Each job keeps a full copy of the field file (1.2 MB here, with `lx1 = 6` and 2128 elements) in addition to the logs and histories. A long chain of a large case fills the disk: delete the old `<case>0.f00001_<t0>_<t1>` once the next ones are good. Keep the newest.
- Only the probe history, `lift_drag`, `total_energy` and `total_enstrophy` are joined. Other files written by the `.usr` file need their own handling in `check_restart.py` (`consolidate_file_patterns`).
- Two job names that agree in their first 8 characters are taken for each other when the queue is searched. The name is `<first letter of the case name><Re>_<end time>`: the cut removes the end of the end time, so for Re = 150 the end times 100 and 1000 both give `1150_100`. One chain never has both in the queue at once; two chains of the same account can.
- Intermediate field files written by `writeInterval` shorter than `single_job_time` are not used and not backed up. They accumulate.
- `check_restart.py` cannot know a job that was cancelled while it waited in the queue: nothing runs and nothing resubmits it. Start it again by hand.
- Tested on one machine with Slurm (one node, partition `LocalQ`) and without Slurm (`squeue` hidden from the path), Linux, gfortran, Intel MPI, 5 ranks. Not tested: PBS, several nodes, a restart with a different number of ranks (Nek5000 supports it, nothing here prevents it), 3D cases.

## Measured result

Case `50/` (Re = 50), 5 ranks, target 20 time units, `single_job_time` = 5 (4 jobs), compared with one run of 20 time units from the same initial field and the same `.par`.

Shipped `.par` (variable time step, `targetCFL = 0.5`):

- Time continues without a jump: the joined history has 664 records from 0.030 to 20.000, strictly increasing, spacing 0.0088 to 0.0303, the same number of records as the single run. At each restart the last record of one job (5.000, 10.000, 15.000) is followed by a record 0.030 later.
- The probe values (u, v at (5, 0, 0)) differ from the single run by at most 2.6e-5 (u) and 4.3e-5 (v), against peak-to-peak amplitudes of 0.055 and 0.45. This includes the error of interpolating the single-run history, whose record times differ, to the times of the chained run; use the fixed-`dt` test below for the restart error itself.
- The chain stopped at 20 and submitted no further job.

Fixed time step (`dt = 0.01`, `variableDt = no`, so the record times of both runs are identical, 667 records, maximum time difference 0), local run without Slurm:

| Quantity | Result |
|---|---|
| First job (0 to 5.01) against the single run | identical to 3e-8 |
| Change of the difference to the single run across a restart (last record before to first record after) | 1.2e-6 to 4.9e-6 in u and 4.5e-6 to 1.1e-5 in v, at the three restarts; the differences add up to the next row |
| Difference at t = 20 (after 3 restarts) | 2.1e-5 in u, 2.5e-5 in v, i.e. 3.8e-4 and 5.6e-5 of the peak-to-peak amplitude |
| Kink in the curve | the second difference of the probe signal at the restarts (4e-5 to 1.4e-4) stays below its largest value in the run (3.8e-4 for u, 7.6e-4 for v); no kink is visible in the signal, only the offset above |
| Same, with `writeDoublePrecision = yes` | the same numbers to 3 digits: the offset is not a rounding error of the single-precision field file |
| Force coefficients with `writeToFieldFile = no` (restart file without pressure) | at the first record after each restart the drag differs from the single run by 7.5e-3 (0.5 % of 1.437), 6.0e-3 at the next, and the difference falls below 1e-4 after 11 records (0.4 time units) |
| Force coefficients with `writeToFieldFile = yes` | largest difference over the run: 4.5e-6 in drag, 1.7e-5 in lift |

The offset of a few 1e-6 per restart comes from the start of the time stepper from a single time level. The cause of the force transient is the missing pressure (a restart file without pressure starts the run with p = 0); with the pressure in the field file the transient is gone. This is why the `.par` files in `50/` and `150/` set `writeToFieldFile = yes`.

Chain on Slurm (`150/`, Re = 150, the committed scripts with `target_end_time` set to 12 and `single_job_time` to 5 (committed: 200 and 50), jobs of 5, 5 and 2 time units, 5 ranks): the viscosity was set to -150 by the script, the joined history has 411 records from 0.030 to 12.000, strictly increasing, with the records at 5.000 and 10.000 followed by 5.029 and 10.029, and the last job was shortened to end at 12.000. Each job was submitted by the previous one; no job was submitted after the target. The kill of a `check_restart.py` started in the same folder by `clean_and_start.sh` was tested with a stand-in process: it was stopped, and a process of the same name in another folder was not.

Timing on this machine (24 cores, 5 ranks used): 4.2 s for a job of 5 time units, 15.2 s for 20 time units in one job (7.5 ms per step). On Slurm the next job started 25 to 30 s after its submission on this machine (observed in four jobs); on a loaded cluster the wait in the queue is added to every job.
