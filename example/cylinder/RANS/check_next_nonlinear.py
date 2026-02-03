#!/usr/bin/env python3
import sys
sys.path.append("../")
from case import *

if __name__ == "__main__":
    root = os.getcwd()

    try:
        with open(file_ft, 'r') as file:
            target_total_time = float(file.read().strip())
    except FileNotFoundError:
        print(f"Not found file: {file_ft}")
    except Exception as e:
        print(f"An error occurred while reading the file: {str(e)}")

    print('')
    print('Target final time is set to:',target_total_time)

    pf = glob.glob('*.par')[0]
    cn = pf.split('/')[-1].split('.')[0]
    print(f"Case name: {cn}")
    print(f"Parameter file: {pf}")

    final_time_his = check_time(hisfile,which_time='final')
    initial_time_his = check_time(hisfile,which_time='initial')
    time_velocity_file = extract_time_from_binary(final_dns_file)

    print(f"Initial time in hisfile: {initial_time_his}")
    print(f"Final time in hisfile: {final_time_his}")
    print(f"Time in velocity file: {time_velocity_file}")

    if time_velocity_file < target_total_time:
        # we need to restart the simulation

        if initial_time_his > 0.0 and final_time_his == time_velocity_file:
            # hisfile is updated
            time_string = '_'+str(initial_time_his)+'_'+str(final_time_his)
            print(f"Making string: {time_string}")
            files = [final_dns_file, hisfile, liftdrag, globenergy, globenstro]
            for file in files:
                if os.path.exists(file):
                    print(f"Copying {file} to {file + time_string}")
                    shutil.copy(file, file + time_string)
                    if file != hisfile and file != final_dns_file:
                        print(f"Removing {file}")
                        os.remove(file)
                    else:
                        print(f"Resetting {file}")
                        reset_his_file(hisfile)
                else:
                    print(f"File {file} does not exist.")

        # job_name = f"{cns}{Re}d"
        # check_job_status(job_name)
        nexttime = time_velocity_file + delta_simu
        if time_velocity_file < delta_simu:
            nexttime = delta_simu

        c_pf(pf,pf, {'GENERAL':{'startfrom':str(final_dns_file)}})
        c_pf(pf,pf, {'GENERAL':{'endTime':str(nexttime)}})
        c_pf(pf,pf, {'GENERAL':{'writeinterval':str(delta_simu)}})

        # resubmit_job(pbs_file,folder_path,job_name)

        #print(f"Submitting another restart.")
        subprocess.call(['bash', 'run.sh'])

    else:    
        print(f"Current time: {final_time_his} == Final time: {target_total_time}")
        print(f"Current time is less than final time. No action taken.")
        sys.exit()

#append_files('1cyl.his_*_*', 'output.txt')