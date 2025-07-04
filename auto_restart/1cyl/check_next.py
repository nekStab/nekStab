#!/usr/bin/env python3
import sys
sys.path.append("../")
from case import *

if __name__ == "__main__":
    root = os.getcwd()
    re = directory_name = os.path.basename(root)
    try:
        with open(file_ft, 'r') as file:
            lines = file.readlines()
            target_total_time = float(lines[0].strip())
            delta_simu = float(lines[1].strip())
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
    if np.isnan(time_velocity_file): 
        print(f"The file {final_dns_file} does not exist.")
        print('Stopping file.')
        sys.exit()
        
    print(f"Initial time in hisfile: {initial_time_his}")
    print(f"Final time in hisfile: {final_time_his}")
    print(f"Time in velocity file: {time_velocity_file}")
    
    if initial_time_his >= 0.0: #and final_time_his == time_velocity_file:
        # hisfile is updated
        if (initial_time_his == 0 and final_time_his == 0):
            time_string = '_'+str(initial_time_his)
        else:   
            time_string = '_'+str(initial_time_his)+'_'+str(final_time_his)
        
        print(f"Making string: {time_string}")
        files = [final_dns_file, hisfile, liftdrag, globenergy, globenstro, 'logfile', 'logerror']
        for file in files:
            if os.path.exists(file):
                print(f"Copying {file} to {file + time_string}")
                shutil.copy(file, file + time_string)
                if file != hisfile and file != final_dns_file:
                    print(f"Removing {file}")
                    os.remove(file)
                elif file == hisfile:
                    print(f"Resetting {file}")
                    reset_his_file(file)
            else:
                print(f"File {file} does not exist.")

    if time_velocity_file < target_total_time:
        # we need to restart the simulation

        nexttime = time_velocity_file + delta_simu
        if time_velocity_file < delta_simu:
            nexttime = delta_simu

        c_pf(pf,pf, {'GENERAL':{'startfrom':str(final_dns_file)}})
        c_pf(pf,pf, {'GENERAL':{'endTime':str(nexttime)}})
        c_pf(pf,pf, {'GENERAL':{'writeinterval':str(delta_simu)}})
        c_pf(pf,pf, {'VELOCITY':{'viscosity':f"-{float(re)}"}})
        
        job_name = f"{cns}{re}_{nexttime}"
        job_name = f"{job_name[:8]:<8}"  # Truncate to 8 

        if not check_job_exist(job_name):# job is up!
            resubmit_job(pbs_file,root,job_name)

        #print(f"Submitting another restart.")
        #subprocess.call(['bash', 'run.sh'])

    else:    
        print(f"Current time: {final_time_his} == Final time: {target_total_time}")
        print(f"Current time is less than final time. No action taken.")

    append_files('lift_drag.dat_*', 'lift_drag.all')
    append_files('total_energy.dat_*', 'total_energy.all')
    append_files('total_enstrophy.dat_*', 'total_enstrophy.all')
    append_files(f'{cn}.his_*_*', f'{cn}.all')
