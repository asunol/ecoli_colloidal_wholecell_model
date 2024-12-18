###### Alp Sunol
# October 2022

# This calculates the MSD of GEMS.

import sys
import subprocess
import time
import os
import math
import numpy as np
import linecache
import copy


def get_value_from_line(name, line_num):  # warning, gets latest value that is a float
    line_old = linecache.getline(name, line_num)
    line = line_old.split(" ")
    #print(line)
    
    val = 0
    for strg in line:
        try:
            val = float(strg)
            
        except ValueError:
            pass
    
    return val 

sim_folder = sys.argv[1] # name of folder that contains output of simulation
sim_name = sys.argv[2] # base name for output file

print(sim_folder)
print(sim_name)

write_name = 'MSD.gems_only.'+sim_name+'.txt'

dump_fldr = sim_folder + 'dump/' # trajectories are in dump folder

R_nuc = 20.0 # each unit in simulations is 7 nm
R_cell = 27.4

V_nuc = 1.333333333*np.pi*R_nuc**3
V_cell = 1.333333333*np.pi*R_cell**3
V_cyt = V_cell - V_nuc
a_part = 1.42857 # radius of particle in simulation units, can change this if overall volume fraction matters (relative volume fraction won't change)
V_part = 1.333333333*np.pi*(a_part)**3



skip_head = 9

file1 = open(write_name, 'w')
file1.close()
print('start')

dump_file = dump_fldr+ 'cell-gems.1.2.lammpstrj' # we output just bGEM positions in lammps as a seperate file (in addition to the file containing all particles)
N_gem = int(get_value_from_line(dump_file, 4))
print('N gems is ', N_gem)


# Simulation time variables - these can change depending on the simulation settings
step_size = int(10000) # write down every this many steps
file_size = int(100000) # how many steps are in a file
start_num = int(1) # start analyzing at first timestep, can throw out when averaging
end_num = int(450) # last file number we want to look at
N_files = int(end_num - start_num) # number of files
steps_per_file = int(file_size/step_size)
print(steps_per_file)

# import all the data:
positions_GEM = np.zeros((N_files*steps_per_file*N_gem,4))
print(positions_GEM.shape)
print('Importing positions')


for num in range(start_num,end_num):
    num_str = str(num)
    #print(num)

    dump_file = dump_fldr+ 'cell-gems.1.' + num_str + '.lammpstrj'

    for step in range(0,steps_per_file):
        tstep = file_size*(num-start_num) + step*step_size
        skips = skip_head + (skip_head+N_gem)*step         
        pos_part =  np.genfromtxt(dump_file, skip_header = skips, max_rows = N_gem)
        #print(pos_part.shape)

        # LAMMPS does not write particle numbers in order, we have to re-order:
        pos_part = pos_part[pos_part[:, 0].argsort()]


        skip_save = (num-start_num)*steps_per_file*N_gem + step*N_gem
        positions_GEM[skip_save:skip_save+N_gem,0] = tstep
        positions_GEM[skip_save:skip_save+N_gem, 1: ] = pos_part[:,2:5]

#print(positions_GEM)
np.save(dump_fldr + '/gems.npy', positions_GEM) # save all bGEM locations in all times as a single file

positions_GEM = np.load(dump_fldr + '/gems.npy')

# Now calculate MSD
tsteps = int(positions_GEM.shape[0]/N_gem)
    
t_lag = np.zeros(tsteps)
MSD = np.zeros(tsteps)
    
# Open and close to reset.  If you want to append to an existing msd, comment this out.
file1 = open(write_name, "w")
file1.close()
print('Positions imported, calculating MSD')   

# Perform mean-squared displacement calculation
for i in range(tsteps-1):
    skip_idx = i+1
    t_lag[i] = step_size*skip_idx
    #print(t_lag[i])
            
    tmax = tsteps - skip_idx
    if tmax > 10000:
        tmax = 10000
        
        
    j_count = 0
    for j in range(0, tmax):
        skip1 = j*N_gem
        skip2 = N_gem*(j+skip_idx)
            
        pos_1 = positions_GEM[skip1:skip1+N_gem,1:]
        pos_2 = positions_GEM[skip2:skip2+N_gem,1:]

        disp = pos_2 - pos_1
            
        disp_squared = np.linalg.norm(disp,axis=1)**2.0
            
        MSD[i] += np.mean(disp_squared)
        j_count += 1
            
    MSD[i]/= j_count
    list_write = str(t_lag[i])+ ' '+ str(MSD[i]) + '\n'
        
    with open(write_name, "a") as file1:
        file1.write(list_write) # write it down as we go in case runtime stops early


## Now calculate MSD inside vs outside seperately
print('Calculating MSD inside vs outside')
pos_gem = positions_GEM
r_gem = np.linalg.norm(pos_gem[:,1:], axis = 1)
in_gem = r_gem>R_nuc
true_gem = in_gem + 0 # vector stating inside or outside the nucleoid

# Append this info to positions of GEMS
pos_gem_big = np.zeros((pos_gem.shape[0],  pos_gem.shape[1]+2))
pos_gem_big[:,0:4] = pos_gem
pos_gem_big[:,4] = r_gem
pos_gem_big[:,5] = true_gem

# New array, this will be in the order of particle by particle, not time by time:
print('Seperating trajectories in vs out')
pos_gem_big_new = np.zeros(pos_gem_big.shape)
for i in range(N_gem):
    pos_gem_i = pos_gem_big[i::N_gem,:]
    tsteps = pos_gem_i.shape[0]
    pos_gem_big_new[i*tsteps:(i+1)*tsteps,:] = pos_gem_i


# Build a list of trajectories that are inside and outside.
# These trajectories can vary in length.  
# Each trajectory ends if bGEM switches from in/out to out/in.
traj_list_inside = list()
traj_list_outside = list()
for gem_id in range(N_gem):
    # split the trajectory if it goes from inside to outside
    gem_jump = np.diff(pos_gem_big_new[(gem_id)*tsteps:(gem_id+1)*(tsteps),5])
    split_id = np.where(np.abs(gem_jump) > 0)[0]
    split_traj = np.split(pos_gem_big_new[(gem_id)*tsteps:(gem_id+1)*(tsteps),:], split_id+1, axis = 0)
    for traj in split_traj:
        traj_new = np.zeros((traj.shape[0],  traj.shape[1]+1))
        traj_new[:,0] = gem_id + 1
        traj_new[:,1:] = traj
        #print(traj_new)
        if traj_new[0,-1] == 0:  # append to inside trajcetory
            traj_list_inside.append(traj_new)
        if traj_new[0,-1] == 1: # append to outside trajectory
            traj_list_outside.append(traj_new)


# Go through each trajectory and do displacements as a function of time.

# First, trajectories in nucleoid:
msd_in = np.zeros(tsteps-1)
msd_in_cnt = np.zeros(tsteps-1)
print('Calculating inside MSD')
for traj in traj_list_inside:
    bout = traj.shape[0]-1

    
    for i in range(bout):
        skip_idx = i

        for j in range(bout-skip_idx):
            pos_1 = traj[j,2:5]
            pos_2 = traj[j+skip_idx,2:5]
            disp = pos_2 - pos_1
            disp_squared = np.linalg.norm(disp)**2.0
            msd_in[skip_idx] += disp_squared
            msd_in_cnt[skip_idx] += 1           
MSD_in = msd_in/msd_in_cnt
write_name_in = 'MSD.gems_only.'+sim_name+'.in.txt'
print(write_name_in)
np.savetxt(write_name_in, MSD_in)
print('Inside MSD done')

# Now, trajectories outside nucleoid:
print('Calculating outside MSD')
msd_out = np.zeros(tsteps-1)
msd_out_cnt = np.zeros(tsteps-1)
for traj in traj_list_outside:
    bout = traj.shape[0]-1
    #hist_traj_inside[bout] += 1
    #print('bout ', bout)
    for i in range(bout):
        skip_idx = i
        #print('skip', skip_idx)      
        for j in range(bout-skip_idx):
            #print(j+skip_idx, 'minus', j)
            pos_1 = traj[j,2:5]
            pos_2 = traj[j+skip_idx,2:5]
            disp = pos_2 - pos_1
            disp_squared = np.linalg.norm(disp)**2.0
            msd_out[skip_idx] += disp_squared
            msd_out_cnt[skip_idx] += 1           
MSD_out = msd_out/msd_out_cnt
write_name_out = 'MSD.gems_only.'+sim_name+'.out.txt'
np.savetxt(write_name_out, MSD_out)
print('Outside MSD done')



## Save Ribosome locations (extra, for other postprocessing if needed):
dump_rib = dump_fldr+ 'cell-ribs.1.2.lammpstrj'
N_rib = int(get_value_from_line(dump_rib, 4))
positions_rib = np.zeros((N_files*steps_per_file*N_rib,4))
print('Saving ribosomes.  N ribs is ', N_rib)
for num in range(start_num,end_num):
    num_str = str(num)

    dump_file = dump_fldr+ 'cell-ribs.1.' + num_str + '.lammpstrj'

    for step in range(0,steps_per_file):
        tstep = file_size*(num-start_num) + step*step_size
        skips = skip_head + (skip_head+N_rib)*step         
        pos_part =  np.genfromtxt(dump_file, skip_header = skips, max_rows = N_rib)

        # LAMMPS does not write particle numbers in order, we have to re-order:
        pos_part = pos_part[pos_part[:, 0].argsort()]


        skip_save = (num-start_num)*steps_per_file*N_rib + step*N_rib
        positions_rib[skip_save:skip_save+N_rib,0] = tstep
        positions_rib[skip_save:skip_save+N_rib, 1: ] = pos_part[:,2:5]

np.save(dump_fldr + '/ribs.npy', positions_rib)
