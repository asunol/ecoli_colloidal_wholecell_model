# Alp Sunol
# This calculates the contact distribution around GEMS

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
write_name = sys.argv[2] # base write name for output file


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

dump_file = dump_fldr+ 'cell-pd.1.2.lammpstrj' # take a sample file to see how many particles there are
N_part = int(get_value_from_line(dump_file, 4))
print('N part is ', N_part)

# Simulation time variables - these can change depending on the simulation settings
step_size = int(10000) # write down every this many steps
file_size = int(100000) # how many steps are in a file
start_num = int(1) # start analyzing at first timestep, can throw out when averaging
end_num = int(450) # last file number we want to look at
N_files = int(end_num - start_num) # number of files
steps_per_file = int(file_size/step_size)
print(steps_per_file)

print('Importing positions')

count = 1
for num in range(start_num,end_num):
    num_str = str(num)
    print(num)

    dump_file = dump_fldr+ 'cell-pd.1.' + num_str + '.lammpstrj' # file name

    for step in range(0,steps_per_file):
        tstep = file_size*(num-start_num) + step*step_size
        skips = skip_head + (skip_head+N_part)*step         
        pos_part =  np.genfromtxt(dump_file, skip_header = skips, max_rows = N_part) # import lammps trajectory file
        #print(pos_part.shape)

        # LAMMPS does not write particle numbers in order, we have to re-order:
        pos_part = pos_part[pos_part[:, 0].argsort()]
        
        save_name = '/particles_'+str(count)+'.npy' # optional, save as a numpy array for quicker access later

        np.save(dump_fldr + save_name, pos_part[:,1:5])
        count += 1


diam_gem = 40  # nm, need to change this based off bGEM size
part_sizes = np.array([0.9, 0.95, 1.0, 1.05, 1.1, diam_gem/14, 26/14, 7/7, 7/7]) # size of all particles in simulation units
cutoff = 1/14  # 1 nm in simulation units
R_nuc = 20.0
N_GEM = 50

snapshot_start = 1
snapshot_end = 2500

bin_edges = np.arange(-.5, 21.5, 1)

# Open files in write mode to clear content before starting the loop
write_name_DNA = write_name + '_DNA_all.txt'
write_name_rib = write_name + '_ribosome.txt'
write_name_negcrowd = write_name + '_negcrowd.txt'
write_name_poscrowd = write_name + '_poscrowd.txt'

with open(write_name_DNA, 'w'), open(write_name_rib, 'w'), open(write_name_negcrowd, 'w'), open(write_name_poscrowd, 'w'):
    pass  # This just opens the files and clears their contents

# Loop through and see how many of each particle type is within the distance cutoff of a bGEM
for i in range(snapshot_start, snapshot_end):
    if i%100==0:
        print(i)
    particle_file = dump_fldr + '/particles_' + str(i) + '.npy'
    pos_part = np.load(particle_file)
    pos_GEM = pos_part[pos_part[:, 0] == 6, :]
    pos_DNA = pos_part[pos_part[:, 0] <= 5, :]
    pos_rib = pos_part[pos_part[:, 0] == 7, :]
    pos_negcrowd = pos_part[pos_part[:, 0] == 8, :]
    pos_poscrowd = pos_part[pos_part[:, 0] == 9, :]
    
    
    # GEM-DNA distances
    gem_DNA_disp = pos_GEM[:, np.newaxis, 1:] - pos_DNA[np.newaxis, :, 1:]
    gem_DNA_distances = np.linalg.norm(gem_DNA_disp, axis=2)
    
    # Calculate thresholds for each DNA particle
    threshold_DNA = cutoff + part_sizes[5]/2 + part_sizes[(pos_DNA[:, 0].astype(int) - 1)]/2
    
    within_threshold_DNA = gem_DNA_distances < threshold_DNA[np.newaxis, :]
    
    # Count and bin for DNA
    num_neighbors_DNA = np.sum(within_threshold_DNA, axis=1)
    hist_DNA, _ = np.histogram(num_neighbors_DNA, bins=bin_edges)
    hist_DNA = hist_DNA / N_GEM
    
    with open(write_name_DNA, 'a') as file_DNA:
        file_DNA.write(' '.join(map(str, hist_DNA)) + '\n')
        
        
    # GEM-Ribosome distances
    gem_rib_disp = pos_GEM[:, np.newaxis, 1:] - pos_rib[np.newaxis, :, 1:]
    gem_rib_distances = np.linalg.norm(gem_rib_disp, axis=2)
    
    threshold_rib = cutoff + part_sizes[5]/2 + part_sizes[6]/2
    
    within_threshold_rib = gem_rib_distances < threshold_rib
    
    num_neighbors_rib = np.sum(within_threshold_rib, axis=1)
    hist_rib, _ = np.histogram(num_neighbors_rib, bins=bin_edges)
    hist_rib = hist_rib / N_GEM
    
    with open(write_name_rib, 'a') as file_rib:
        file_rib.write(' '.join(map(str, hist_rib)) + '\n')
    
    # GEM-Negative Crowder distances
    gem_negcrowd_disp = pos_GEM[:, np.newaxis, 1:] - pos_negcrowd[np.newaxis, :, 1:]
    gem_negcrowd_distances = np.linalg.norm(gem_negcrowd_disp, axis=2)
    
    threshold_negcrowd = cutoff + part_sizes[5]/2 + part_sizes[7]/2
    
    within_threshold_negcrowd = gem_negcrowd_distances < threshold_negcrowd
    
    num_neighbors_negcrowd = np.sum(within_threshold_negcrowd, axis=1)
    hist_negcrowd, _ = np.histogram(num_neighbors_negcrowd, bins=bin_edges)
    hist_negcrowd = hist_negcrowd / N_GEM
    
    with open(write_name_negcrowd, 'a') as file_negcrowd:
        file_negcrowd.write(' '.join(map(str, hist_negcrowd)) + '\n')
    
    # GEM-Positive Crowder distances
    gem_poscrowd_disp = pos_GEM[:, np.newaxis, 1:] - pos_poscrowd[np.newaxis, :, 1:]
    gem_poscrowd_distances = np.linalg.norm(gem_poscrowd_disp, axis=2)
    
    threshold_poscrowd = cutoff + part_sizes[5]/2 + part_sizes[8]/2
    
    within_threshold_poscrowd = gem_poscrowd_distances < threshold_poscrowd
    
    num_neighbors_poscrowd = np.sum(within_threshold_poscrowd, axis=1)
    hist_poscrowd, _ = np.histogram(num_neighbors_poscrowd, bins=bin_edges)
    hist_poscrowd = hist_poscrowd / N_GEM
    
    with open(write_name_poscrowd, 'a') as file_poscrowd:
        file_poscrowd.write(' '.join(map(str, hist_poscrowd)) + '\n')
