###### Alp Sunol
# Janruary 2022

# Calculates volume fraction of each species inside and out of the nucleoid.

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
sim_name = sys.argv[2] # base name for output files
n_types = int(sys.argv[3])


write_name = 'phi.'+sim_name+'.' # name stem of output files

dump_fldr = sim_folder + 'dump/'




R_nuc = 20.0 # simulation units, where 1 unit is 14 nm
R_cell = 27.4

V_nuc = 1.333333333*np.pi*R_nuc**3
V_cell = 1.333333333*np.pi*R_cell**3
V_cyt = V_cell - V_nuc # cytoplasm outside the nucleoid

gem_radius = float(sys.argv[4])/28.0  # input is the diameter in nanometers

rib_radius = 13.0/14.0
a_part = np.array([gem_radius, rib_radius, 0.5, 0.5]) # volume of particle types we want to track, which are bGEM, ribosome, and negative/positive crowders

V_part = 1.333333333*np.pi*(a_part)**3

skip_head = 9

names_list = ['GEM', 'ribosome', 'neg_prot', 'pos_prot']

# Create files
for write_name_specific in names_list:
    file1 = open(write_name + write_name_specific + '.txt', 'w')
    file1.close()


print('start')

dump_file = dump_fldr+ 'cell-pd.1.2.lammpstrj' # can use any of the lammpstrj files, this is just to get number of particles
N_tot = get_value_from_line(dump_file, 4) # number of particles


# Loop through each file.  
# We only look at first output in each file to ensure locations of bGEMs are not very highly correlated each time we record.
for num in range(1,500, 1):
    num_str = str(num)
    print(num)

    dump_file = dump_fldr+ 'cell-pd.1.' + num_str + '.lammpstrj'

    pos_part =  np.genfromtxt(dump_file, skip_header = skip_head, max_rows = N_tot)

    N_diff = np.zeros(n_types)
    N_nuc = np.zeros(n_types)
    N_cyt = np.zeros(n_types)
    for i in range(n_types):
        # First 5 types of particles are DNA which are always in nucleoid.
        pos_diff = pos_part[pos_part[:,1] == 6+i] # diffusers
        N_diff[i] = pos_diff.shape[0]

        r_diff = np.linalg.norm(pos_diff[:,2:5], axis = 1)

        N_nuc[i] = np.sum(r_diff<R_nuc) # diffusers in nucleoid
        N_cyt[i] = N_diff[i] - N_nuc[i]
    
    phi_cyt =  (N_cyt)*V_part/V_cyt
    phi_nuc =  (N_nuc)*V_part/V_nuc
    phi_tot = N_diff*V_part/V_cell

    # For each particle type, write the volume fraction in, out, and total
    for i in range(n_types):
        strng_write = num_str + ' ' + str(phi_nuc[i]) + ' ' +  str(phi_cyt[i]) + ' ' +  str(phi_tot[i])
        strng_write += '\n'
        name_of_file = write_name + names_list[i] + '.txt'
        #print(name_of_file)
        #print(strng_write)
        with open(name_of_file, "a") as file_write:
            file_write.write(strng_write)
    