# Alp M. Sunol
# Started Janruary 2023
# Last edited August 2024

# This generates a random walk inside E. coli.
# Can explicitely define seperate diffusion coefficients inside and outside the nucleoid.
# Can define entry/exist probablilities from the nucleoid.

import sys
import subprocess
import time
import os
import math
import numpy as np
import linecache
import random
import copy



### Functions

"""
Function takes cartesian coordinates pos = Nx3([x,y,z]) as well as cell and nucleoid size parameters 
and determines location inside cell.

input: pos - x coming out of page, z going right, y going up; matrix of particles
       nucleoid and cell size properties
output: vector of all particle with the numbers for each particle denoting:
        0 - inside nucleoid
        1 - outside nucleoid, inside cell
        2 - outside cell
""" 
def cell_region(pos, l_cell, r_cell, r_cap, l_nuc, r_nuc, R_part, gap_min, gap_max): 
    x = pos[:, 0]
    y = pos[:, 1]
    z = pos[:, 2]
    
    # First, analyze particle in the caps:
    cyl_or_cap = np.abs(z) > l_cell  # true if potentially in cap region
    cap_idx = np.where(cyl_or_cap)[0]
    pos_cap = pos[cap_idx, :]  # position vector of positions in cap
    pos_cap[:, 2] -= np.sign(pos_cap[:, 2]) * l_cell  # re-center such that it's a sphere with center where cap starts
    r_coord_cap = np.linalg.norm(pos_cap, axis=1)  # radius within cap coordinates

    cap_out = r_coord_cap > r_cap - R_part  # true if outside the cell
    loc_cap = cap_out + 1  # 2 if outside cell, 1 if inside cell 

    # Now, analyze the ones in the main cylinder
    cap_or_cyl = np.abs(z) < l_cell  # true if not in cap region
    cyl_idx = np.where(cap_or_cyl)[0]
    pos_cyl = pos[cyl_idx, :]
    # Convert to relevant cylindrical coordinates
    r_coord_cyl = np.linalg.norm(pos_cyl[:, :2], axis=1)
    z_coord_cyl = pos_cyl[:, 2]

    # Define the gap region where particles are not considered inside the nucleoid
    gap_region = (z_coord_cyl > gap_min) & (z_coord_cyl < gap_max)

    # Assign outside the cell:
    out_cyl = np.where(r_coord_cyl > r_cell - R_part)[0]
    loc_cyl = np.empty(r_coord_cyl.shape, dtype=np.int64)
    loc_cyl[out_cyl] = int(2)

    # Assign inside cell but outside nucleoid
    in_cyt1 = np.where((r_nuc < r_coord_cyl) & (r_coord_cyl < r_cell - R_part))
    loc_cyl[in_cyt1] = int(1)
    in_cyt2 = np.where((r_nuc > r_coord_cyl) & ((np.abs(z_coord_cyl) > l_nuc) | gap_region))
    loc_cyl[in_cyt2] = int(1)

    # Assign inside nucleoid (excluding gap region):
    in_nuc = np.where((r_nuc > r_coord_cyl) & (np.abs(z_coord_cyl) < l_nuc) & ~gap_region)
    loc_cyl[in_nuc] = int(0)

    # Combine cap and cylinder data
    loc = np.empty(x.shape)
    loc[cap_idx] = loc_cap
    loc[cyl_idx] = loc_cyl
    
    return loc


"""
Function initializes random locations of N particles throughout the cell. Ratio inside vs outside the nucleoid
can be set. The distributions inside and outside the cell are uniformly distributed.

input: N - number of particles
       pen - density inside divided by outside (pen = N_in/V_in / N_out/V_out)
       nucleoid and cell size properties
       gap_min - minimum z-coordinate of the gap
       gap_max - maximum z-coordinate of the gap
output: position vector in cartesian coordinates that are distributed as specified with the input
"""    
def init_distribution(N, pen, l_cell, r_cell, r_cap, l_nuc, r_nuc, R_part, gap_min, gap_max): 
    
    
    
    # Determine number in each region
    V_caps = (4/3) * np.pi * (r_cap**3)  # there are two half caps
    V_cyl = 2 * l_cell * np.pi * (r_cell**2)
    V_cell = V_caps + V_cyl
    
    # Subtract the volume of the gap from the nucleoid volume
    V_gap = 2 * (gap_max - gap_min) * np.pi * (r_nuc**2)
    V_nuc = 2 * l_nuc * np.pi * (r_nuc**2) - V_gap

    V_in = V_nuc
    V_out = V_cell - V_nuc

    if pen == 0:
        N_in = 0
    else:
        N_in = int(N * (1 + V_out/V_in/pen)**(-1))
    if pen > 1e4:
        N_in = N
    N_out = N - N_in
    print('inside:', N_in, 'outside:', N_out)
      
    
    
    # Determine number in each region outside the nucleoid
    # Region 1: in the caps 
    V1 = (4/3) * np.pi * (r_cap**3)  # Volume of the caps
    
    # Region 2: outside the nucleoid in z direction but not in caps
    V2 = 2 * (l_cell - l_nuc) * np.pi * (r_cell**2)  # Volume of region 2

    # Region 3: around the nucleoid ends in r direction + gap region
    # Calculate V3, which now includes the gap (Va) and the region between nucleoid and cell wall (Vb)
    Va = (gap_max - gap_min) * np.pi * (r_cell**2)  # Volume of the gap region (cylinder)

    # Volume of the regions between the nucleoid and cell wall but outside the gap
    Vb =  ((l_nuc - gap_max) + (gap_min+l_nuc) ) * np.pi * (r_cell**2 - r_nuc**2) 
    
    V3 = 0  # Initialize V3
    
    # Add the gap volume to V3 if the particle fits in the gap
    if 2 * R_part <= (gap_max - gap_min):
        print('Particles fit in gap region')
        V3 += Va

    # Add the volume outside the nucleoid to V3 if the particle fits between the nucleoid and cell wall
    if r_cell - 2*R_part>r_nuc:
        print('Particles fit around nucleoid')
        V3 += Vb

    # Determine number of particles in each region based on available volume
    V_total = V1 + V2 + V3
    N_out_1 = int(N_out * V1 / V_total)  # Number of particles in the caps
    if V3 == 0:
        N_out_2 = N_out - N_out_1
        N_out_3 = 0
    else:
        N_out_2 = int(N_out * V2 / V_total)  # Number of particles in region 2
        N_out_3 = N_out - N_out_1 - N_out_2  # Number of particles in region 3
        


    # Further split N_out_3 into particles that go into the gap (Na) and those that go into the nucleoid-cell wall region (Nb)
    Na, Nb = 0, 0
    mode = 0 # none in either
    if V3 > 0:
        if 2 * R_part <= (gap_max - gap_min):
            if r_cell - 2 * R_part < r_nuc:
                Na = N_out_3
            else:
                Na = int(N_out_3 * Va / V3)
        if r_cell - 2 * R_part > r_nuc:
            Nb = N_out_3 - Na
    #print(f"V_cell (entire cell): {V_cell}")
    #print(f"V_in (nucleoid): {V_nuc}")
    #print(f"V_out (cytoplasm): {V_out}")
    #print(f"V_out_avail (cytoplasmic regions GEM can fit into): {V_total}")
    #print(f"V ratio): {V_in/V_out") 
    #print(f"V1 (caps): {V1}")
    #print(f"V2 (region 2): {V2}")
    #print(f"V3 (region 3): {V3}")
    print(f"N_in (inside nucleoid): {N_in}")
    print(f"N_out (outside nucleoid): {N_out}")
    print(f"N_out_1 (caps): {N_out_1}")
    print(f"N_out_2 (region 2): {N_out_2}")
    print(f"N_out_3 (region 3): {N_out_3}")
    print(f"N3a (gap): {Na}")
    print(f"N3b (nucleoid-cell wall region): {Nb}")

    # Now, we place the particles.
    # First, randomly distribute the ones inside the Nucleoid:
    z_in_initial = np.random.uniform(-l_nuc, l_nuc, (N_in, 1))

    # Transform z_in_initial to avoid the gap
    z_in = np.zeros_like(z_in_initial)

    # For particles below the gap (lower_mask)
    lower_mask = z_in_initial < 0
    z_in[lower_mask] = z_in_initial[lower_mask] * (l_nuc + gap_min) / l_nuc + gap_min

    # For particles above the gap (upper_mask)
    upper_mask = z_in_initial >= 0
    z_in[upper_mask] = z_in_initial[upper_mask] * (l_nuc - gap_max) / l_nuc + gap_max

    r_in = r_nuc * np.sqrt(np.random.uniform(0, 1, (N_in, 1)))
    theta_in = np.random.uniform(0, 2*np.pi, (N_in, 1))
    x_in = r_in * np.sin(theta_in)
    y_in = r_in * np.cos(theta_in)

    pos_in = np.hstack((x_in, y_in, z_in))

    # Region 1: in the caps 
    r_out_1 = (r_cap - R_part) * (np.random.uniform(0, 1, (N_out_1, 1)))**(1/3)
    theta_out_1 = 2 * np.pi * np.random.uniform(0, 1, (N_out_1, 1))
    phi_out_1 = np.arccos(np.random.uniform(-1, 1, (N_out_1, 1)))
    x_out_1 = r_out_1 * np.sin(phi_out_1) * np.cos(theta_out_1)
    y_out_1 = r_out_1 * np.sin(phi_out_1) * np.sin(theta_out_1)
    z_out_1 = r_out_1 * np.cos(phi_out_1)
    # Now split the sphere into the caps:
    z_out_1 += np.sign(z_out_1) * l_cell 
    pos_out_1 = np.hstack((x_out_1, y_out_1, z_out_1))

    # Region 2: outside the nucleoid in z direction but not in caps
    pos_neg_rand = 2 * np.random.randint(0, 2, size=(N_out_2, 1)) - 1
    z_out_2 = pos_neg_rand * np.random.uniform(l_nuc, l_cell, (N_out_2, 1))
    r_out_2 = (r_cell - R_part) * np.sqrt(np.random.uniform(0, 1, (N_out_2, 1)))
    theta_out_2 = np.random.uniform(0, 2 * np.pi, (N_out_2, 1))
    x_out_2 = r_out_2 * np.sin(theta_out_2)
    y_out_2 = r_out_2 * np.cos(theta_out_2)
    pos_out_2 = np.hstack((x_out_2, y_out_2, z_out_2))

    # Region 3: around the nucleoid but higher in the z direction and gap
    pos_out_3a, pos_out_3b = None, None
    if Na > 0:  # Place particles in the gap
        z_out_3a = np.random.uniform(gap_min, gap_max, (Na, 1))
        r_out_3a = (r_cell - R_part) * np.sqrt(np.random.uniform(0, 1, (Na, 1)))
        theta_out_3a = np.random.uniform(0, 2 * np.pi, (Na, 1))
        x_out_3a = r_out_3a * np.sin(theta_out_3a)
        y_out_3a = r_out_3a * np.cos(theta_out_3a)
        pos_out_3a = np.hstack((x_out_3a, y_out_3a, z_out_3a))
    if Nb > 0:  # Place particles in the nucleoid-cell wall region
        z_out_3b = np.random.uniform(-l_nuc, l_nuc, (Nb, 1))
        r_out_3b = np.random.uniform(r_nuc + R_part, r_cell - R_part, (Nb, 1))
        theta_out_3b = np.random.uniform(0, 2 * np.pi, (Nb, 1))
        x_out_3b = r_out_3b * np.sin(theta_out_3b)
        y_out_3b = r_out_3b * np.cos(theta_out_3b)
        pos_out_3b = np.hstack((x_out_3b, y_out_3b, z_out_3b))
    # Combine region 3 positions
    if Na > 0 and Nb > 0:
        pos_out_3 = np.vstack((pos_out_3a, pos_out_3b))
    elif Na > 0:
        pos_out_3 = pos_out_3a
    elif Nb > 0:
        pos_out_3 = pos_out_3b
    else:
        pos_out_3 = np.empty((0, 3))  # Ensure pos_out_3 is defined as an empty array if neither is used


    # Combine all positions
    pos_init = np.vstack((pos_in, pos_out_1, pos_out_2, pos_out_3))
    
    return pos_init



'''
Function generates a random walk throughout E. coli.

input: pos_init - matrix of initial positions of format Nx3(x,y,z)
       dt - timestep of simulation
       t_tot - total time of simulation
       D_in - diffusion coefficient inside the nucleoid
       D_out - diffusion coefficient outside the nucleoid (rest of the cytoplasm)
       pen_in - inside penetration probability (probability of going outside to in, 1 goes in, 0 stays out)
       pen_out - outside penetration probability (probability of going inside to out, 1 goes out, 0 stays in)
       inputs to cell_region (l_cell, r_cell, r_cap, l_nuc, r_nuc, R_part)
       
output: pos - position vector of shape (tsteps, N_part, dim)
'''
def ecoli_walk(pos_init, dt, t_tot, D_in, D_out, pen_in, pen_out, l_cell, r_cell, r_cap, l_nuc, r_nuc, R_part, gap_min, gap_max):
    
    tsteps = int(t_tot / dt)
    N_part = pos_init.shape[0]
    tot_out = np.zeros((tsteps + 1,))  # how many are in which region
    dim = 3
    print('Initialize displacements')
    rand_disp = np.random.normal(0, 1, (tsteps, N, dim))  # create all random displacements
    pos = np.zeros((tsteps + 1, N_part, dim))
    t = np.zeros(tsteps + 1)
    
    pos[0, :, :] = pos_init  # initialize
    tot_out[0] = np.sum(cell_region(pos_init, l_cell, r_cell, r_cap, l_nuc, r_nuc, R_part, gap_min, gap_max))

    print('Evolving time forward')
    for idx, disp in enumerate(rand_disp):
        percentage_done = 100.0 * idx / tsteps
        # print(idx, percentage_done, tsteps)
        if np.isclose((100.0 * idx / tsteps) % 10.0, 0, atol=100.0 / tsteps):
            print(100.0 * idx / tsteps, '% done')
            
        pos_old = pos[idx, :, :]  # position at time t
        region_old = cell_region(pos_old, l_cell, r_cell, r_cap, l_nuc, r_nuc, R_part, gap_min, gap_max)  # where is it

        D_vec = np.zeros((N, 1))  # vector of diffusion coefficients for separate nucleoid and cytoplasm diffusivity
        D_vec[np.where(region_old == 0)] = D_in
        D_vec[np.where(region_old == 1)] = D_out

        pos_new = pos_old + disp * np.sqrt(2 * D_vec * dt)
        region_new = cell_region(pos_new, l_cell, r_cell, r_cap, l_nuc, r_nuc, R_part, gap_min, gap_max)

        change_region = region_new - region_old

        # Doesn't change region:
        idx_same = np.where(change_region == 0)
        pos[idx + 1, idx_same, :] = pos_new[idx_same, :]

        # Tries to go outside the cell:
        idx_outside_cell = np.where(region_new == 2)[0]
        pos[idx + 1, idx_outside_cell, :] = pos_old[idx_outside_cell, :]  # reject the movement, time still goes forward

        # Tries to go into the nucleoid
        idx_into_nuc = np.where((region_new == 0) & (change_region == -1))[0]
        pen_decision_in = np.random.uniform(size=idx_into_nuc.shape)
        conf_out_nuc = np.where((pen_decision_in > pen_in))  # stay confined outside nucleoid
        pen_in_nuc = np.where((pen_decision_in < pen_in))  # penetrate into the nucleoid
        pos[idx + 1, idx_into_nuc[conf_out_nuc], :] = pos_old[idx_into_nuc[conf_out_nuc], :]
        pos[idx + 1, idx_into_nuc[pen_in_nuc], :] = pos_new[idx_into_nuc[pen_in_nuc], :]

        # Tries to go outside the nucleoid
        idx_out_nuc = np.where((region_new == 1) & (change_region == 1))[0]
        pen_decision_out = np.random.uniform(size=idx_out_nuc.shape)
        conf_in_nuc = np.where((pen_decision_out > pen_out))  # stay confined inside nucleoid
        pen_out_nuc = np.where((pen_decision_out < pen_out))  # penetrate out of the nucleoid
        pos[idx + 1, idx_out_nuc[conf_in_nuc], :] = pos_old[idx_out_nuc[conf_in_nuc], :]
        pos[idx + 1, idx_out_nuc[pen_out_nuc], :] = pos_new[idx_out_nuc[pen_out_nuc], :]

        tot_out[idx + 1] = np.sum(cell_region(pos[idx + 1, :, :], l_cell, r_cell, r_cap, l_nuc, r_nuc, R_part, gap_min, gap_max))

        t[idx + 1] = t[idx] + dt
        
    return pos, t, tot_out



'''
Function calculates mean-squared displacement from timeseries position data.

input: pos - position vector of shape (tsteps, N_part, dim)
       dt - time between two steps (simulation timestep)
       lag_step - how many simulation timesteps we want a lag time output for 
       t_end_frac - fraction of the time of simulation that we want MSD calcs for
output: lag time vector t_lag and MSD vector
'''

def calc_msd(pos, dt, lag_step=1, t_end_frac=1.0):
    tsteps = int((t_end_frac * (pos.shape[0] - 1) / lag_step))  # number of timesteps of output
    print(tsteps)
    t_lag = np.zeros(tsteps)
    MSD = np.zeros(tsteps)

    # Calculate MSD of random walk:
    for i in range(tsteps):
        percentage_done = 100 * i / tsteps
        # print(percentage_done, i)
        if np.isclose((percentage_done) % 10.0, 0, atol=100 / tsteps):
            print(percentage_done, '% done')

        skip_idx = lag_step * (i + 1)
        t_lag[i] = dt * skip_idx
        # print(t_lag[i])

        tmax = tsteps * int(lag_step) - skip_idx

        if tmax > 10000:
            tmax = 10000

        j_count = 0
        for j in range(0, tmax + 1):
            skip1 = j
            skip2 = j + skip_idx
            # print(i, skip1, skip2)

            pos_1 = pos[skip1, :, :]
            pos_2 = pos[skip2, :, :]

            disp = pos_2 - pos_1

            disp_squared = np.linalg.norm(disp, axis=1)**2.0

            MSD[i] += np.mean(disp_squared)
            j_count += 1

        MSD[i] /= j_count
        # print(MSD[i])
    return t_lag, MSD

#########################################################################################
# Enter in nucleoid and cell size parameters.
# l_cell = 1100.0 # nm, cylindrical length of cell not including caps (divided by 2)
# r_cell = 400.0 # nm, cylindrical radius of cell
# r_cap = 400.0 # nm, radius pole caps

# l_nuc = 775.0 # nm, cylindrical length of nucleoid not including caps (divided by 2)
# r_nuc = 280.0 # nm, cylindrical radius of nucleoid

# Stationary phase:
l_cell = 700.0  # nm, cylindrical length of cell not including caps (divided by 2)
r_cell = 400.0  # nm, cylindrical radius of cell
r_cap = 400.0   # nm, radius pole caps

l_nuc = 650     # nm, cylindrical length of nucleoid not including caps (divided by 2)
r_nuc = 300.0   # nm, cylindrical radius of nucleoid
gap_max = 0     # nm, z value of where gap between nucleoid bilobe ends (a positive value), no gap in stationary phase
gap_min = -1*gap_max     # nm, z value of where gap between nucleoid bilobe starts (a negative value)

# Exponential phase:
#l_cell = 1000.0  # nm, cylindrical length of cell not including caps (divided by 2)
#r_cell = 350.0  # nm, cylindrical radius of cell
#r_cap = 350.0   # nm, radius pole caps

#l_nuc = 950     # nm, cylindrical length of nucleoid not including caps (divided by 2)
#r_nuc = 300.0   # nm, cylindrical radius of nucleoid
#r_nuc = 390.0   # nm, cylindrical radius of nucleoid if we want bGEMs to not be able to move from caps to the gap

# The below parameters are zero if the nucleoid is not bi-lobed
#gap_max = 200     # nm, z value of where gap between nucleoid bilobe ends (a positive value)
#gap_min = -1*gap_max     # nm, z value of where gap between nucleoid bilobe starts (a negative value)


V_cell = (4.0 / 3.0) * np.pi * (r_cap**3) + 2 * l_cell * np.pi * (r_cell**2)  # nm^3
V_nuc = 2 * l_nuc * np.pi * (r_nuc**2)  # nm^3

# GEM parameters:
#R_gem = 10.0  # nm, radius of GEM
R_gem = float(sys.argv[1])
print('gem radius = ', R_gem)

N = 500

# Define penetrability.  If 0, everything is outside.  If >1e4, everything is inside.
pen = 1.0  # overall penetrability of initial condition, also defines the ratios for dynamics
pen = float(sys.argv[2])

print('penetrability= ', pen)
pos_init = init_distribution(N, pen, l_cell, r_cell, r_cap, l_nuc, r_nuc, R_gem, gap_min, gap_max)

loc = cell_region(pos_init, l_cell, r_cell, r_cap, l_nuc, r_nuc, R_gem, gap_min, gap_max)
np.sum(loc)

pen_in = 0  # probability of penetrating into the nucleoid, will redefine later
pen_out = 0 # probability of penetrating out of nucleoid, will redefine later
# Below is the probablility of penetrationg for the higher species
# if everything stays where it is initially, this is zero
# if the higher penetrability is free to move into the new region, this is 1
pen_higher = float(sys.argv[3]) 


dt = .001 # seconds
t_tot = 40.0  # total time in s
D = 67411.95381724922 # nm^2/s
D_out = D
D_ratio = float(sys.argv[4])
D_in = D_ratio * D_out

if pen >= 1.0:
    pen_in = pen_higher  # probability of penetrating into the nucleoid
    pen_out = pen_in / (pen * D_ratio**1)  # probability of penetrating outside the nucleoid
    if pen>1e4:
        pen_out = 0 
else:
    pen_out = pen_higher  # probability of penetrating into the nucleoid
    pen_in = pen_out * (pen * D_ratio**1)  # probability of penetrating outside the nucleoid

print('pen in ', pen_in, 'pen out ', pen_out)
exp_name = sys.argv[5]

print(exp_name)

# Do the random walk:
print('Parameters for cell: ')
print(l_cell, r_cell, r_cap, l_nuc, r_nuc)
print('Parameters for random walk: ')
print(D_in, D_out, pen_in, pen_out, R_gem)
print('')
print('Performing random walk')

pos, t, tot_out = ecoli_walk(pos_init, dt/1, t_tot/1, D_in, D_out, pen_in, pen_out, l_cell, r_cell, r_cap, l_nuc, r_nuc, R_gem, gap_min, gap_max)

# Save everything
save_directory = sys.argv[6]
np.savetxt(save_directory+'t_' + exp_name + '.txt', t)
np.savetxt(save_directory+'tot_out_' + exp_name + '.txt', tot_out)

import os
try:
    os.makedirs(save_directory+exp_name)
except FileExistsError:
    # directory already exists
    pass

for i in range(N):
    pos_n = pos[:, i, :]
    np.savetxt(save_directory+exp_name + '/' + 'pos_' + str(i) + '_' + exp_name + '.txt', pos_n)

print('')
# Calculate MSD:
print('Calculating MSD')
t_lag, MSD = calc_msd(pos, dt, lag_step=5, t_end_frac=0.30)

np.savetxt(save_directory+'t_lag_' + exp_name + '.txt', t_lag)
np.savetxt(save_directory+'MSD_' + exp_name + '.txt', MSD)
