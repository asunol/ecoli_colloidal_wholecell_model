These files are for the event-based long-time simulations in a confined spherocylinder. 

The MSD files are the mean squared displacements, t is the time of the simulation, t_lag is the time lags at which MSD was calculated, and tot_out is the number of bGEMs that remain outside the nucleoid.  Each simulation generated 500 trajectories.

The ending 'free' means the probability of going in or out is 1 (whichever is higher), and the other is 1/penatrability.  For example, the 20nm have a penatrability of 2.5, meaning they are 2.5 times more likely to be inside (normalized by volume).  Therefore, the probability that it goes into the nucleoid is 1.0, and the probability that a move goes out going accepted is 1.0/2.5.

The ending 'confined' means the particles are confined to inside or outside the nucleoid depending on their initial position.  So a particle inside remains inside and outside remains outside.  Their initial position matches the time_in/time_total in the experiments.  For example, roughly 60% of 20nm are inside the nucleoid.

The ending 'in' means all particles are inside the nucleoid region.

All simulations except "in" obey the time_in/time_total results from experiments.

The following are the diffusion coefficients determeined from these random walks:
D_50nm = 39654.090480734834 #nm^2/s
D_15q = 12256.718875863495 #nm^2/s
D_18q = 25927.67454509585 #nm^2/s
D_40nm = 67411.95381724921 #nm^2/s
D_20nm = 79075.6056507322 #nm^2/s

The following are the nucleoid size parameters that were used and gave correct plateuas:
l_cell = 700.0  # nm, cylindrical length of cell not including caps (divided by 2)
r_cell = 400.0  # nm, cylindrical radius of cell
r_cap = 400.0 # nm, radius pole caps

l_nuc = 650  # nm, cylindrical length of nucleoid not including caps (divided by 2)
r_nuc = 300.0  # nm, cyclindrical radius of nucleoid

We also ran additional simulations to represent exponentially growing cells.  These also had an additional gap size, r_gap, that allowed us to represent a bilobe nucleoid structure.  These were only run with 50 nm bGEMs to compare to experimental data, and bGEMs were completely excluded for the nucleoid for both runs.  The parameters are:
D_50nm_exp = 
l_cell = 1000.0  # nm, cylindrical length of cell not including caps (divided by 2)
r_cell = 400.0  # nm, cylindrical radius of cell
r_cap = 400.0 # nm, radius pole caps

l_nuc = 950  # nm, cylindrical length of nucleoid not including caps (divided by 2)
r_nuc = 300.0 or 390.0  # nm, cyclindrical radius of nucleoid, 390.0 blocks bGEMs from migrating from the interlobe to the cell caps and vice-versa.
r_gap = 400.0 # nm
D_50nm_exp = 51855.349090191696 #nm^2/s

All trajectories and MSD calculations were generated in the file ecoli_long_time_diffusion.py.  We ran simulations scaled at the 40 nm -7q time and rescaled the time by the appropriate diffusion coefficient to convert MSD to the appropriate value.  For more detail please see python file and documentation as well as the plotting notebook.