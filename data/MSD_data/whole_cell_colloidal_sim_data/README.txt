MSD of bGEMs for colloidal whole-cell simulations.  The .in and .out files are from trajectories that remain completely inside and outside the nucleoid during the duration for at least the duration of the specified lag time.  The simulation timestep is output in the file that doesn't have '.in' or '.out', which is the MSD of the entire trajectory, whether the bGEM is inside or outside the nucleoid and whether it switched regions during the trajectory or not.  

Outputs are in LAMMPS lennard-jones units.  For unit conversions, please see plotting notebook.

MSD was calculated with the file MSD_gem.py