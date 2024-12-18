## MSD Calculation for bGEMS

### Description:
This script calculates the mean-squared displacement (MSD) of bGEMS from the colloidal whole-cell simulations. It analyzes their trajectories over time, computes their MSD, and separates the analysis into particles that are **inside** and **outside** a defined nucleoid region.

Additionally, it saves ribosome positions for further postprocessing if needed.

### Inputs:
1. **sim_folder**: Path to the folder containing simulation data, specifically LAMMPS trajectory files (`dump/` subfolder).
2. **sim_name**: Base name for output files.

### Usage:
Run the script using Python with two command-line arguments:

```
python MSD_gem.py <simulation_folder> <simulation_name>
```

- `simulation_folder`: Directory containing LAMMPS trajectory files (`dump/` subfolder).
- `simulation_name`: Name prefix for output files.

Example:
```
python MSD_gem.py ./simulation_data/ run1
```

### Outputs:
The script generates the following files:
1. `MSD.gems_only.<sim_name>.txt`: Overall MSD for all GEMS.
2. `MSD.gems_only.<sim_name>.in.txt`: MSD for GEMS **inside** the nucleoid.
3. `MSD.gems_only.<sim_name>.out.txt`: MSD for GEMS **outside** the nucleoid.
4. `gems.npy`: Numpy file storing GEM positions at all timesteps.
5. `ribs.npy`: Numpy file storing ribosome positions at all timesteps.

### Key Parameters (change these in script if necessary):
- `R_nuc = 20.0`: Nucleoid radius (simulation units).
- `step_size = 10000`: Timestep interval for analysis.
- `file_size = 100000`: Number of steps per LAMMPS dump file.
- `start_num` and `end_num`: Range of LAMMPS dump files to analyze.
- `a_part = 1.42857`: GEM particle radius in simulation units.

### Notes:
- Input trajectory files should be in LAMMPS dump format and named as `cell-gems.1.X.lammpstrj` for GEMs and `cell-ribs.1.X.lammpstrj` for ribosomes. Code can be adjusted if output is all particle types by filtering for the desired particles.
- Ensure consistent particle numbering across snapshots for accurate reordering.
- The script assumes particle types are identified in the trajectory files.

### Dependencies:
- Python 3.x
- NumPy


### Author:
The code was written by Alp M. Sunol.

For additional assistance, contact:  
**Alp M. Sunol**  
Email: asunol@seas.harvard.edu  
or  
**Roseanna N. Zia**  
Email: rzia@missouri.edu