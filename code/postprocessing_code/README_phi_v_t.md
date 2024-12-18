## Volume Fraction Calculation for Each Species

### Description:
This script calculates the **volume fraction** of different particle species both inside and outside the nucleoid region from the colloidal whole-cell simulation trajectories. It processes particle positions from LAMMPS dump files and computes their volume fractions in the nucleoid, nucleoid-free cytoplasm, and the entire cell.

### Inputs:
1. **sim_folder**: Path to the folder containing simulation data (LAMMPS dump files in `dump/` subfolder).
2. **sim_name**: Base name for output files.
3. **n_types**: Number of particle types to track.
4. **gem_diameter**: Diameter of GEM particles (in nanometers).

### Usage:
Run the script using Python with four command-line arguments:

```
python phi_v_t.py <simulation_folder> <simulation_name> <n_types> <gem_diameter>
```

- `simulation_folder`: Directory containing LAMMPS trajectory files (`dump/` subfolder).
- `simulation_name`: Prefix for output files.
- `n_types`: Total number of particle species to analyze.
- `gem_diameter`: Diameter of GEMs in nanometers.

Example:
```
python phi_v_t.py ./simulation_data/ run1 4 40.0
```

### Outputs:
Four text files are generated, one for each particle type:
1. `phi.<sim_name>.GEM.txt`: Volume fractions for GEM particles.
2. `phi.<sim_name>.ribosome.txt`: Volume fractions for ribosome particles.
3. `phi.<sim_name>.neg_prot.txt`: Volume fractions for negatively crowding proteins.
4. `phi.<sim_name>.pos_prot.txt`: Volume fractions for positively crowding proteins.

Each file contains the following columns for each simulation step:
- Step number
- Volume fraction **inside** the nucleoid
- Volume fraction **outside** the nucleoid (cytoplasm)
- Total volume fraction in the cell

### Key Parameters (change these in script if necessary):
- `R_nuc = 20.0`: Nucleoid radius (simulation units, where 1 unit = 14 nm).
- `R_cell = 27.4`: Cell radius in simulation units.
- `V_nuc`, `V_cyt`, `V_cell`: Computed volumes of the nucleoid, cytoplasm, and cell.
- `gem_radius`: GEM radius derived from input diameter.
- `a_part`: Radii of tracked particle types [GEM, ribosome, negative crowders, positive crowders].

### Notes:
- Input trajectory files should be in LAMMPS dump format and named as `cell-pd.1.X.lammpstrj`.
- Particle types are expected to follow sequential IDs starting from 6 (bGEMs).  First 5 are DNA particles. Can adjust this part if necessary.
- The script assumes particles are tagged correctly and consistently across snapshots.

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
