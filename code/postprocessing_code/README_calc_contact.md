## Contact Distribution Calculation Around bGEMS

### Description:
This script calculates the contact distribution of particles around bGEMs (large GEMs) from colloidal whole-cell simulation data. It analyzes particle positions over multiple simulation snapshots, computes distances between bGEMs and other particle types (DNA, ribosomes, and positive and negative crowder molecules), and bins the contact distributions into histograms. Output files store the results for further analysis.

### Inputs:
1. **sim_folder**: Path to the folder containing simulation data, specifically LAMMPS trajectory files in the `dump` subfolder.
2. **write_name**: Base name for the output files where contact distributions will be written.

### Usage:
Run the script using Python with two command-line arguments:

```
python calc_contact.py <simulation_folder> <output_file_base>
```

- `simulation_folder`: Directory containing LAMMPS trajectory files (`dump/` subfolder).
- `output_file_base`: Prefix for output files.

Example:
```
python calc_contact.py ./simulation_data/ results/contact_output
```

### Outputs:
Four text files are generated based on the particle type:
1. `<output_file_base>_DNA_all.txt`: Contact distribution for DNA particles.
2. `<output_file_base>_ribosome.txt`: Contact distribution for ribosome-sized particles.
3. `<output_file_base>_negcrowd.txt`: Contact distribution for negatively crowding particles.
4. `<output_file_base>_poscrowd.txt`: Contact distribution for positively crowding particles.

Each file contains binned contact distribution histograms for each snapshot.

### Key Parameters (change these in script if necessary):
- `R_nuc = 20.0`: Nucleus radius (simulation units).
- `R_cell = 27.4`: Cell radius (simulation units).
- `a_part = 1.42857`: Particle radius in simulation units.
- `diam_gem = 40`: Diameter of the bGEM in nm.
- `cutoff = 1/14`: Contact distance threshold in simulation units.
- `snapshot_start` and `snapshot_end`: Range of snapshots to analyze.

### Notes:
- Input trajectory files should be in LAMMPS dump format and named as `cell-pd.1.X.lammpstrj`.
- Ensure all necessary numpy files (`particles_X.npy`) are generated before running the contact analysis.
- The script assumes a consistent particle ID scheme across snapshots.

### Dependencies:
- Python 3.x
- NumPy

### Author:
The code was written by Alp M. Sunol.

For additional assistance, contact:  
**Alp M. Sunol**  
Email: asunol@seas.harvard.edu  
or\
**Roseanna N. Zia**\
Email: rzia@missouri.edu
