# E. coli Colloidal Whole-Cell Model

This repository contains the code and postprocessed data associated with the computational models in the publication:  
**Macromolecular interactions and geometrical confinement determine the 3D diffusion of ribosome-sized particles in live *Escherichia coli* cells**  
Authors: Valverde-Mendez*, Sunol*, et. al   
Accepted in *PNAS (2024)*  


## Key Components

1. **Code**
   - `event_based_longtime_simulation_code`: Simulates long-time molecular trajectories in confined environments.
   - `postprocessing_code`: Includes scripts to compute contact distributions and MSD.
   - `packmol_code`: Packmol files to create initial configurations for whole-cell simulations.


2. **data**
   - `MSD_data`: Mean Squared Displacement (MSD) data for colloidal whole-cell model trajectories, long-time event-based simulation trajectories, and experimental trajectories (experimental data was obtained by Diana Valverde-Mendez and is included here only for plotting the comparison to simulations)
   - `coordination_number`: Contact distribution data from colloidal whole-cell simulations.
   - `nucleoid_voronoi_distribution`: Void width distributions from Voronoi analysis.
   - `bGEM_in_v_out`: localization data on whether a bGEM is inside or outside the nucleoid

3. **JupyterNotebooks**
   - Jupyter notebooks for reproducing plots and visualizations in the manuscript and supplement.

## Repository Structure
```plaintext
.
├── code/                           # Source code for simulations and analysis
│   ├── event_based_longtime_simulation_code/   # Event-based simulation code
│   ├── packmol_code/               # Packmol setup scripts
│   └── postprocessing_code/        # Scripts for contact distribution analysis, localization, and MSD
│
├── data/                           # Input data and results
│   ├── bGEM_in_v_out/              # bGEM localization data
│   ├── coordination_number/        # Particle coordination number analysis
│   ├── MSD_data/                   # MSD analysis results
│   │   ├── event_based_long_time_sim_data/     # Long-time simulation MSD results
│   │   ├── Experimental_Data/      # Experimental MSD comparison
│   │   └── whole_cell_colloidal_sim_data/      # Colloidal whole-cell MSDs
│   └── nucleoid_voronoi_distribution/          # Void size data of nucleoid
│
└── JupyterNotebooks/           # Jupyter notebooks for visualization and analysis

```


## Dependencies

To run the scripts, the following dependencies are required:
- Python (>= 3.8)
- NumPy
- Pandas
- Matplotlib
- Jupyter Notebook

## Usage

All associated scripts and data have README.md or README.txt files with detailed instructions and descriptions.

## Citation
If you use any of the data or code, please cite both the associated manuscript and this repository:  

**Alp M. Sunol, Jennifer L. Hofmann, Roseanna N. Zia**  
**E. coli Computational Whole-Cell Model**  
Database: GitHub  
Deposited on: December 13, 2024.  


---

*For questions or issues, please contact [Roseanna N. Zia](mailto:rzia@missouri.edu).*