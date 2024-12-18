# README: Packmol Files for Initial Configuration Generation

This folder contains Packmol scripts and associated PDB files used to generate initial configurations for our whole-cell colloidal simulations.

## System Description  
- **Nucleoid**: Spherical, self-assembled at a volume fraction of 0.10 with 10kT interactions.  
- **bGEMs**: 40 nm in size in this example.  
- **Polysomes**: Consist of 6 ribosomes in a linear shape. Although initially linear, the polysomes are allowed to relax into their equilbrium confirmation during the actual simulation before recording any final data.  

## **Files**  
The following files are included in this folder:  

- **rib.pdb**: Ribosome structures.  
- **polysome.pdb**: Linear polysome structures consisting of 6 ribosomes.  
- **frozen_nucleoid_10kT_f10.pdb**: Example initial nucleoid configuration frozen at 10kT and 0.10 volume fraction.  
- **GEM.pdb**: File for bGEM structures.  
- **crowder_positive.pdb**: Positively charged crowder molecules.  
- **crowder_negative.pdb**: Negatively charged crowder molecules.  
- **cell_final_f10_10kT_40nm.inp**: Packmol input script for generating the initial configuration.  


## Workflow  
1. **Preliminary Simulations (Without bGEMs)**:  
   - Initial configurations were generated using Packmol.
   - These configurations were further equilibrated in LAMMPS to determine the equilibrium distribution of molecules inside and outside the nucleoid.  

2. **Adding bGEMs**:  
   - The system was re-initialized with bGEMs using the equilibrium distributions obtained from the preliminary simulations.  
   - The initial distribution of bGEMs was chosen to match experimental heatmaps and occupation times.  

3. **Dynamic Simulations**:  
   - Dynamic simulations were run in LAMMPS to ensure the system reached equilibrium before any measurements were recorded for final data analysis.  
   - **Note**: The relative abundances of molecules inside and outside the nucleoid in the Packmol scripts thus may not exactly match the values observed during final equilibrium simulations.

## Placement of Molecules  
- Without explicit placement, Packmol tends to position molecules outside the nucleoid.  
- To avoid this, specific amounts of each molecule were explicitly placed inside the nucleoid, while the rest were placed generally within the cell.  
- As an extra precaution, molecules not in the nucleoid could also have been explicitly placed outside the nucleoid; however, we tried this and found this step unnecessary as it did not affect the final distribution that packmol generates.  
