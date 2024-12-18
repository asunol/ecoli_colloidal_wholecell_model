# Event-based long-time diffusion simulation in *E. coli*

This code simulates the random walk of particles inside an *E. coli* cell, taking into account the cell's geometry and the ability to define separate diffusion coefficients and entry/exit probabilities for particles inside and outside the nucleoid region. It provides functionality for initializing particle distributions, evolving their trajectories, and calculating mean-squared displacements (MSD) over time.

---

## **How to Use**

### Required Input Parameters
The script is executed with the following command-line arguments:

```bash
python ecoli_random_walk.py <R_gem> <pen> <pen_higher> <D_ratio> <exp_name> <save_directory>
```

- `<R_gem>`: Radius of the particle (GEM) in nanometers.
- `<pen>`: Overall penetrability ratio of the nucleoid versus cytoplasm. 
  - If `pen = 0`, all particles are outside the nucleoid.
  - If `pen > 1e4`, all particles are inside the nucleoid.
- `<pen_higher>`: Penetrability probability for high penetrability regions.
- `<D_ratio>`: Diffusion coefficient ratio inside the nucleoid relative to outside.
- `<exp_name>`: A unique name for the experiment to label the output files.
- `<save_directory>`: The directory path where the simulation outputs will be saved.

### Code Features
- **Cell and Nucleoid Geometry**:
  - The cell is represented as a cylinder with hemispherical caps.
  - The nucleoid is modeled as a cylindrical region with optional bilobe gaps.
- **Random Walk Simulation**:
  - Separate diffusion coefficients for nucleoid and cytoplasm.
  - Probability-based entry/exit rules for particles transitioning between regions.
- **MSD Calculation**:
  - Computes the mean-squared displacement for a fraction of the total simulation time.

### Output
- Time series of particle positions saved in individual text files.
- Aggregate data for lag times (`t_lag`) and MSD values (`MSD`).
- A count of particles in each region over time (`tot_out`).
- All outputs are saved in `<save_directory>/<exp_name>`.

### Example Command
```bash
python ecoli_long_time_diffusion.py 20.0 1.0 0.5 0.9 experiment1 ./results/
```

---

## **Description of Functions**

1. **`cell_region`**:
   - Determines whether particles are inside the nucleoid, cytoplasm, or outside the cell based on coordinates.

2. **`init_distribution`**:
   - Initializes the positions of particles throughout the cell with user-defined penetration and geometry.

3. **`ecoli_walk`**:
   - Simulates the time evolution of particles using a random walk approach. Handles transitions between regions based on penetration probabilities.

4. **`calc_msd`**:
   - Computes the mean-squared displacement of particles over time for specified lag steps.

---

## **Outputs**

- **Position Data**: 
  - `pos_<i>_<exp_name>.txt`: Particle position data for particle `i` in the simulation.
- **Region Totals**: 
  - `tot_out_<exp_name>.txt`: Number of particles in each region (nucleoid, cytoplasm, or outside).
- **MSD Data**:
  - `t_lag_<exp_name>.txt`: Lag times.
  - `MSD_<exp_name>.txt`: Mean-squared displacement data.
- **Time Data**:
  - `t_<exp_name>.txt`: Time steps for the simulation.

---

## **Important Notes**

- This code was developed to simulate particle dynamics inside an *E. coli* cell and requires fine-tuning for specific experimental setups.
- Make sure to adjust the **cell and nucleoid geometry parameters** as well as the **diffusion coefficient** in the script to match your desired conditions.
- The code assumes an isotropic diffusion process and models boundary interactions using hard constraints on movement.

---
## Author:
The code was written by Alp M. Sunol.
For additional assistance, contact:  
**Alp M. Sunol**  
Email: asunol@seas.harvard.edu  
or\
**Roseanna N. Zia**\
Email: rzia@missouri.edu