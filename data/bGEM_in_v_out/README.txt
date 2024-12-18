These files contain data for the volume fraction inside, outside, and total for the bGEMs throughout the simulation.

They were generated using the postprocessing script phi_v_t.py.

The data format is:
snapshot phi_in phi_out phi_total

We have 200 snapshots for 20nm and negatively charged 40 nm bGEM variants.

For 50 nm and positevely charged bGEMs, we ran the simulation and analysis for longer (500 snapshots) to collect more appropriate statistics, as their motion was slower so they naturally need to diffuse for more time to sample enough of the cell volume.

Each snapshot is taken every 100000 simulation timesteps such that the bGEMs diffuse a significant amount within this time window.