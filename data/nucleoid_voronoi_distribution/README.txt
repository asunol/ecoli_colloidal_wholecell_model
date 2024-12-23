These files contain the voronoi void width distributions of the nucleoid model throughout the paper.

VoroHist.txt has the histogram for of the void sizes (column 1) for the heterogeneous (5 kT, column 2) and homogeneous (10kT, column 4) nucleoid structures (as described in the SI to the paper) at a volume fraction of 0.10. The other columns are a 6 kT (column 3) nucleoid and then the same set of interaction strneghts at volume fraction 0.20 (columns 5-7) and 0.30 (columns 8-10), which are not presented in the main manuscript.

For the nucleoid used in the whole-cell simulations (10 kT, phi_n = 0.10), we also calculated the void-width distribution at different times as described in the SI. For these, the void width of each individual void was stored, and the files are named edge_voids.tstep.txt. Please see the python notebook for the conversion of these files into a histogram.

All calculations were performed using methods described in Ryu and Zia, J. Rheol, 2022, which utilizes the library Voro++: https://math.lbl.gov/voro++/
