### Instructions for reproducing the numerical results in "Noisy subspace clustering via matching pursuits", 2016, by Michael Tschannen and Helmut Bölcskei

Required software: matlab, make, pdflatex

Required data and scripts:
- SSC implementation: Download "SSC_ADMM_v1.1.zip" from http://www.ccs.neu.edu/home/eelhami/codes/SSC_ADMM_v1.1.zip and move the extracted folder "SSC_ADMM_v1.1" to "ssc_mps"
- NSN implementation: Download "greedysc_codes.zip" from http://dhpark22.github.io/files/greedysc_codes.zip and move the extracted folder "greedysc_codes" to "ssc_mps"


#### Tables 1 and 2:

Run the script "faceclustering.m" in Matlab.

#### Figures 2--4:

Cd to the folder "sscompmp_supp" and type
- make fig2 for Figure 2
- make fig3 for Figure 3
- make fig4 for Figure 4
- make all for Figures 2--4

To generate the figures without simulation, move the files from "simulated_data" to the parent folder and generate the figures following the instructions above.

