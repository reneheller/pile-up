# pile-up
A gnuplot script for Monte Carlo simulations of disk and stellar tidal torques acting on hot Jupiters

To run this script through a terminal, first open gnuplot (e.g. bytyping "gnuplot" into your terminal and then pressing the RETURN key). Within the gnuplot environment, run the script via

gnuplot>i=1; call "pile-up.gp" 1000 pdf

This will generate 1000 curves of the disk and stellar tidal torques acting on a Jupiter-mass planet as a function of distance to its young, solar-type host star.

Details are described in the research paper by René Heller (2018), "Formation of hot Jupiters through Disk Migration and Evolving Stellar Tides", Astronomy & Astrophysics. This gnuplot script was used to generate Figure 3b in this paper.
