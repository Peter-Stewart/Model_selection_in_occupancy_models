This repository contains code to fully reproduce the simulations and analyses in our paper ___Model selection in occupancy models: inference versus prediction___. 

**File List:**

* **Fig2_simulation.R** - code to reproduce the simulations, analyses and graphs in Figure 2 of our manuscript
* **Occ_Mbias.R** - code to run the occupancy simulations with M-bias present in the occupancy process only
* **Det_Mbias.R** - code to run the occupancy simulations with M-bias present in the detection process only
* **Occ_Det_Mbias.R** - code to run the occupancy simulations with M-bias present in both the occupancy and detection processes
* **Figures.R** - code to analyse the results of the three occupancy simulations and produce figures 4-6 in the main text and S1-S6 in the supplementary material

The repository also contains a **Data** folder containing the results from the three occupancy simulations. These data are included only for convenience, as the code above is sufficient to replicate them in full. Metadata for these results files can be found in the **metadata.md** file contained in the **Data** folder.
