# Boltzmann Solver

This repository contains semi-classical methods for studying coupled spin and charge transport through exact solutions of the Boltzmann equation. These solutions are achieved via the exact diagonalization of the relaxation matrix in the collision integral. The provided codes were used to study the conversion between charge currents and spin signals and their lifetimes in tellurium crystals, as detailed in the paper
titled "Efficient spin accumulation carried by slow relaxons in chiral tellurium"
([arXiv link](https://arxiv.org/abs/2407.01187)). 

## Main Script

The main script exactly solves the Boltzmann transport equation for applied electric fields in chiral tellurium crystals (see "Exact Boltzmann transport approach" section in the paper for details). These solutions allow us to calculate the time evolution of the electron distribution function and its relaxon spectral decomposition, thus determining current-induced spin and orbital momentum accumulation and their lifetime. To run the code in MATLAB, use the following command:

```
find_EF_dependence
```

## Plotting and Analyzing Observables

The output data is saved in a "run_name.mat" file, which includes observables such as:
- Charge currents
- Current-induced spin and orbital angular momentum density
- Charge-to-spin conversion efficiency
- Spin relaxation time
- Total induced deviation in the electron distribution function
- Current-induced magnetization

These observables are visualized in "Observables_currents.png", "Observables_SL.png", and "Observables.png". 

The relaxon spectral decomposition of the induced shift in the electron distribution function as a function of the chemical potential is plotted in "Gamma_vs_EF.png". Additionally, the relaxon spectral decomposition and the fits of times dependencies (e.g., of the induced spin density at each intermediate value of the chemical potential) are stored in the "data_temp_observables" folder.

Please note that calculations in the paper were performed using the Habroc
cluster at the University of Groningen with a denser grid in momentum space, as described in the Methods section of the paper. The results
presented in the paper can be reproduced using the data and plotting procedures in the "Data_and_Plotting_procedures" folder. For instance, run "plot_FIG_2_3_and4_for_paper2.m" to plot observables.

## Additional Script for a single run using an effective kp-model

The main script works with ab initio wavefunctions obtained using the PAOFLOW package. However, the codes also support working with electron states from toy tight-binding models and low-energy kp models. To choose the model, modify the "model_ind" descriptor according to the description in the "Boltzmann_Solver.m" script. For example, you can reproduce the results for the kp model by using "inputfile_Boltzmann_kp.sh" as the input file in the "Single_run_for_kp_model" folder. This script needs to be run on a cluster due to the higher memory requirements for the exact diagonalization of the relaxation matrix.

## Contact

If you have any questions or require further clarification regarding
the procedure, feel free to reach out to the corresponding authors via
email. An extended Python version will be available soon.
