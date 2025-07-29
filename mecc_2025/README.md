# Safe and Efficient, Data-driven Connected Cruise Control
This repo host the parameters and the source code for manuscript that is accepted to MECC 2025

# Requirement
MATLAB 2024b and above

# Usage
Once you clone the repo, you can try to re-generate the results shown on the manuscript. You are also welcome to use our code to test on the data you collected on your own. (The original data used on the manuscript will be available upon request)

To regenerate the paper result, please first go to the result folder. 

Second, to generate the bar plot, please use [`bar_plot.m`](/result/bar_plot.m); 

Third, to generate the speed profiles, please use [`plotter.m`](/result/plotter.m)

To generate the detailed results from optimization:

For energy simulation: please first go to energy_simulation folder, then use [`brute_force_search.m`](/energy_simulation/brute_force_search.m); 

For spectrum simulation: please first go to spectrum_analysis folder, then use [`spectrum_analysis.m`](/spectrum_analysis/spectrum_analysis.m); 

# Contact:
[Chaozhe He](https://www.buffalo.edu/~chaozheh/) \
Assistant Professor \
Department of Mechanical and Aerospace Engineering \
University at Buffalo, The State University of New York
