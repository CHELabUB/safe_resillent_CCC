# Safe and Efficient, Data-driven Connected Cruise Control
This repo host the parameters and the source code for paper that is accepted to [MECC 2025](https://mecc2025.a2c2.org/)

[[arXiv]](https://arxiv.org/abs/2507.22227).

## Citation

```
@InProceedings{Xiao2025MECCSafeEfficientCCC,
      title={Safe and Efficient Data-driven Connected Cruise Control}, 
      author={Haosong Xiao and Chaozhe R. He},
      year={2025},
      eprint={2507.22227},
      archivePrefix={arXiv},
      primaryClass={eess.SY},
      url={https://arxiv.org/abs/2507.22227}, 
}
```
## Requirement
MATLAB 2024b and above

## Usage
Once you clone the repo, you can try to re-generate the results shown on the paper. 
You are also welcome to use our code to test on the data you collected on your own. 
The original data used in the paper may be available upon reasonable request.

To regenerate the paper result, please first go to the result folder. 

Second, to generate the bar plot, please use [`bar_plot.m`](/result/bar_plot.m); 

Third, to generate the speed profiles, please use [`plotter.m`](/result/plotter.m)

To generate the detailed results from optimization:

For energy simulation: please first go to energy_simulation folder, then use [`brute_force_search.m`](/energy_simulation/brute_force_search.m); 

For spectrum simulation: please first go to spectrum_analysis folder, then use [`spectrum_analysis.m`](/spectrum_analysis/spectrum_analysis.m); 

## Contact
[Chaozhe He](https://www.buffalo.edu/~chaozheh/) \
Assistant Professor \
Department of Mechanical and Aerospace Engineering \
University at Buffalo, The State University of New York
