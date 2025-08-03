# spectrum analysis
This folder contains the details of spectrum analysis. 

# Requirement
MATLAB 2024b and above

# Usage
Once you clone the repo, you can try to re-generate the results shown on the manuscript. 

To generate the detailed results of the spectrum analysis: please first go to spectrum analysis folder, then use [`spectrum_analysis.m`](spectrum_analysis.m), the structure of the generated data will be shown as : 

```text
perturbation_eqpt_results/
  └── run_path (e.g., spectrum_3_leading_vehicles)
        ├── run_info.mat                # Summary for all 6 profiles (big tables)
        ├── run_folder (e.g., profile1)
        │     ├── profile1_data_summary.mat
        │     ├── 0001.mat              # nn = 1, run details for first beta combination
        │     ├── 0002.mat              # nn = 2, run details for second beta combination
        │     ├── 0003.mat              # nn = 3, run details for third beta combination
        │     └── ...                   # ...
        │     └── nn.mat                # nn = n^3, last beta combination
        └── run_folder (e.g., profile6)
              ├── profile6_data_summary.mat
              ├── 0001.mat              # nn = 1, run details for first beta combination
              └── ...                   # ...
              └── nn.mat                # nn = n^3, last beta combination
```
