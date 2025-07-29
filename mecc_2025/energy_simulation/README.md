# energy simluation
This folder contians the simulator of efficient CCC. 

# Requirement
MATLAB 2024b and above

# Usage
Once you clone the repo, you can try to re-generate the results shown on the manuscript. 

To generate the detailed results for the energy optimization: please first go to energy_simulation folder, then use [`brute_force_search.m`](brute_force_search.m), the structure of the generated data will be shown as : 

```text
result/
  └── run_path (e.g., CCC_car_filter_1)
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


