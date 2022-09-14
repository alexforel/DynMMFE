# Dynamic MMFE
This code produces the results with synthetic data of the paper "Dynamic stochastic lot sizing with forecast evolution in
rolling-horizon planning".

This code has been developed on `julia 1.6.0`.

## How to reproduce the results:
The scripts that generate data and results are included in the root folder. The tables in the paper are obtained directly by running the scripts, which produces `.tex` files. The figures have been generated in tikz. The scripts provide their inputs.

The scripts are organized as follows:
* Figures 2, 3 and 7: the data is generated by running `illustrate_mmfe_and_demand_results.jl`,
* Figures EC.3.1, EC.3.2, and EC.3.3: the data is generated by running `run_repeated_simulations_synthetic_data.jl`,
* Tables 2, 3, EC.4.1, EC.4.2: the tables are generated by running `run_repeated_simulations_synthetic_data.jl`,
* Tables 4 is generated by running `run_sensitivity_analysis_correlation.jl`,
* Figures EC.5.1 and EC.5.2: the data is generated by running `run_sensitivity_analysis_capacity.jl`,
* Figures EC.5.3 and EC.5.4: the data is generated by running `run_sensitivity_analysis_scenario_structure.jl`,
* Tables EC.6.1 is generated by running `run_free_return_experiment.jl`.
