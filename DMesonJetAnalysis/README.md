# How to run analysis on fast simulations

The command to execute the analysis is:

 `./runDataAnalysis.py FastSimAnalysis.yaml --gen powheg --proc PROC --ts TS --stage X`

where `PROC` is either `charm` or `beauty`, `TS` is the timestamp of the simulation (assigned automatically when the jobs are submitted to the grid, it is the number at the end of the directory where you find the output files on the grid), and `X` is the merging stage.
In order to make it work you'll have to setup `FastSimAnalysis.yaml` so that it points to the correct location where you downloaded the results of the simulation. This is set in the field `input_path`. The script then expect the following directory structure:

`{input_path}/FastSim_powheg_{PROC}_{TS}/stage_{X}/output`

Inside output the script expects to find many directories 001, 002,.. and inside each of them a file named `AnalysisResults_FastSim_powheg_{PROC}_{TS}.root`

The script will produce automatically a lot of pdf files and a root file will be saved in `{input_path}/FastSim_powheg_{PROC}_{TS}/stage_{X}/output/FastSimAnalysis_powheg_{PROC}_{TS}.root`
