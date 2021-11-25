# Internalisation

Repository for the preprint "Identifying cell-to-cell variability in internalisation using flow cytometry" available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.11.24.469957v1).

## Data

All data used in the analysis is available in `Data`. Due to GitHub size limitations, the script `ProcessData.jl` must be run to process individual gated `.fcs` files into a `CSV` file used for analysis.


## Running the code
 
### Installation

Code to produce all results is contained within this repository, including the Julia modules `Inference` and `Model`. To download, first clone this Github repo. Next, add the two module folders to your `LOAD_PATH`:
```
push!(LOAD_PATH,"/path/to/repo/Module/Inference")
push!(LOAD_PATH,"/path/to/repo/Module/Model")
```
Next, run the following code (press `]` to enter the `Pkg` REPL) to install all required packages and activate the project
```
(v1.6) pkg> activate .
(Spheroids) pkg> instantiate
```

### Results

Scripts to produce the results are available in the `Results` folder. For example, to produce the main results, run `MainResult_Inference.jl` (takes approximately 24 h to perform ABC on a standard, quad-core machine). For convenience, this script saves its output to `MainResult.jld2`.

Code used to produce each figure in the main document and supplementary material is available in the `Figures` folder. For example, to reproduce Figure 4, run `Figures/Fig4/Fig4.jl`.
