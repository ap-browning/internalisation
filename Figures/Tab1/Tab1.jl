#=
    Table 2
=#

using Inference
using Model

using JLD2
using Plots
using StatsPlots
using Random
using Latexify

include("../figure_defaults.jl")

##############################################################
## Load model and results
##############################################################

include("../../Results/DeterministicModel.jl")

##############################################################
## Make table
##############################################################

res = latexify(tab1,env=:tabular)
write("$(@__DIR__)/tab1.tex",res)