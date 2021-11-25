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

include("../../Results/MainResult_Setup.jl")
@load "Results/MainResult.jld2"

##############################################################
## Make table
##############################################################

sigdigits = 4
tab2 = DataFrame(
    "Parameter"     => param_names,
    "Lower"         => minimum.(prior.v),
    "Upper"         => maximum.(prior.v),
    "Estimate"      => round.(θ;sigdigits)
)
res = latexify(tab2,env=:tabular)
write("$(@__DIR__)/tab2.tex",res)