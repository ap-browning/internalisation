#=
    Figure 4
=#

using Inference
using Model

using JLD2
using Plots
using StatsPlots
using Random

include("../figure_defaults.jl")

##############################################################
## Load model and results
##############################################################

include("../../Results/MainResult_Setup.jl")
@load "Results/MainResult.jld2"

##############################################################
## Compute results
##############################################################

# Simulate model
Y = simulate_model(θ,T,param_dist,E;N=100_000,η)

##############################################################
## Figure 4(a,b)
##############################################################

fig4ab = plot_fit_univariate(X,Y,T)

##############################################################
## Figure 4(c,d)
##############################################################

fig4cd = plot_fit_bivariate(X,Y,T)

##############################################################
## Save...
##############################################################

savefig(fig4ab,"$(@__DIR__)/fig4ab.svg")
savefig(fig4cd,"$(@__DIR__)/fig4cd.svg")