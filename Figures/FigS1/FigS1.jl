#=
    Figure S1
=#

using Model

using CSV
using DataFrames
using DataFramesMeta
using Plots
using StatsPlots

include("../figure_defaults.jl")

##############################################################
## Run/load analysis
##############################################################

include("../../Results/EstimateBinding.jl")

##############################################################
## Fig S1
##############################################################

figS1 = plot(log_γ, LLopt .- maximum(LLopt),
    ylim=(-3.0,0.5),
    frange = -4.0,
    fα = 0.2,
    lw = 2.0,
    c  = :black,
    box = :on,
    widen=false,
    legend=:none,
    xlabel="log(γ)",
    ylabel="Profile Loglikelihood"
)
hline!(figS1,[-1.92],ls=:dash,lw=2.0,c=:red)

##############################################################
## Save...
##############################################################

savefig("$(@__DIR__)/FigS1.svg")