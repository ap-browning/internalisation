#=
    Figure 5
=#

using Inference
using Model

using ColorSchemes
using JLD2
using Latexify
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
## Table S1
##############################################################

tabS1 = DataFrame(summarystats(C))[:,[1,2,3,6,7]]
tabS1.bestfit = θ

# Reorder
tabS1 = tabS1[!,[1,6,2,3,4,5]]


# Rounding
for i in 1:ncol(tabS1)
    if typeof(tabS1[1,i]) <: Number
        tabS1[:,i] = round.(tabS1[:,i],sigdigits=4)
    end
end

res = latexify(tabS1,env=:tabular)
write("$(@__DIR__)/tabS1.tex",res)
