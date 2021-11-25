#=
    MainResult.jl

    Perform ABC-SMC, then ABC-MCMC, to produce data for the main results.

    WARNING: This script takes ~1 day to run on a 4-core machine. For convenience, the 
    results are saved in the `JLD2` file `MainResult.jld2`.
    
=#

using Inference
using Model

using CSV, DataFrames, DataFramesMeta
using Distributions
using JLD2
using Plots
using StatsBase
using StatsPlots


##############################################################
## Setup problem
##############################################################

include("SupportingResult_OnlyCy5_Setup.jl")

##############################################################
## ABC SMC
##############################################################

P = abc_smc(dist,prior,maxsteps=30;n=1000,ptarget=0.005,α=0.75)

##############################################################
## Determine discrepancy function
##############################################################

# The threshold is chosen so the "best fit particle" will be accepted ~50% of the time
D = [dist(minimum(P).θ) for i = 1:100]
ε = quantile(D,0.5)
K = d -> logistic(d,ε,20.0)

##############################################################
## ABC MCMC
##############################################################

# Pilot chains (1,000,000 iterations, skip=100)
@time C_pilot = abc_mcmc(dist,prior,
    minimum(P).θ,
    2.38^2 * cov(P) / ndims(P),
    K,
    n=10000,param_names=param_names,nchains=4,thinning=100
)

# Tuned chains (10,000,000 iterations, skip=100)
@time C = abc_mcmc(dist,prior,
    minimum(P).θ,
    2.38^2 * cov(Matrix(C_pilot)') / ndims(P),
    K,
    n=100000,param_names=param_names,nchains=4,thinning=100
)


##############################################################
## Find best fit
##############################################################

Θ = Matrix(C[1:1000:end])
D = [mean(dist(θ[:]) for i = 1:100) for θ in eachcol(Θ)]
θ = Θ[:,findmin(D)[2]]

##############################################################
## Save results
##############################################################

@save "Results/SupportingResult_OnlyCy5.jld2" P ε C_pilot C θ