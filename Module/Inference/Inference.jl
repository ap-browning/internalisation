#=

    Inference.jl

    A Julia module code to perform ABC SMC with CDF/correlation matching.

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=#
module Inference

    using AbstractMCMC
    using AdvancedMH    
    using DataFrames
    using DataFramesMeta
    using Distributions
    using Interpolations
    using LinearAlgebra
    using MCMCChains
    using .Threads
    using StatsBase
    using Statistics
    using StatsPlots
    using Printf
    
    export abc_smc, distances, locations, abc_mcmc, acceptance_rate, ess, logistic, logit
    export makediscrepancy, makediscrepancy_univariate
    export plot_posteriors, plot_chains
    export pooled_data, pooled_std, pooled_cov, pooled_cor
    export Particle

    include("abc.jl")
    include("particles.jl")
    include("discrepancies.jl")
    include("plotting.jl")
    include("pool_data.jl")

end