#=

    Model.jl

    A Julia module containing code to solve and simulate the model.

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=#
module Model

    using ColorSchemes
    using Distributions
    using Interpolations
    using KernelDensity
    using LinearAlgebra
    using MCMCChains
    using Plots
    using Random
    using SpecialFunctions
    using Statistics
    using StatsBase
    using StatsFuns
    using StatsPlots
       
    export solve_ode_model, solve_ode_model_all_vars, solve_alt_ode_model
    export GammaAlt, LogNormalAlt, GaussianCopula, TCopula, Copula, MvDependent, quantile_interpolated, vec_to_cor, marginalize
    export plot_fit_univariate, plot_fit_bivariate, plot_distribution, plot_distributions, plot_biv_distribution, plot_biv_distributions
    export simulate_model, simulate_model_noiseless, parameters_to_antibody
    export col_Q̄,col_Ū,col_Q,col_U,col_Q̄Ū,col_QU
    export pdf_ci

    col_Q̄ = "#F59640"  # Quenchable, quenched sample
    col_Ū = "#1BA7AD"  # Unquenchable, quenched sample
    col_Q = "#F54D3A"  # Quenchable, unquenched sample
    col_U = "#2339AC"  # Unquenchable, unquenched sample
    col_Q̄Ū = "#c56f88"
    col_QU = "#a8216b"

    include("deterministic.jl")
    include("distributions.jl")
    include("plotting.jl")
    include("statistical.jl")
    
end