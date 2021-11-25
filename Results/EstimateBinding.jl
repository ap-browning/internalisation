#=

    EstimateBinding.jl

    Perform inference using the deterministic model where the antibody binding
    rate, γ, is also estimated (this model is non-identifiable).

=#

using Model

using CSV, DataFrames, DataFramesMeta
using Distributions
using ForwardDiff
using JLD2
using LinearAlgebra
using Optim
using Roots
using StatsBase
using StatsPlots

##############################################################
## Run DeterministicModel.jl
##############################################################

include("DeterministicModel.jl")

##############################################################
## Fit ODE model to data
##############################################################

# Likelihood
function loglike_alt(θ)
    loglike = 0.0
    for t in T
        # Solve ODE model
        A,I = solve_alt_ode_model(t,θ[1:4])
        # Predicted unquenched signal
        m_unquenched = A + e / θ[5]
        # Predicted quenched signal
        m_quenched = I + (1 - η) * (A - I) + e / θ[5]
        # Data ~ N(model,σ²)
        loglike += loglikelihood(Normal(θ[5]*m_unquenched,σ),@subset(data_summary,:Time .== t,(!).(:Quenched)).Signal_Cy5)
        loglike += loglikelihood(Normal(θ[5]*m_quenched,σ),@subset(data_summary,:Time .== t,:Quenched).Signal_Cy5)
    end
    return loglike
end

# Obtain maximum likelihood estimate (just about meaningless, since model is non-identifiable)
θalt = optimize(p -> -loglike_alt(p),ones(5)).minimizer

##############################################################
## Profile (log) binding rate γ
##############################################################

log_γ = range(-1.0,5.0,length=50)
LLopt = [-optimize(p -> -loglike_alt([γ;p]),collect(θdet)).minimum for γ in exp.(log_γ)]

# Get 95% CI
pll = γ -> -optimize(p -> -loglike_alt([γ;p]),collect(θdet)).minimum - maximum(LLopt)
find_zero(γ -> pll(γ) + quantile(Chisq(1),0.95) / 2, [0.0,1.0])