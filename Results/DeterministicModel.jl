#=

    DeterministicModel.jl

    Perform inference using the deterministic model. I.e., take the 'typical' approach
    and neglect parameter heterogeneity.

=#

using Inference
using Model

using CSV, DataFrames, DataFramesMeta
using Distributions
using ForwardDiff
using JLD2
using LinearAlgebra
using Optim
using StatsBase
using StatsPlots

##############################################################
## Data
##############################################################

# Load all data
data = CSV.read("Data/CSV/DualMarker.csv",DataFrame)

# Load summarised data
data_summary = @subset(CSV.read("Data/CSV/DualMarker_Summary.csv",DataFrame), :Temp .== 37.0, :Time .> 0.0)

# Observation times
T = sort(unique(data_summary.Time))

# Get pooled standard deviation
σ = pooled_std(data_summary,[:Quenched,:Time],[:Signal_Cy5])[1]

# Get quenching efficiency
η = 1 - (@df @subset(data,:Time .== 120.0,:Temp .== 4.0,:Quenched) geomean(:Signal_Cy5)) / 
        (@df @subset(data,:Time .== 120.0,:Temp .== 4.0,(!).(:Quenched)) geomean(:Signal_Cy5))

# Get autofluorescence
e = @df @subset(data,:Time .== 0.0, :Quenched .== false, :Temp .== 37.0) geomean(:Signal_Cy5)

##############################################################
## Fit ODE model to data
##############################################################

# Likelihood
function loglike(θ)
    loglike = 0.0
    for t in T
        # Solve ODE model
        A,I = solve_ode_model(t,θ[1:3])
        # Predicted unquenched signal
        m_unquenched = A + e / θ[4]
        # Predicted quenched signal
        m_quenched = I + (1 - η) * (A - I) + e / θ[4]
        # Data ~ N(model,σ²)
        loglike += loglikelihood(Normal(θ[4]*m_unquenched,σ),@subset(data_summary,:Time .== t,(!).(:Quenched)).Signal_Cy5)
        loglike += loglikelihood(Normal(θ[4]*m_quenched,σ),@subset(data_summary,:Time .== t,:Quenched).Signal_Cy5)
    end
    return loglike
end

# Obtain maximum likelihood estimate
θdet = optimize(p -> -loglike(p),ones(4)).minimizer

# Convert to named tuple for easy manipulation
θdet = NamedTuple{(:λ,:β,:p,:α)}(θdet)

##############################################################
## Confidence intervals
##############################################################

# Use observed FIM to approximate confidence intervals
H = ForwardDiff.hessian(loglike,collect(θdet))
se = sqrt.(diag(inv(-H)))
ci = [θdet[i] .+ [-1.0,1.0] * se[i] * quantile(Normal(),0.975) for i = 1:4]

# Tabulate estimates
tab1 = DataFrame(
    Parameter = [:λ,:β,:p,:α],
    Estimate  = round.(collect(θdet),sigdigits=3),
    CI95      = [round.(ci[i],sigdigits=3) for i = 1:length(ci)]
)

