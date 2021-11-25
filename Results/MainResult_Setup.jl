#=

    MainResult_Setup.jl

    Setup data, model, etc to produce the main result.

=#


using Inference
using Model

using CSV, DataFrames, DataFramesMeta
using Distributions
using StatsBase


##############################################################
## Data
##############################################################

# Load data
data = CSV.read("Data/CSV/DualMarker.csv",DataFrame)

# Autofluorescence distribution
E = [@subset(data,:Time .== 0.0, :Quenched .== false, :Temp .== 37.0).Signal_Cy5,
     @subset(data,:Time .== 0.0, :Quenched .== false, :Temp .== 37.0).Signal_BDP]

# Observation times
T = sort(unique(data.Time))[2:end]

# Format data to match `simulate_model` output
X = [[[@subset(data, :Time .== t, :Quenched .== qu, :Temp .== 37.0)[:,ch] 
        for ch in [:Signal_Cy5,:Signal_BDP]] 
        for qu in [true,false]] 
        for t  in T]

# Obtain quenching efficiency
η = 1 - (@df @subset(data,:Time .== 120.0,:Temp .== 4.0,:Quenched) mean(:Signal_Cy5)) / 
        (@df @subset(data,:Time .== 120.0,:Temp .== 4.0,(!).(:Quenched)) mean(:Signal_Cy5))

        
##############################################################
## Setup model
##############################################################

# Parameter distributions
function param_dist(ξ)
    μR,σR,μλ,σλ,ωλ,μβ,σβ,ωβ,ρRλ,ρRβ,ρ̄λβ = ξ

    # Marginal distributions
    R = LogNormalAlt(1.0,μR,σR)
    λ = GammaAlt(μλ,σλ,ωλ)
    β = GammaAlt(μβ,σβ,ωβ)

    # Construct copula from correlation matrix
    ρλβ = ρRλ * ρRβ + ρ̄λβ  * sqrt((1 - ρRλ^2) * (1 - ρRβ^2))
    C = GaussianCopula([ρRλ,ρRβ,ρλβ])

    # Return multivariate distribution (use interpolated quantile for λ and β)
    MvDependent(C,[R,λ,β],q=[false,true,true])
end

# Model
model = θ -> simulate_model(θ,T,param_dist,E;N=1000,η)


##############################################################
## Setup inference
##############################################################

# Pre-compute discrepancy metric
disc = makediscrepancy(X,w=[1.0,80.0])

# ABC distance metric
dist = θ -> (disc ∘ model)(θ)

# Prior
prior = Product([
    Uniform(5000,12000),# α₁
    Uniform(20.0,80.0), # α₂
    Uniform(0.0,10.0),   # σ₁
    Uniform(0.0,3.0),   # σ₂
    Uniform(0.0,1.0),   # μR
    Uniform(0.0,1.0),   # σR
    Uniform(0.0,0.5),   # μλ
    Uniform(0.0,0.3),   # σλ
    Uniform(-2.0,1.0),  # ωλ
    Uniform(0.0,0.20),  # μβ
    Uniform(0.0,0.10),  # σβ
    Uniform(-2.0,1.0),  # ωβ
    Uniform(-1.0,1.0),  # ρRλ
    Uniform(-1.0,1.0),  # ρRβ
    Uniform(-1.0,1.0),  # ρ̄λβ
    Uniform(0.0,0.15),  # p,
])

# Parameter names
param_names = ["α₁","α₂","σ₁","σ₂","μR","σR","μλ","σλ","ωλ","μβ","σβ","ωβ","ρRλ","ρRβ","ρ̄λβ","p"]