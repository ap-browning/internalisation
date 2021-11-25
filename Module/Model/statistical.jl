#=

    statistical.jl

    Contains code to simulate the statistical model.

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=#
"""
    simulate_model(θ,t,d,E;N=1000)

Simulate the full statistical model for 2N cells, returning `[[Q̄,Ū],[Q,U]]` where
`[Q̄,Ū]` are paired measurements of the quenched sample, and `[Q,U]` of the unquenched
sample.

Parameter values are given as a vector `θ` where
    `α₁,α₂,σ₁,σ₂ = θ[1:4]`
are the observation process parameters,
    `d(θ[5:end-1])`
provides a multivariate distribution parameterised by hyperparameters `θ[5:end-1]`, and
    `p = θ[end]`
is the proportion recycled.
"""
function simulate_model(θ::Vector,t::Number,d::Function,E;N::Int=1000,η::Number=1.0)

    # Parameters
    α₁,α₂,σ₁,σ₂ = θ[1:4]        # Observation process parameters
    ξᵢ = rand(d(θ[5:end-1]),2N) # Sample from parameter distribution
    p  = θ[end]                 # Proportion recycled
    Rᵢ = @view ξᵢ[1,:]
    λᵢ = @view ξᵢ[2,:]
    βᵢ = @view ξᵢ[3,:]
        
    # Sample Aᵢ and Iᵢ
    Aᵢ = similar(λᵢ)
    Iᵢ = similar(λᵢ)
    for i = 1:2N
        Aᵢ[i],Iᵢ[i] = solve_ode_model(t,[λᵢ[i],βᵢ[i],p])
    end

    # Sample errors
    idx1 = sample(1:length(E[1]),N)
    idx2 = sample(1:length(E[2]),N)

    V1 = max.(0.0,Aᵢ[1:N] .* Rᵢ[1:N])    # Total antibody 
    V2 = max.(0.0,Aᵢ[N+1:end] .* Rᵢ[N+1:end])    # Total antibody
    V̄2 = max.(0.0,Iᵢ[N+1:end] .* Rᵢ[N+1:end] .+ (1 - η) * (Aᵢ[N+1:end] - Iᵢ[N+1:end]))    # Total internal antibody

    # Unquenched
    Q = α₁ * V1 + σ₁ * sqrt.(α₁ * V1) .* randn(N) + E[1][idx1]
    U = α₂ * V1 + σ₂ * sqrt.(α₂ * V1) .* randn(N) + E[2][idx1]

    # Quenched
    Q̄ = α₁ * V̄2 + σ₁ * sqrt.(α₁ * V̄2) .* randn(N) + E[1][idx2]
    Ū = α₂ * V2 + σ₂ * sqrt.(α₂ * V2) .* randn(N) + E[2][idx2]

    return [[Q̄,Ū],[Q,U]]

end
simulate_model(θ::Vector,T::Vector,d::Function,E;kwargs...) = [simulate_model(θ,t,d,E;kwargs...) for t in T]


function simulate_model_noiseless(θ::Vector,T::Vector,d::Function,E=nothing;N::Int=1000,η::Number=1.0)

    # Parameters
    ξᵢ = rand(d(θ[5:end-1]),N) # Sample from parameter distribution
    p  = θ[end]                 # Proportion recycled
    Rᵢ = @view ξᵢ[1,:]
    λᵢ = @view ξᵢ[2,:]
    βᵢ = @view ξᵢ[3,:]
    
    function sample_quenched_quenchable_signal(t)
        Aᵢ = similar(λᵢ)
        Iᵢ = similar(λᵢ)
        for i = 1:N
            Aᵢ[i],Iᵢ[i] = solve_ode_model(t,[λᵢ[i],βᵢ[i],p])
        end
        return Iᵢ .* Rᵢ .+ (1 - η) * (Aᵢ - Iᵢ)
    end

    return [sample_quenched_quenchable_signal(t) for t in T]

end

function parameters_to_antibody(t::Number,λ::Vector,β::Vector,p::Number)
    Aᵢ = similar(λ)
    Iᵢ = similar(λ)
    @simd for i = 1:length(Aᵢ)
        Aᵢ[i],Iᵢ[i] = solve_ode_model(t,[λ[i],β[i],p])
    end
    return Aᵢ,Iᵢ
end
parameters_to_antibody(T::Vector,λ::Vector,β::Vector,p::Number) = [parameters_to_antibody(t,λ,β,p) for t in T]