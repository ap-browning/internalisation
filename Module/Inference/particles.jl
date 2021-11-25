#=

    particles.jl

    Implements the type `Particle` used in ABC SMC.

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=# 

"""
    Particle(θ::Vector,d::Number,w::Number)

SMC particle with location `θ`, distance `d`, and weight `w`.
"""
mutable struct Particle
    θ::Vector   # Location
    d::Number   # Discrepancy
    w::Number   # Weight
end
Particle(θ,d) = Particle(θ,d,1.0)


##############################################################
## Extend methods to particles
##############################################################

import Base: sort!
import Distributions: insupport
import StatsBase: sample, weights
import Statistics: cov, mean


weights(P::Vector{Particle})  = weights([p.w for p in P])
distances(P::Vector{Particle}) = [p.d for p in P]
locations(P::Vector{Particle}) = hcat([p.θ for p in P]...)
insupport(d::Distribution,p::Particle) = insupport(d,p.θ)


"""
    sort!(p::Vector{Particle})

Sort a vector of `Particle`s, `P` by distance `d`.
"""
sort!(P::Vector{Particle}) = permute!(P,sortperm(distances(P)))


"""
    cov(p::Vector{Particle})

Weighted covariance matrix of particles `P`.
"""
cov(P::Vector{Particle}) = cov(permutedims(locations(P)),weights(P))


mean(P::Vector{Particle}) = mean(locations(P),weights(P),dims=2)[:]

"""
    sample(prior,dist,N=1)

Sample particles from the `prior` distribution and attach distance using function `dist`.
"""
function StatsBase.sample(prior::Distribution,dist::Function,w::Number=1)
    θ = rand(prior)
    d = dist(θ)
    Particle(θ,d,w)
end
StatsBase.sample(prior::Distribution,dist::Function,n::Int) = [sample(prior,dist,1/n) for i in 1:n]

"""
    copy(p::Particle)

Copies a particle.
"""
Base.copy(p::Particle) = Particle(p.θ,p.d,p.w)

"""
    findmin(p::Particle)

Returns (pᵢ,i) where pᵢ is the particle with the smallest discrepancy.
"""
function Base.findmin(p::Vector{Particle})
    D = [p[i].d for i = 1:length(p)]
    idx = findmin(D)[2]
    p[idx],idx
end

"""
    minimum(p::Particle)

Finds the particle with the smallest discrepancy.
"""
Base.minimum(p::Vector{Particle}) = findmin(p)[1]

"""
    findmax(p::Particle)

Returns (pᵢ,i) where pᵢ is the particle with the largest discrepancy.
"""
function findmax(p::Vector{Particle})
    D = [p[i].d for i = 1:length(p)]
    idx = findmin(-D)[2]
    p[idx],idx
end

"""
    maximum(p::Particle)

Finds the particle with the largest discrepancy.
"""
Base.maximum(p::Vector{Particle}) = findmax(p)[1]


Base.ndims(P::Vector{Particle}) = ndims(P[1])
Base.ndims(p::Particle) = length(p.θ)

ess(P::Vector{Particle}) = 1 / sum(weights(P).^2)