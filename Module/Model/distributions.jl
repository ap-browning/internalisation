#=

    distributions.jl

    Contains code to create alternatively parameterised distributions, including copulae.

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=#

import Base: rand, minimum, maximum
import Distributions: pdf, cdf
import Statistics: mean, var, std, quantile
import StatsBase: skewness, sample

abstract type Copula end

##############################################################
## Alternative parameterisation of the gamma distribution
##############################################################
struct GammaAlt{T<:Real} <: ContinuousUnivariateDistribution
    μ::T
    σ::T
    ω::T
    d::ContinuousUnivariateDistribution  # Underlying (truncated) Gamma distribution
    GammaAlt{T}(μ,σ,ω,d) where {T} = new{T}(μ,σ,ω,d) 
end
struct GammaAltNegative{T<:Real} <: ContinuousUnivariateDistribution
    μ::T
    σ::T
    ω::T
    d::ContinuousUnivariateDistribution  # Underlying (truncated) Gamma distribution
    GammaAltNegative{T}(μ,σ,ω,d) where {T} = new{T}(μ,σ,ω,d) 
end

"""
    GammaAlt(μ,σ,ω)

Construct a truncated (x > 0) Gamma distribution `d` where the mean, standard deviation
and skewness of the untruncated distribution are given by μ, σ and ω.

"""
function GammaAlt(μ::T,σ::T,ω::T) where {T <: Real}
    μ > 0.0 || error("Mean must be positive.")
    σ > 0.0 || error("Standard deviation must be positive.")
    if ω < 0
        α = 4/ω^2
        θ = -σ * ω / 2
        d = Truncated(Gamma(α, θ) - α*θ - μ,-Inf,0.0)
        GammaAltNegative{T}(μ,σ,ω,d)
    else
        α = 4/ω^2
        θ = σ * ω / 2
        d = Truncated(Gamma(α, θ) - α*θ + μ,0.0,Inf)
        GammaAlt{T}(μ,σ,ω,d)
    end
end

#### Evaluation
rand(rng::AbstractRNG, d::GammaAltNegative) = -rand(rng,d.d)
pdf(d::GammaAltNegative,x::Real) = pdf(d.d,-x)
logpdf(d::GammaAltNegative,x::Real) = logpdf(d.d,-x)
cdf(d::GammaAltNegative,x::Real) = 1 .- cdf(d.d,-x)
quantile(d::GammaAltNegative,p::AbstractArray) = -quantile(d.d,1.0 .- p)
quantile(d::GammaAltNegative,p::Number) = -quantile(d.d,1.0 .- p)
minimum(d::GammaAltNegative) = 0.0
maximum(d::GammaAltNegative) = Inf

rand(rng::AbstractRNG, d::GammaAlt) = rand(rng,d.d)
pdf(d::GammaAlt,x::Real) = pdf(d.d,x)
logpdf(d::GammaAlt,x::Real) = logpdf(d.d,x)
cdf(d::GammaAlt,x::Real) = cdf(d.d,x)
quantile(d::GammaAlt,p::AbstractArray) = quantile(d.d,p)
quantile(d::GammaAlt,p::Number) = quantile(d.d,p)
minimum(d::GammaAlt) = 0.0
maximum(d::GammaAlt) = Inf


##############################################################
## Alternative parameterisation of the Log Normal distribution
##############################################################
"""
    LogNormalAlt(μ,μ₁,σ₁)

Construct a truncated (x > 0) Log-Normal distribution `d` where the mean of the 
untruncated distribution is `μ`. `μ₁` and `σ₁` are the standard Log Normal parameters
in the unshifted distribution.
"""
LogNormalAlt(μ,μ₁,σ₁) = Truncated(LogNormal(μ₁,σ₁) - exp(μ₁ + σ₁^2/2) + μ,0,Inf)


##############################################################
## Gaussian copula
##############################################################
"""
    GaussianCopula(P)

Construct a Gaussian copula with correlation matrix `P`.
"""
struct GaussianCopula{N} <: Copula
    P::Matrix           # Covariance matrix
    L::LowerTriangular  # Cholesky decomposition of `P`
end
GaussianCopula(P::Matrix) = GaussianCopula{size(P,1)}(P,cholesky(P).U')
GaussianCopula(p) = (P = vec_to_cor(p); GaussianCopula{size(P,1)}(P,cholesky(P).U'))

rand(rng::AbstractRNG,C::GaussianCopula{N}) where {N} = normcdf.(C.L * randn(rng,N))
rand(rng::AbstractRNG,C::GaussianCopula{N},n::Int) where {N} = normcdf.(C.L * randn(rng,N,n))
rand(C::GaussianCopula{N},n::Int) where {N} = normcdf.(C.L * randn(N,n))
sample(C::GaussianCopula{N},n::Int=1) where {N} = normcdf.(C.L * randn(N,n))

pdf(C::GaussianCopula{N},u::Vector) where {N} = (x = norminvcdf.(u); 1 / sqrt(det(C.P)) * exp.(-0.5 * x' * (inv(C.P) - I) * x))

marginalize(C::GaussianCopula{N};dims=[1,2]) where {N} = GaussianCopula(C.P[dims,dims])


##############################################################
## t copula
##############################################################
"""
    TCopula(P,ν=1)

Construct a t-copula with correlation matrix `P` and degrees of freedom `ν`. 

WARNING: CURRENTLY UNTESTED

"""
struct TCopula{N} <: Copula
    P::Matrix           # Covariance matrix
    ν::Real             # df parameter
    L::LowerTriangular  # Cholesky decomposition of `P`
end
TCopula(P::Matrix,ν::Real=1) = TCopula{size(P,1)}(P,ν,cholesky(P).U')
TCopula(p,ν::Real=1) = (P = vec_to_cor(p); TCopula{size(P,1)}(P,ν,cholesky(P).U'))

rand(rng::AbstractRNG,C::TCopula{N}) where {N} = tdistcdf.(C.ν,C.L * randn(rng,N) / sqrt.(rand(Chisq(C.ν))' / C.ν))
rand(rng::AbstractRNG,C::TCopula{N},n::Int) where {N} = tdistcdf.(C.ν,C.L * randn(rng,N,n) ./ sqrt.(rand(Chisq(C.ν),n)' / C.ν))
rand(C::TCopula{N},n::Int) where {N} = tdistcdf.(C.ν,C.L * randn(N,n) ./ sqrt.(rand(Chisq(C.ν),n)' / C.ν))
sample(C::TCopula{N},n::Int=1) where {N} = tdistcdf.(C.ν,C.L * randn(N,n) ./ sqrt.(rand(Chisq(C.ν),n)' / C.ν))

Γ = gamma
function pdf(C::TCopula{N},u::Vector) where {N}
    x = tdistinvcdf.(C.ν,u);
    Γ((C.ν + N)/2) / 
    (Γ(C.ν/2) * C.ν^(N/2) * π^(N/2) * sqrt(det(C.P))) * 
    (1 + x' * inv(C.P) * x / C.ν)^(-(C.ν + N) / 2) / 
    prod(tdistpdf.(C.ν,x))
end

##############################################################
## Copula-MVDistribution
##############################################################
"""
    MvDependent(C,v)

Construct dependent multivariate distribution where the marginals are described by
elements of `v` (a vector), and the dependence structure is described a copula `C`.
"""
struct MvDependent{N}
    C::Copula
    v::Vector
    q::Vector{Bool}
end
MvDependent(C::Copula,v::Vector;q=falses(length(v))) = MvDependent{length(v)}(C,v,q)

function sample(d::MvDependent{N},n=1) where {N}
    U = rand(d.C,n)
    for i = 1:N
        U[i,:] = d.q[i] ? quantile_interpolated(d.v[i],U[i,:]) : quantile(d.v[i],U[i,:])
    end
    return U
end
rand(d::MvDependent,n::Int=1) = sample(d,n)

marginalize(d::MvDependent;dims=[1,2]) = MvDependent(
    marginalize(d.C;dims),
    d.v[dims],
    q=d.q[dims]
)

pdf(d::MvDependent,x::Vector) = pdf(d.C,max.(0.0,cdf.(d.v,x))) * prod(pdf.(d.v,x))

##############################################################
## Helpful functions
##############################################################
"""
    vec_to_cor(p)

Convert number of vector `p` to a correlation matrix.

Example:

    vec_to_cor([0.1,0.2,0.3])
    3×3 Matrix{Float64}:
     1.0  0.1  0.2
     0.1  1.0  0.3
     0.2  0.3  1.0
"""
function vec_to_cor(p::Vector)
    n = 0.5 * (1 + sqrt(1 + 8length(p)))
    isinteger(n) || error("Check length of input.")
    n = Int(n)
    P = zeros(n,n)
    P[(triu(ones(n,n)) - I) .== 1.0] = p
    P += P' + I
end
vec_to_cor(p::Real) = vec_to_cor([p])


"""
    quantile_interpolated(d,u)

Approximate quantile function using interpolation. Output is comparable to 
`quantile(d,u)`, with a considerable speed-up if `length(u)` is large and 
the quantile function of `d` is expensive.
"""
function quantile_interpolated(d::Distribution,u::Vector;pts::Int=100)
    length(u) > 2 || "`u` must contain more than 2 elements."
    x = range(quantile(d,[extrema(u)...])...,length=pts)
    y = cdf(d,x)
    LinearInterpolation(y,x,extrapolation_bc=Line()).(u)
end


Distributions.cdf(x::Vector) = invperm(sortperm(x)) / length(x) .- 1 / 2length(x)
