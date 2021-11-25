#=

    discrepancies.jl

    Contains code to assess discrepancy between data and model simulation using Cramér–von Mises distance.

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=# 

"""
    makecdf(x::Vector;n=1000)

Creates a linear interpolation approximation of the distribution function of `x`.

"""
function makecdf(x::Vector;n=1000)
    n_x = length(x)
    xp = range(minimum(x),maximum(x) + 1e-3,length=n)
    p = [count(x .< xpᵢ) for xpᵢ in xp] / n_x
    itp = LinearInterpolation(xp, p, extrapolation_bc = Interpolations.Flat())
    return itp
end


"""
    cvmises(cdf::AbstractInterpolation,y::Vector)

Computes Cramér–von Mises distance between distribution defined by the distribution function `cdf`
(given as an interpolation) and the univariate sample `y`.

"""
function cvmises(cdf::AbstractInterpolation,y::Vector)
    n = length(y)
    cdfs = cdf(sort(y))
    return sqrt(sum( ((2(1:n) .- 1.0) / 2n - cdfs).^2 ) / n + 1 / 12n^2)
end


"""
    addist(cdf::AbstractInterpolation,y::Vector)

Computes Anderson-Darling distance between distribution defined by the distribution function `cdf`
(given as an interpolation) and the univariate sample `y`.

"""
function addist(cdf::AbstractInterpolation,y::Vector)
    n = length(y)
    cdfs = min.(max.(cdf(sort(y)),1e-9),1-1e-9)
    sqrt(-n - sum( (2(1:n) .- 1.0) / n .* (log.(cdfs) + log.(1.0 .- cdfs[end:-1:1]) )  ))
end

"""
    ksdist(cdf::AbstractInterpolation,y::Vector)

Computes Kolmogorov–Smirnov distance between distribution defined by the distribution function `cdf`
(given as an interpolation) and the univariate sample `y`.

"""
function ksdist(cdf::AbstractInterpolation,y::Vector)
    n = length(y)
    cdfs = cdf(sort(y))
    δp = maximum((1:n) / n - cdfs)
    δn = -minimum((0:n-1) / n - cdfs)
    return max(δn, δp)
end


"""
kuiperdist(cdf::AbstractInterpolation,y::Vector)

Computes Kuiper's distance between distribution defined by the distribution function `cdf`
(given as an interpolation) and the univariate sample `y`.

"""
function kuiperdist(cdf::AbstractInterpolation,y::Vector)
    n = length(y)
    cdfs = cdf(sort(y))
    δp = maximum((1:n) / n - cdfs)
    δn = -minimum((0:n-1) / n - cdfs)
    return δn + δp
end


"""
    makediscrepancy(X;w=[0.5,0.5])

Create discrepancy function `dist(Y)` using data `X` where:
    X[i] = [[Q̄,Ū],[Q,U]]
and
    dist(Y) = ΣᵢΣⱼ( 2w[2] * cor(X[i][j]...) + w[1] * Σₖ addist(X[i][j][k]) )

By default, `w = [0.5,0.5]` so Anderson-Darling distance is weighted equally with correlation.

"""
function makediscrepancy(X::Vector{Vector{Vector{Vector{Float64}}}},f=addist;w=[0.5,0.5])
    Xcdf = [[makecdf.(X[i][j]) for j = 1:length(X[1])] for i = 1:length(X)]
    Xcor = [[cor(X[i][j]...) for j = 1:length(X[1])] for i = 1:length(X)]
    function dist(Y::Vector{Vector{Vector{Vector{Float64}}}})
        d = 0.0
        for i = 1:length(X), j = 1:length(X[1])
            d += 2 * w[1] * f(Xcdf[i][j][1],Y[i][j][1]) +    # Quenchable
                 1 * w[1] * f(Xcdf[i][j][2],Y[i][j][2]) +    # Unquenchable
                 1 * w[2] * norm(Xcor[i][j] - cor(Y[i][j]...))  # Correlations
        end
        return d
    end
    return dist
end


"""
    makediscrepancy_univariate(X;w=[0.5,0.5])

Create discrepancy function `dist(Y)` using data `X` where:
    X[i] = [[Q̄,Ū],[Q,U]]
and
    dist(Y) = ΣᵢΣⱼ addist(X[i][j][1])

Compared to `makediscrepancy`, only information relating to the quenchable probe is considered.

"""
function makediscrepancy_univariate(X::Vector{Vector{Vector{Vector{Float64}}}},f=addist)
    Xcdf = [[makecdf.(X[i][j]) for j = 1:length(X[1])] for i = 1:length(X)]
    function dist(Y::Vector{Vector{Vector{Vector{Float64}}}})
        d = 0.0
        for i = 1:length(X), j = 1:length(X[1])
            d += w[1] * f(Xcdf[i][j][1],Y[i][j][1])   # Quenchable
        end
        return d
    end
    return dist
end
