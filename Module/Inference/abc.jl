#=

    abc.jl

    Contains code to conduct approximate Bayesian computation (ABC)

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=# 

##############################################################
## ABC Sequential Monte-Carlo (SMC)
##############################################################
"""
    abc_smc(dist,prior;n=100,α=0.5,S₀=100,c=0.01,ptarget=0.005,maxsteps=20,f=2.0)

Perform ABC SMC.
"""
function abc_smc(dist,prior;n=100,kwargs...)
    P = sample(prior,dist,n)    # Initialise using prior
    print("Initialised $n particles. $(round(minimum(P).d,sigdigits=4)) ≤ d ≤ $(round(maximum(P).d,sigdigits=4)).\n")
    return abc_smc(P,dist,prior;kwargs...)
end
function abc_smc(P,dist,prior;step=1,α=0.5,ptarget=0.01,maxsteps=30,f=2.0)

    # Loop until target reached
    p = 1.0
    while p > ptarget && step ≤ maxsteps
        @printf("Step %-3s (ε = %-8s\b...",
            "$step",
            "$(round(sort(distances(P))[(Int ∘ floor)(α * length(P))],sigdigits=4))",
        )
        
        p = abc_smc_step!(P,dist,prior,α,f)
        step += 1

        @printf("\b\b\b dmin = %-8s p = %-8s ESS = %-8s)\n",
            "$(round(minimum(P).d,sigdigits=4))",
            "$(round(p,sigdigits=2))",
            "$(round(1 / sum(weights(P).^2),sigdigits=5))"
        )     

    end

    # Finish
    return P
end

"""
    abc_smc_step!(P,dist,prior,α,f=2.0)

Perform a single step of ABC SMC, where α * n particles are dicarded.

Proposal is estimated as f * Σ, where Σ is the (weighted) covariance matrix of the non-dicarded particles.

"""
function abc_smc_step!(P,dist,prior,α,f=2.0)
    N = length(P); Nα = (Int ∘ floor)(α * N)

    sort!(P)                # Sort particles
    ε = P[N-Nα].d           # Current threshold   
    Σ = f*cov(P[1:N-Nα])    # Proposal covariance
    prop = MvNormal(Σ)      # Proposal

    # Resample and update particles
    tries = 0
    @threads for j = (N - Nα + 1) : N
        θs,d = similar(P[1].θ),Inf
        while d > ε
            θs = copy(sample(P[1:N-Nα],weights(P[1:N-Nα])).θ) + rand(prop)
            if insupport(prior,θs)
                d = dist(θs)
                tries += 1
            end
        end
        P[j] = Particle(θs,d,pdf(prior,θs) / sum(p.w * pdf(prop,θs - p.θ) for p in P[1:N-Nα]))
    end

    # Reweight
    w = weights(P)
    w[1:N-Nα] /= sum(w[1:N-Nα])         # Normalise
    w[N-Nα+1:N] /= sum(w[N-Nα+1:N])     # Normalise
    w[1:N-Nα] *= (N - Nα) / N           # Re-weight
    w[N-Nα+1:N] *= Nα / N               # Re-weight
    [p.w = wᵢ for (p,wᵢ) in zip(P,w)]

    # Report acceptance probability (exclude trials where we propose a value outside the prior)
    return Nα / tries
end


##############################################################
## ABC Markov-chain Monte-Carlo (MCMC)
##############################################################
"""
    abc_mcmc(dist,prior,θ₀,Σ,ε;n=1000,nchains=nthreads(),...)

Perform ABC MCMC using the Metropolis-Hastings algorithm (implemented through `AdvancedMH.jl`).
"""
function abc_mcmc(dist,prior,θ₀,Σ,K::Function;n=1000,nchains=Threads.nthreads(),kwargs...)
    parallel_args = nchains == 1 ? n : (MCMCThreads(),n,nchains)
    sample(
        DensityModel(θ -> insupport(prior,θ) ? log(K(dist(θ))) : -Inf),
        RWMH(MvNormal(Σ)),
        parallel_args...;
        init_params=θ₀,
        chain_type=Chains,
        kwargs...
    )
end
abc_mcmc(dist,prior,θ₀,Σ,ε::Number;kwargs...) = abc_mcmc(dist,prior,θ₀,Σ,d -> d ≤ ε ? 1.0 : 0.0;kwargs...)


##############################################################
## ABC MCMC Discrepancy Metrics
##############################################################

logistic(d,ε,k) = 1 / (1 + exp(k * (d - ε)))
logit(p,ε,k) = ε + log((1 - p) / p) / k

##############################################################
## MCMCChains.Chains helper functions
##############################################################

# Convert MCMC chain output to a matrix
Base.Matrix(C::Chains) = permutedims(hcat([Matrix(C[:,i,:])[:] for i = 1:size(C,2)-1]...))

acceptance_rate(C::Chains) = mean(hcat([cumsum(C[1:end-1,1,i] .!= C[2:end,1,i],dims=1)[:] ./ (1:size(C,1)-1) for i = 1:size(C,3)]...),dims=2)[:]

MCMCChains.missing_datetime(T::Type) = missing

StatsBase.sample(C::Chains,n) = (CM = Matrix(C); [CM[:,i] for i in sample(1:size(CM,2),n)])
StatsBase.sample(C::Chains) = sample(C,1)[1]