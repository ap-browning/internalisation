#=
    Figure 7
=#

using Inference
using Model

using JLD2
using KernelDensity
using Plots
using StatsFuns
using StatsPlots
using Random

include("../figure_defaults.jl")

##############################################################
## Load model and results
##############################################################

include("../../Results/MainResult_Setup.jl")
@load "Results/MainResult.jld2"

include("../../Results/DeterministicModel.jl")


##############################################################
## Fig 6ab: Plot statistics at best fit (log scale)
##############################################################

function Wpdf_fraction(θ,t)
    d = param_dist(θ[5:end-1])
    ξ = sample(d,100_000)
    Z = parameters_to_antibody(t,ξ[2,:],ξ[3,:],θ[end])
    W = Z[2] ./ Z[1]
    return x -> pdf(kde(W),x)
end
function Wpdf_amount(θ,t)
    d = param_dist(θ[5:end-1])
    ξ = sample(d,100_000)
    Z = parameters_to_antibody(t,ξ[2,:],ξ[3,:],θ[end])
    W = Z[2]
    return x -> pdf(kde(W),x)
end

# Plot fraction internalised
lt = range(log10(0.01),log10(200.0),length=500)
x  = range(0.0,1.0,length=200)
f  = hcat([Wpdf_fraction(θ,tᵢ)(x) for tᵢ in 10.0.^(lt)]...)

# Normalise pdf by mode
for i = 1:length(lt)
    f[:,i] /= maximum(f[:,i])
end

# Plot amount internalised
lt = range(log10(0.01),log10(200.0),length=500)
x2  = range(0.0,2.2,length=200)
f2  = hcat([Wpdf2(θ,tᵢ)(x2) for tᵢ in 10.0.^(lt)]...)

# Normalise pdf by mode
for i = 1:length(lt)
    f2[:,i] /= maximum(f2[:,i])
end

# Plot
fig6a = contourf(lt,x,f,c=:RdPu_9,lw=0.0,ylabel="Fraction internal")
fig6b = contourf(lt,x2,f2,c=:RdPu_9,lw=0.0,ylabel="Relative amount internalised",ylim=(0.0,2.2))

# Add deterministic solution
x₁ = [solve_ode_model(t,collect(θdet))[2] / solve_ode_model(t,collect(θdet))[1] for t in 10.0.^(lt)]
x₂ = [solve_ode_model(t,collect(θdet))[2]  for t in 10.0.^(lt)]

plot!(fig6a,lt,x₁,lw=3.0,c=:black,label="")
plot!(fig6a,lt,x₁,lw=2.0,c=:white,label="")
plot!(fig6b,lt,x₂,lw=3.0,c=:black,label="")
plot!(fig6b,lt,x₂,lw=2.0,c=:white,label="")

fig6ab = plot(fig6a,fig6b,xticks=(-2:1:2,10.0.^(-2:1:3)),colorbar_title="Probability density")
plot!(fig6ab,size=(800,200),axis=:x,xlabel="Time [min]",widen=:false,grid=:y)

##############################################################
## Fig 6c: Plot with uncertainty
##############################################################

Θ = sample(C[1:100:end],1000)
x = [range(0.0,0.8,length=200); range(0.8,1.0,length=200)]

# Get m,l,u at various time points
t₁ = [1.0;5.0:5.0:30.0]
t₂ = 60.0:30.0:180.0
t  = [t₁;t₂]
L  = Array{Any,1}(undef,length(t))
U  = similar(L)
M  = similar(L)
for (i,tᵢ) in enumerate(t)
    # Calculate lower and upper ci
    l,u = pdf_ci(θ -> Wpdf(θ,tᵢ),x,Θ)
    m = Wpdf(θ,tᵢ)(x)
    # Normalise by mode
    L[i] = l / maximum(m)
    U[i] = u / maximum(m)
    M[i] = m / maximum(m)
end

# Deterministic model
T₁ = range(0.0,35.0,length=200)
T₂ = range(0.0,210.0,length=200)
X₁ = [solve_ode_model(t,collect(θdet))[2] / solve_ode_model(t,collect(θdet))[1] for t in T₁]
X₂ = [solve_ode_model(t,collect(θdet))[2] / solve_ode_model(t,collect(θdet))[1] for t in T₂]

fig6c1 = plot(X₁,T₁,lw=2.0,c=:red,ls=:dash)
for (i,tᵢ) in Iterators.reverse(enumerate(t₁))
    plot!(fig6c1,x,tᵢ .+ 3.33 * U[i],frange=tᵢ .+ 3.33 * L[i],c=RGBA(0.5,0.5,0.5,0.4),lw=0.0)
    plot!(fig6c1,x,tᵢ .+ 3.33 * M[i],c=:black,lw=2.0)
end
plot!(fig6c1,
    yticks = [1; 5:5:30],
    ylim = (0.0,35.0)
)

fig6c2 = plot(X₂,T₂,lw=2.0,c=:red,ls=:dash)
plot!(fig6c2,x,30.0 .+ 15.0 * U[length(t₁)],frange=30.0 .+ 15.0 * L[length(t₁)],c=RGBA(0.5,0.5,0.5,0.4),lw=0.0)
plot!(fig6c2,x,30.0 .+ 15.0 * M[length(t₁)],c=:black,lw=2.0)
for (i,tᵢ) in Iterators.reverse(enumerate(t₂))
    plot!(fig6c2,x,tᵢ .+ 15.0 * U[i+length(t₁)],frange=tᵢ .+ 15.0 * L[i+length(t₁)],c=RGBA(0.5,0.5,0.5,0.4),lw=0.0)
    plot!(fig6c2,x,tᵢ .+ 15.0 * M[i+length(t₁)],c=:black,lw=2.0)
end
plot!(fig6c2,
    yticks = 30:30:180,
    xlim=(0.9,1.0),
    ylim = (0.0,210.0)
)

fig6c = plot(fig6c1,fig6c2,
    xlabel = "Fraction internal antibody",
    ylabel = "Time [min]",
    xrotation = 90,
    yrotation = 90,
    xguidefontrotation = 180,
    ymirror = true,
    axis = :x,
    grid = :x,
    xflip = true,
    size=(400,600),
    bottom_margin=5Plots.mm,
    legend=:none
)

##############################################################
## Fig 6de: Proportion on surface
##############################################################

# Univariate (with uncertainty)
fun = p -> p[3] / sum(p[2:3])
fig6d = plot_distribution(θ -> param_dist(θ[5:end-1]),fun,C,θ,xlabel="Proportion on surface")
vline!(fig6d,[θdet.β / (θdet.λ + θdet.β)],lw=2.0,c=:red,ls=:dash)
plot!(fig6d,
    axis=:x,
    yticks=[],
    grid=:x,
)

# Bivariate (at best fit)
p = rand(param_dist(θ[5:end-1]),100_000)
S = fun.(eachcol(p))
R = p[1,:]
fig6e = density2d(S,R,
            xlabel = "Proportion on surface",
            ylabel = "R",
            fill   = true,
            lw     = 0.0,
            c      = :RdPu_9,
            levels = 8,
            box    = :on,
            ylim = (0.0,3.0),
            xlim = (0.0,1.0)
        )
scatter!(fig6e,[θdet.β / (θdet.λ + θdet.β)],[1.0],c=:white,widen=false,msw=0.0)

# Fig 6de
fig6de = plot(fig6d,fig6e,layout=grid(2,1),size=(300,500),legend=:none,widen=false)
add_plot_labels!(fig6de,offset=3)

##############################################################
## Save...
##############################################################
savefig(fig6ab,"$(@__DIR__)/fig6ab.svg")
savefig(fig6c,"$(@__DIR__)/fig6c.svg")
savefig(fig6de,"$(@__DIR__)/fig6de.svg")