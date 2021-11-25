#=
    Figure 9
=#

using Inference
using Model

using JLD2
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

##############################################################
## "Noiseless" model: Q̄ at t = 10.0 and 120.0
##############################################################

model2 = θ -> simulate_model_noiseless(θ,[10.0,120.0],param_dist;N=1000,η)

##############################################################
## Figure 7(a)
##############################################################

Z = model2(θ)
ρ = cor([norminvcdf.(u) for u in cdf.(Z)]...)
fig7a = scatter(cdf.(Z)...,
            α      = 0.3,
            msw    = 0.0,
            c      = :black,
            legend = :none,
            box    = :on,
            xlabel = "F_20",
            ylabel = "F_120"
        )
annotate!(fig7a,0.55,0.1,"ρ = $(round(ρ,sigdigits=3))")

##############################################################
## Figure 7(b)
##############################################################

Θ = sample(C[1:100:end],1000)
ρ = [cor([norminvcdf.(u) for u in cdf.(model2(θ))]...) for θ in Θ]
fig7b = stephist(ρ,normalize=:pdf,frange=0.0,
            legend = :none,
            grid   = :x,
            axis   = :x,
            yticks = [],
            xlabel = "cor(Q̄(10),Q̄(120))",
            c      = :grey,
            widen  = false,
            xlim   = (0.49,1.01)
        )

##############################################################
## Figure 7(c)
##############################################################

idx = ρ .> 0.9
Θfilt = Θ[idx]

plts = [(   stephist(hcat(Θ...)[i,:],xlabel=param_names[i],frange=0.0,c=:grey,normalize=:pdf,label="All",α=0.8);
            stephist!(hcat(Θfilt...)[i,:],frange=0.0,c=:purple,normalize=:pdf,α=0.8,label="ρ > 0.9")
) for i = [10,11]]

fig7c = plot(plts...,xlim=(0.0,0.1),layout=grid(2,1),axis=:x,yticks=[],widen=false)

##############################################################
## Figure 7
##############################################################

fig7 = plot(fig7a,fig7b,fig7c,layout=grid(1,3),size=(700,250),bottom_margin=5Plots.mm)
add_plot_labels!(fig7)

##############################################################
## Save...
##############################################################

savefig(fig7,"$(@__DIR__)/fig7.svg")