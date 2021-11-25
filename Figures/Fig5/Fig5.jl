#=
    Figure 6
=#

using Inference
using Model

using JLD2
using Plots
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
## Figure 5
##############################################################

# Parameter limits
lims = [[0.0,3.0],[0.0,0.5],[0.0,0.2]]

# Dimension combinations for bivariate plots
dims = [[1,2],[1,3],[2,3]]

plts1 = [plot_distribution(θ -> param_dist(θ[5:end-1]).v[i],sample(C,2000),θ;lims=lims[i],xlabel=["R","λ","β"][i]) for i = 1:3]
vline!(plts1[1],[1.0],lw=2.0,c=:red,ls=:dash)
vline!(plts1[2],[θdet.λ],lw=2.0,c=:red,ls=:dash)
vline!(plts1[3],[θdet.β],lw=2.0,c=:red,ls=:dash)

plts2 = [plot_biv_distribution(param_dist(θ[5:end-1]);
            dims=dims[i],
            lims=lims[dims[i]],
            labs=["R","λ","β"][dims[i]],
            c=:RdPu_9,levels=10,box=:on
        ) for i = 1:3]
scatter!(plts2[1],[1.0],[θdet.λ],c=:white,msw=0.0)
scatter!(plts2[2],[1.0],[θdet.β],c=:white,msw=0.0)
scatter!(plts2[3],[θdet.λ],[θdet.β],c=:white,msw=0.0)

fig5 = plot(plts1[1],plts2[1],plts1[2],plts2[2:3]...,plts1[3],layout=@layout([a _ _; b c _; d e f]))

plot!(fig5,
    widen=false,
    legend=:none,
    size=(600,600)
)
[plot!(fig5,yticks=[],axis=:x,subplot=i) for i = [1,3,6]]
plot!(fig5,yticks=[],subplot=5)

##############################################################
## Save...
##############################################################

savefig(fig5,"$(@__DIR__)/fig5.svg")