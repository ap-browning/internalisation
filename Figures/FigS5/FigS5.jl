#=
    Figure Snew
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

##############################################################
## Compute results
##############################################################

# Simulate model at 100 samples from the posterior
Θ = sample(C[1:100:end],100)
Y = [simulate_model(θ,T,param_dist,E;N=100_000,η) for θ in Θ]

##############################################################
## Figure 4(a,b)
##############################################################

#fig4ab = plot_fit_univariate(X,Y,T)
qlim=[0.0,3e4]
ulim=[0.0,150.0]

using KernelDensity

function vertical_density_multi!(plt,X1,X2,lim,i;fill=false,c=:blue)
    x = range(lim...,length=200)
    y1 = pdf(kde(X1),x)
    sc = maximum(y1)
    plot!(plt,i .+ 0.5 * y1/sc,x,fillrange=i,lc=c,fc=Model.brighten_colour(c,0.7),lw=1.0)
    for k = 1:length(X2)
        y2 = pdf(kde(X2[k]),x)
        plot!(plt,i .+ 0.5 * y2/sc,x,c="#333333",α=0.1,lw=1.0)
    end
end

plt1 = plot(legend=:none)
for i = 1:7
    vertical_density_multi!(plt1,X[i][2][1],[y[i][2][1] for y in Y],qlim,i+0.2,c=col_Q)
    vertical_density_multi!(plt1,X[i][1][1],[y[i][1][1] for y in Y],qlim,i,c=col_Q̄)
end
plot!(plt1,xticks=(1:7,T),box=:off,axis=:all,widen=true)

plt2 = plot(legend=:none)
for i = 1:7
    vertical_density_multi!(plt2,X[i][2][2],[y[i][2][2] for y in Y],ulim,i+0.2,c=col_U)
    vertical_density_multi!(plt2,X[i][1][2],[y[i][2][2] for y in Y],ulim,i,c=col_Ū)
end
plot!(plt2,xticks=(1:7,Int.(T)),box=:off,axis=:all,widen=true)

plot!(plt1,
    yformatter=y->y/1e4,
    ylabel="Cy5 (10⁴)",
    xticks=(1:7,[])
)

plot!(plt2,
    ylabel="BDP",
    xlabel="Time [min]"
)

plt = plot(plt1,plt2,layout=grid(2,1),size=(700,350))

plot!(plt,
    axis=:y,
    bottom_margin=2Plots.mm
)

savefig(plt,"$(@__DIR__)/figS5.svg")