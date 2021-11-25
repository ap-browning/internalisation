#=
    Figure 5
=#

using Inference
using Model

using ColorSchemes
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
## Posteriors
##############################################################

porder = [:μλ :σλ  :ωλ  :μβ  :σβ  :ωβ;
          :μR  :σR  :p :ρRλ :ρRβ :ρ̄λβ ;
          :α₁  :α₂  :σ₁  :σ₂ :σ₂ :σ₂]   # Plot σ₂ thrice to make grid even

layout = grid(3,6)

fig3 = plot_posteriors(C,lims=extrema.(prior.v);porder,layout,nbins=40)

plot!(fig3,subplot=13,xformatter=y->Int(y/1000),xlabel="α₁ ('000)")
plot!(fig3,subplot=2,xticks=0.0:0.1:0.3)
add_plot_labels!(fig3,offset=2)
plot!(fig3,
    grid=:none,
    axis=:x,
    size=(800,400)
)

##############################################################
## Inset (joint distribution of σ₁ and σ₂)
##############################################################

cs = ColorScheme(range(colorant"white",colorant"black",length=200))

fig3s = histogram2d(C[:,3,:][:],C[:,4,:][:],
            c=palette(cs),
            bins=20,
            xlabel="σ₁",
            ylabel="σ₂",
            box=:on,
            xlim=(0.0,10.0),
            ylim=(0.0,0.8),
            clim=(0.0,4500),
            size=(300,300))

##############################################################
## ABC demonstration
##############################################################

# Choose a "poor" parameter combination
θpoor = 0.5*θ + 0.5*rand(prior)

# Simulate model at best fit, and at a poor parameter combination
Y = model(θ)
Ypoor = model(θpoor)

# Data and model at t = 30 min
t = 30.0
x = X[findfirst(T .== 10.0)][1]
y = Y[findfirst(T .== 10.0)][1]
ypoor = Ypoor[findfirst(T .== 10.0)][1]

# Cy5 plots
fig3_p1_a1 = density(x[1],lw=2.0,c=:black,label="Data")
density!(fig3_p1_a1,y[1],lw=2.0,c=:red,label="Model")

fig3_p1_a2 = density(x[1],lw=2.0,c=:black,label="Data")
density!(fig3_p1_a2,ypoor[1],lw=2.0,c=:red,label="Model")

fig3_p1_a = plot(fig3_p1_a1,fig3_p1_a2,layout=grid(2,1),link=:x,axis=:x,yticks=[],xlim=(0.0,2e4),xformatter=x->x/1e4,xlabel="Cy5 (10⁴)")

# BDP plots
fig3_p1_b1 = density(x[2],lw=2.0,c=:black,label="Data")
density!(fig3_p1_b1,y[2],lw=2.0,c=:red,label="Model")

fig3_p1_b2 = density(x[2],lw=2.0,c=:black,label="Data")
density!(fig3_p1_b2,ypoor[2],lw=2.0,c=:red,label="Model")

fig3_p1_b = plot(fig3_p1_b1,fig3_p1_b2,layout=grid(2,1),link=:x,axis=:x,yticks=[],xlim=(0.0,150),xlabel="BDP")

# Joint distribution plots
cs1 = ColorScheme(range(colorant"white",colorant"black",length=200))

fig3_p1_c1 = density2d(x...,c=palette(cs1),fill=true,lw=0.0)
density2d!(fig3_p1_c1,y...,c=:red)

fig3_p1_c2 = density2d(x...,c=palette(cs1),fill=true,lw=0.0)
density2d!(ypoor...,c=:red,levels=30)

fig3_p1_c = plot(fig3_p1_c1,fig3_p1_c2,layout=grid(2,1),link=:all,xlim=(0.0,1.5e4),ylabel="BDP",ylim=(0.0,100.0),box=:on,xformatter=x->x/1e4,xlabel="Cy5 (10⁴)")

# ABC demo plot
fig3ab = plot(fig3_p1_a,fig3_p1_b,fig3_p1_c,size=(500,300),layout=grid(1,3))

##############################################################
## Save...
##############################################################

savefig(fig3,"$(@__DIR__)/fig3.svg")
savefig(fig3s,"$(@__DIR__)/fig3s.svg")
savefig(fig3ab,"$(@__DIR__)/fig3ab.svg")
