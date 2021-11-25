#=
    Figure 1
=#

using Model # To get colours

using CSV
using ColorSchemes
using DataFrames
using DataFramesMeta
using Plots
using StatsPlots

include("../figure_defaults.jl")

##############################################################
## Data (at t = 10 min)
##############################################################

data = @subset(CSV.read("Data/CSV/DualMarker.csv",DataFrame), :Temp .== 37.0, :Time .== 10.0)

##############################################################
## Univariate
##############################################################

p1 = @df @subset(data,(!).(:Quenched)) density(:Signal_Cy5  / 1e4,c=col_Q,lw=2.0,frange=0.0,fα=0.4,label="Not quenched")
@df @subset(data,:Quenched) density!(p1,:Signal_Cy5 / 1e4,c=col_Q̄,lw=2.0,frange=0.0,fα=0.4,ls=:dash,label="Quenched")
plot!(p1,xlim=(0.0,1.5),widen=:false,yticks=[],xlabel="Cy5 (10⁴)")

p2 = @df @subset(data,(!).(:Quenched)) density(:Signal_BDP,c=col_U,lw=2.0,frange=0.0,fα=0.4,label="Not quenched")
@df @subset(data,:Quenched) density!(p2,:Signal_BDP,c=col_Ū,lw=2.0,frange=0.0,fα=0.4,ls=:dash,label="Quenched")
plot!(p2,xlim=(0.0,100.0),widen=:false,yticks=[],xlabel="BDP")

fig1di = plot(p1,p2,layout=grid(2,1),axis=:x,grid=:x)

##############################################################
## Bivariate plots
##############################################################

cs1 = ColorScheme(range(colorant"white",parse(Colorant,col_QU),length=200))
cs2 = ColorScheme(range(colorant"white",parse(Colorant,col_Q̄Ū),length=200))

p3 = @df @subset(data,:Quenched) density2d(:Signal_Cy5,:Signal_BDP,c=palette(cs2),lc=col_Q̄Ū,lw=1.0,fill=true)

p4 = @df @subset(data,(!).(:Quenched)) density2d(:Signal_Cy5,:Signal_BDP,c=palette(cs1),lc=col_QU,lw=1.0,fill=true)

fig1dii = plot(p3,p4,layout=grid(2,1),
    box=:on,xlim=(0.0,1e4),ylim=(0.0,100),
    xformatter=x->x / 1e4,
    xlabel="Cy5 (10⁴)",
    ylabel="BDP",grid=:all)

##############################################################
## Save
##############################################################

plot(fig1di,fig1dii,layout=grid(1,2),size=(500,500))
savefig("$(@__DIR__)/Fig1.svg")