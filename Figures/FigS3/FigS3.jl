#=
    Figure S3
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
## Data
##############################################################

data = @subset(CSV.read("Data/CSV/DualMarker.csv",DataFrame), :Temp .== 37.0, :Time .> 0.0)

##############################################################
## Univariate
##############################################################

T = sort(unique(data.Time))

plts_q = [plot() for t in T]
plts_u = [plot() for t in T]

for (i,t) in enumerate(T)

    @df @subset(data,:Time .== t,(!).(:Quenched)) density!(plts_q[i],:Signal_Cy5 / 1e4,c=col_Q,lw=2.0,frange=0.0,fα=0.4,label="Not quenched")
    @df @subset(data,:Time .== t,:Quenched) density!(plts_q[i],:Signal_Cy5 / 1e4,c=col_Q̄,lw=2.0,frange=0.0,fα=0.4,ls=:dash,label="Quenched")
    plot!(plts_q[i],xlim=(0.0,3.0),widen=:false,yticks=[],xlabel="Cy5 (10⁴)",axis=:x,title="$(Int(t)) min")

    @df @subset(data,:Time .== t,(!).(:Quenched)) density!(plts_u[i],:Signal_BDP,c=col_U,lw=2.0,frange=0.0,fα=0.4,label="Not quenched")
    @df @subset(data,:Time .== t,:Quenched) density!(plts_u[i],:Signal_BDP,c=col_Ū,lw=2.0,frange=0.0,fα=0.4,ls=:dash,label="Quenched")
    plot!(plts_u[i],xlim=(0.0,120.0),widen=:false,yticks=[],xlabel="BDP",axis=:x)

    if i > 1
        plot!(plts_q[i],legend=:none)
        plot!(plts_u[i],legend=:none)
    end

end

figS3 = plot(plts_q...,plts_u...,layout=grid(2,7),size=(800,250))

##############################################################
## Save...
##############################################################

savefig("$(@__DIR__)/FigS3.svg")