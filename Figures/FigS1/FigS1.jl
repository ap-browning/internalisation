#=
    Figure 3
=#

using Model # To get colours

using CSV
using DataFrames
using DataFramesMeta
using Plots
using StatsPlots

##############################################################
## Data
##############################################################

data = CSV.read("Data/CSV/DualMarker.csv",DataFrame)

# Get autofluorescence information
E = [@subset(data,:Time .== 0.0, :Quenched .== false, :Temp .== 37.0).Signal_Cy5,
     @subset(data,:Time .== 0.0, :Quenched .== false, :Temp .== 37.0).Signal_BDP]

##############################################################
## Figure 3
##############################################################

fig3 = density2d(E[1],E[2],
    colorbar=false,
    xlim=(0.0,20.0),
    ylim=(0.0,20.0),
    box=:on,axis=:all,aspect_ratio=:equal,lw=0.0,c=:RdPu_9,levels=8,fill=true,
    xlabel="Cy5 autofluoresence",
    ylabel="BDP autofluoresence"
)

##############################################################
## Save
##############################################################

savefig("$(@__DIR__)/Fig3.svg")