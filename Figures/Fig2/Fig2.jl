#=
    Figure 2
=#

using Model

using CSV
using DataFrames
using DataFramesMeta
using Plots
using StatsPlots

include("../figure_defaults.jl")

##############################################################
## Run/load analysis
##############################################################

include("../../Results/DeterministicModel.jl")

##############################################################
## Fig 2b (model v data)
##############################################################

fig2b = @df @subset(data_summary,(!).(:Quenched)) scatter(:Time,:Signal_Cy5, 
                msc=col_Q,c=col_Q,label="Data") 
fig2b = @df @subset(data_summary,:Quenched) scatter!(fig2b,:Time,:Signal_Cy5,
                msc=col_Q̄,c=col_Q̄,label="Data (Quenched)") 

# Plot solution
t = range(0.0,180.0,length=200)
X = hcat([solve_ode_model(tᵢ,collect(θdet)) for tᵢ in t]...)
m_unquenched = X[1,:] .+ e / θdet.α
m_quenched = X[2,:] + (1 - η) * (X[1,:] - X[1,:]) .+ e / θdet.α
plot!(fig2b,t,θdet.α * m_unquenched,c=col_Q,label="Model")
plot!(fig2b,t,θdet.α * m_quenched,c=col_Q̄,label="Model (Quenched)")
plot!(fig2b,
    legend=:bottomright,
    ylim=(0.0,1.2e4),
    yformatter=y->y/1e4,
    ylabel="Cy5 (10⁴)",
    yticks=0:2e3:1.2e4
)

##############################################################
## Fig 2c (model all variables)
##############################################################

fig2c = plot()

X = hcat([solve_ode_model_all_vars(tᵢ,collect(θdet)) for tᵢ in t]...) 
cols = [:gray,:blue,:red,:orange]
labs = ["T(t)","S(t)","E(t)","F(t)"]
stys = [:dash,:solid,:solid,:solid]
[plot!(fig2c,t,X[i,:],c=cols[i],label=labs[i],ls=stys[i]) for i = 1:4]
plot!(fig2c,
    ylim=(0.0,1.0),
    ylabel="Relative count"
)

##############################################################
## Fig 2
##############################################################

# Figure 2
fig2 = plot(fig2b,fig2c,
    xticks=0:30:180,
    xlim=(-5.0,185.0),
    xlabel="Time [min]",
    box=:on,axis=:all,widen=true,lw=2.0,
    size=(800,300)
)
add_plot_labels!(fig2)

##############################################################
## Save
##############################################################

savefig("$(@__DIR__)/Fig2.svg")