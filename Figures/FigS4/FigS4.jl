#=
    Figure S4
=#
using Optim

xlim = [0.0,0.5]
xgrd = range(xlim...,length=200)

μ = 0.25
σ = [0.01,0.05,0.1]

ω = [-1.5,-0.5,0.0001,0.5,1.5]

σᵢ = σ[1]

plts = [plot() for i = 1:length(σ)]

for (i,σᵢ) in enumerate(σ)

    pltsᵢ = [plot() for i = 1:length(ω)]

    for (j,ωⱼ) in enumerate(ω)

        λ = GammaAlt(μ,σᵢ,ωⱼ)

        # Find mode
        m = optimize(x -> -pdf(λ,x),xlim...).minimizer

        if ωⱼ > 0.0
            xmap = x -> x .> m ? 2m - x : x
        elseif ωⱼ < 0.0
            xmap = x -> x .< m ? 2m - x : x
        else
            xmap = x -> x
        end

        # Plot symmetric distribution
        plot!(pltsᵢ[j],xgrd,pdf(λ,xmap.(xgrd)),lw=2.0,c=:red,α=0.5,ls=:dash)
        plot!(pltsᵢ[j],xgrd,pdf(λ,xgrd),lw=2.0,c=:black)

    end

    plts[i] = plot(pltsᵢ...,layout=grid(1,length(ω)))
    if i < length(σ)
        plot!(plts[i],xticks=(0.0:0.1:0.5,[]))
    else
        plot!(plts[i],xticks=0.0:0.1:0.5)
    end
    plot!(plts[i],subplot=1,ylabel="σ = $σᵢ")

end

figS4 = plot(plts...,layout=grid(length(σ),1),yticks=[],axis=:x,legend=:none,size=(700,400))
[plot!(figS4,subplot=i,title="ω = $ωᵢ") for (i,ωᵢ) in enumerate(ω)]
[plot!(figS4,subplot=i,xlabel="λ") for i = 11:15]
figS4

##############################################################
## Save...
##############################################################

savefig("$(@__DIR__)/FigS4.svg")