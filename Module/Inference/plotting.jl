#=

    plotting.jl

    Contains code to produce inference related plots.

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=#

"""
    plot_posteriors(Θ)

Produce density plot of posteriors. Optionally provide names of variables (a vector of `String`s or `Symbol`s) and 
`Matrix{Symbol}` to specify the plot porder and porder.

"""
function plot_posteriors(Θ::Matrix,W=ones(size(Θ,2));
        idx = trues(size(Θ,2)),
        param_names = ["p$i" for i = 1:size(Θ,1)],
        porder = nothing,
        lims = nothing,
        nbins = nothing,
        layout=nothing,
        kwargs...
    )

    n = length(param_names)
    param_names = Symbol.(param_names)
    porderv = porder === nothing ? param_names : filter(x -> x != :_,permutedims(porder)[:])

    # Convert to named tuple for plotting
    Θ2 = NamedTuple{Tuple(param_names)}([Θ[i,:] for i = 1:n])
    lims2 = lims === nothing ? nothing : NamedTuple{Tuple(param_names)}(lims)

    # Loop through porder and create plots
    plts = [stephist(Θ2[porderv[i]][idx],weights=W[idx],xlabel=porderv[i],
        xlim = lims === nothing ? :auto : lims2[porderv[i]],
        frange = 0.0,
        c = :grey,
        bins = (lims === nothing) || (nbins === nothing) ? :auto : range(lims2[porderv[i]]...,length=nbins),
        kwargs...
    ) for i = 1:length(porderv)]

    # Put together
    Plots.plot(plts...,
        layout=layout === nothing ? length(porderv) : layout,
        xtickfontsize=6,
        yticks=[],
        size=(700,400),
        legend=:none,
        widen=false
    )

end

plot_posteriors(P::Vector{Particle};kwargs...) = plot_posteriors(locations(P),weights(P);kwargs...)
plot_posteriors(C::Chains;burnin=0,skip=1,kwargs...) = plot_posteriors(Matrix(C[burnin+1:skip:end]);param_names=C.name_map.parameters,kwargs...)


function plot_chains(C::Chains;
        porder = nothing,
        lims = nothing
    )  

    n = size(C,2) - 1
    param_names = C.name_map.parameters
    porderv = porder === nothing ? param_names : permutedims(porder)[:]

    # Convert to named tuple for plotting
    lims2 = lims === nothing ? nothing : NamedTuple{Tuple(param_names)}(lims)

    # Loop through porder and create plots
    plts = [plot(C[:,porderv[i],1],xlabel=porderv[i],
        ylim = lims === nothing ? :auto : lims2[porderv[i]]) for i = 1:n]

    # Put together
    Plots.plot(plts...,
        porder=porder===nothing ? length(param_names) : grid(size(porder)...),
        xtickfontsize=6,
        ytickfontsize=6,
        size=(700,400),
        legend=:none
    )

end

