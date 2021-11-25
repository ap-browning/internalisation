#=

    plotting.jl

    Contains code to produce model related plots.

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=# 
"""
    density2d(x,y,...)

Create 2D kernel density plot.
"""
@userplot density2d
@recipe function f(kc::density2d; levels=10, clip=((-3.0, 3.0), (-3.0, 3.0)))
    x,y = kc.args

    x = vec(x)
    y = vec(y)

    k = KernelDensity.kde((x, y))

    legend --> false

    @series begin
        seriestype := contourf
        colorbar := false
        (collect(k.x), collect(k.y), k.density')
    end

end

cs1 = ColorScheme(range(colorant"white",parse(Colorant,col_QU),length=200))
cs2 = ColorScheme(range(colorant"white",parse(Colorant,col_Q̄Ū),length=200))

function plot_fit_univariate(X,Y,T;qlim=[0.0,3e4],ulim=[0.0,150.0])

    plt1 = plot(legend=:none)
    for i = 1:7
        vertical_density!(plt1,X[i][2][1],Y[i][2][1],qlim,i+0.2,c=col_Q)
        vertical_density!(plt1,X[i][1][1],Y[i][1][1],qlim,i,c=col_Q̄)
    end
    plot!(plt1,xticks=(1:7,T),box=:off,axis=:all,widen=true)

    plt2 = plot(legend=:none)
    for i = 1:7
        vertical_density!(plt2,X[i][2][2],Y[i][2][2],ulim,i+0.2,c=col_U)
        vertical_density!(plt2,X[i][1][2],Y[i][1][2],ulim,i,c=col_Ū)
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

end




function plot_fit_bivariate(X,Y,T;qlim=[0.0,2.4e4],ulim=[0.0,120.0],levels=5,kwargs...)

    plts_Q = [plot() for i = 1:7]
    plts_Q̄ = [plot() for i = 1:7]
    for i = 1:7
        density2d!(plts_Q[i],X[i][2]...,xlim=qlim,ylim=ulim,c=palette(cs1),lw=0.0,fill=true,levels=levels)
        density2d!(plts_Q[i],Y[i][2]...,xlim=qlim,ylim=ulim,c="#333333",lw=0.5,levels=levels)
        density2d!(plts_Q̄[i],X[i][1]...,xlim=qlim,ylim=ulim,c=palette(cs2),lw=0.0,fill=true,levels=levels)
        density2d!(plts_Q̄[i],Y[i][1]...,xlim=qlim,ylim=ulim,c="#333333",lw=0.5,levels=levels)
    end
    
    plts3 = plot(plts_Q...,box=:on,axis=:on,xticks=(0:0.8e4:2.4e4,[]),yticks=(0.0:40.0:120.0,[]),layout=grid(1,7),margin=-1Plots.mm)
    plts4 = plot(plts_Q̄...,box=:on,axis=:on,xticks=(0:0.8e4:2.4e4,[]),yticks=(0.0:40.0:120.0,[]),layout=grid(1,7),margin=-1Plots.mm)
    
    
    plot(plts3,plts4;layout=grid(2,1),size=(700,200),bottom_margin=-2Plots.mm,left_margin=-3Plots.mm,kwargs...)
    
end


##############################################################
## R, λ, β uncertainty plots
##############################################################

"""
    plot_distribution(d,Θ::Vector{Vector})
    plot_distribution(d,C::Chains)

Plots the distribution `d(θ)` with 95% credible interval.
"""
function plot_distribution(
    d::Function,
    Θ::Vector{Vector{Float64}},
    θ=nothing;
    lims = quantile(d(mean(Θ)),[0.0001,0.9999]),
    n_x = 200,
    kwargs...
)

    x = range(lims...,length=n_x)
    f = hcat([pdf(d(θ),x) for θ in Θ]...)
    l = [quantile(f[i,:],0.025) for i in 1:n_x]
    u = [quantile(f[i,:],0.975) for i in 1:n_x]

    plt = plot(x,u,frange=l,c=:black,α=0.2,lw=0.0,label="")
    if θ === nothing && (!).(isempty(kwargs))
        #plot!(plt,kwargs...)
    else
        plot!(plt,x,pdf(d(θ),x),c=:black,label="",lw=2.0;kwargs...)
    end

    return plt

end
plot_distribution(d::Function,C::Chains,θ=nothing;n=200,kwargs...) = plot_distribution(d,sample(C,n),θ;kwargs...)

function plot_distribution(d::Function,fun::Function,Θ,θ=nothing;n_p=10_000,kwargs...)
    dfun = θ -> kde([fun(p) for p in eachcol(rand(d(θ),n_p))])
    return plot_distribution(dfun,Θ,θ;lims=[0.0,1.0],kwargs...)
end


"""
    plot_distributions(d,Θ::Vector{Vector})
    plot_distributions(d,C::Chains)

Plots the distributions `d(θ).v[1], d(θ).v[2], d(θ).v[3]` with 95% credible intervals.
"""
function plot_distributions(d::Function,Θ::Vector{Vector{Float64}},θ=nothing;
    lims = [quantile(d(mean(Θ)).v[i],[0.0001,0.9999]) for i = 1:3],
    names = ["R","λ","β"],
    kwargs...
)
    return plot(
        [plot_distribution(θ -> d(θ).v[i],Θ,θ;lims=lims[i],xlabel=names[i],kwargs...) for i = 1:3]...,
        layout=grid(1,3)
    )
end
plot_distributions(d::Function,C::Chains,θ=nothing;n=200,kwargs...) = plot_distributions(d,sample(C,n),θ;kwargs...)




##############################################################
## Bivariate parameter plots
##############################################################

function plot_biv_distribution(d::MvDependent;
    dims = [1,2],
    lims = [quantile(d.v[i],[0.0001,0.9999]) for i in dims],
    labs = ["R","λ","β"][dims],
    n_x = 200,
    kwargs...
)
    X = range(lims[1]...,length=n_x)
    Y = range(lims[2]...,length=n_x)
    f = [pdf(marginalize(d;dims),[x,y]) for x in X, y in Y]
    return plot(X,Y,f';xlabel=labs[1],ylabel=labs[2],st=:contourf,lw=0.0,colorbar=false,kwargs...)
end

function plot_biv_distributions(d::MvDependent;
    lims = [quantile(d.v[i],[0.0001,0.9999]) for i in 1:3],
    labs = ["R","λ","β"],
    kwargs...
)
    dims=[[1,2],[1,3],[2,3]]
    return plot([
        plot_biv_distribution(d;
            dims=dims[i],
            lims=lims[dims[i]],
            labs=labs[dims[i]],
            kwargs...
        ) for i = 1:3]..., layout=grid(1,3))
end

##############################################################
## Helpful functions
##############################################################

function vertical_density!(plt,X1,X2,lim,i;fill=false,c=:blue)
    x = range(lim...,length=200)
    y1 = pdf(kde(X1),x)
    y2 = pdf(kde(X2),x)
    sc = maximum(y1)
    plot!(plt,i .+ 0.5 * y1/sc,x,fillrange=i,lc=c,fc=brighten_colour(c,0.7),lw=1.0)
    plot!(plt,i .+ 0.5 * y2/sc,x,c="#333333",lw=1.0)
end

function brighten_colour(hex,α)
    col = parse(Colorant,hex)
    col = convert(RGB,col)
    rgb = [col.r, col.g, col.b]
    rgb = 1 - α .+ α * rgb
    col = RGB(rgb...)
end

function pdf_ci(f::Function,x,Θ::Vector{Vector{Float64}})
    F = hcat([f(θ)(x) for θ in Θ]...)
    l = [quantile(F[i,:],0.025) for i = 1:length(x)]
    u = [quantile(F[i,:],0.975) for i = 1:length(x)]
    return l,u
end