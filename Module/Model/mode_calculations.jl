function find_mode(C::Chains;n=1000,k=10)
    Θ = sample(C,n)
    dk = [sort([norm(Θ[j] - Θ[i]) for j = 1:n])[k+1] for i = 1:n]
    return Θ[findmin(dk)[2]]
end