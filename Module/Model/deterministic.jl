#=

    deterministic.jl

    Contains code to solve the deterministic ODE model

           β     λ    pβ
        T --→ S --→ E --→ S + F

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=# 

"""
    solve_ode_model(t,[λ,β,p])

Solve the ODE model using an analytical solution obtained in Mathematica.
Returns [A(t),I(t)].

"""
function solve_ode_model(t,θ::Vector)
    λ,β,p = θ
    A = (exp(-t*β)*(exp(-t*((-1+p)*β+λ))*(-1+p)*p*β^3*λ-(-1+p)*(β-λ)*λ*(p*β+λ)^2+exp(t*β)*((-1+p)*β+λ)*(λ^2*(β+λ)+p^2*β*(-λ^2+β^2*(1+t*λ)+β*λ*(1+t*λ))+p*λ*(-λ^2+β^2*(1+t*λ)+β*λ*(1+t*λ)))))/((β+λ)*((-1+p)*β+λ)*(p*β+λ)^2)
    S = (exp(-t*(p*β+λ))*β*((-1+p)*β*λ+exp(t*(p*β+λ))*p*(β+λ)*((-1+p)*β+λ)-exp(t*((-1+p)*β+λ))*(-1+p)*λ*(p*β+λ)))/((β+λ)*((-1+p)*β+λ)*(p*β+λ))
    I = A - S
    return [A,I]
end
solve_ode_model(t,λ,β,p) = solve_ode_model(t,[λ,β,p])


"""
    solve_ode_model_all_vars(t,[λ,β,p])

Solve the ODE model using matrix exponentiation and return all variables.
Returns [T(t),S(t),E(t),F(t)].

"""
function solve_ode_model_all_vars(t,θ::Vector)
    λ,β,p = θ
    M = [-β 0 0 0; β -λ p*β 0; 0 λ -p*β 0; 0 0 p*β 0]
    x₀ = [λ/(λ+β),β/(λ+β),0,0]
    return exp(M*t) * x₀
end
solve_ode_model_all_vars(t,λ,β,p) = solve_ode_model_all_vars(t,[λ,β,p])


"""
    solve_alt_ode_model(t,[γ,λ,β,p])

Solve the alternate ODE model where antibody binds to receptors at constant rate γ.
"""
function solve_alt_ode_model(t,θ::Vector)
    γ,λ,β,p = θ
    M = [-β 0 0 0 0; β -γ 0 p*β 0; 0 γ -λ 0 0; 0 0 λ -p*β 0; 0 0 0 p*β 0]
    x₀ = [λ/(λ+β),β/(λ+β),0,0,0]
    x = exp(M*t) * x₀
    A = sum(x[3:5])
    I = sum(x[4:5])
    return [A,I]
end
solve_alt_ode_model(t,γ,λ,β,p) = solve_alt_ode_model(t,[γ,λ,β,p])
