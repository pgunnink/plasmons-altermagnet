using Unitful
using Parameters
import Unitful.ħ
using HCubature
# first do linear
@with_kw mutable struct ParamsMagnons
    K = 0.00u"meV"
    J = 5u"meV"
    δJ = 4u"meV"
    a = 5u"Å"
    Δ = 2sqrt(K * (8J + K))
    Js = 8Δ^-1 * J^2 * a^2
    α = 4a^2 * δJ
    T = 100u"K"
    η = 0.001u"GHz"
    # vs = 
end
fBE(x, p) = 1 / (exp(x / (Unitful.k * p.T)) - 1)
# ωk(k, σ, p) = p.Δ + p.Js * (k[1]^2 + k[2]^2 + k[3]^2) + σ * p.α * (k[1]^2 - k[2]^2)


ωk(k, σ, p) = p.vs * sqrt(k[1]^2 + k[2]^2 + k[3]^2) + σ * p.α * (k[1]^2 - k[2]^2)

# kT(p) = 

A(k, p) = 2p.K + 8p.J
B(k, p) = p.J * (8 - p.a^2 * (k[1]^2 + k[2]^2 + k[3]^2))
ωk(k, σ, p) = sqrt(
    A(k, p)^2 - B(k, p)^2
) + σ * p.α * (k[1]^2 - k[2]^2)

ϕ(k, p) = 0.5atanh(-B(k, p) / A(k, p))

u(k, p) = cosh(ϕ(k, p))
v(k, p) = -sinh(ϕ(k, p))

U(k, p) = [
    u(k, p)^2*v(k, p)^2 (u(k, p)^4+v(k, p)^4)/2
    (u(k, p)^4+v(k, p)^4)/2 u(k, p)^2*v(k, p)^2
]

# function χ0(ω, q, σ, p; rtol=1e-4)
#     k_unit = u"Å^-1"
#     χ_unit = u"meV^-1"
#     integrand = x -> (fBE(ωk(x .* k_unit, σ, p), p) - fBE(ωk(x .* k_unit .+ q, σ, p), p)) / (
#                          ħ * ω + ωk(x .* k_unit, σ, p) - ωk(x .* k_unit .+ q, σ, p) + 1im * ħ * p.η
#                      ) |> χ_unit |> ustrip

#     edge = 0.05pi / p.a |> k_unit |> ustrip
#     a = (-edge, -edge, -edge)
#     b = (edge, edge, edge)

#     (2π)^-2 .* hcubature(integrand, a, b, rtol=rtol) .* p.a^3 .* k_unit^3 .* u"meV^-1 * Å^-3" .|> u"meV^-1 * Å^-3"
# end


function χ0(ω, q, σ, p; N=100)
    k_unit = u"Å^-1"
    χ_unit = u"meV^-1"
    integrand = x -> (fBE(ωk(x .* k_unit, σ, p), p) - fBE(ωk(x .* k_unit .+ q, σ, p), p)) / (
                         ħ * ω + ωk(x .* k_unit, σ, p) - ωk(x .* k_unit .+ q, σ, p) + 1im * ħ * p.η
                     ) |> χ_unit |> ustrip

    edge = 0.05pi / p.a |> k_unit |> ustrip
    a = (-edge, -edge, -edge)
    b = (edge, edge, edge)
    res = 0 + 0im
    krange = range(-edge, stop=edge, length=N)
    dk = krange[2] - krange[1]
    for x in krange
        for y in krange
            for z in krange
                res += integrand([x, y, z])
            end
        end
    end
    (res .* dk^3 .* k_unit^3 .* χ_unit) .|> u"meV^-1 * Å^-3"
end
