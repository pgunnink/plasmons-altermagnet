using DrWatson
@quickactivate
using Revise
using LinearAlgebra
using Parameters
using StaticArrays
using HCubature
using Unitful
using ProgressMeter

import Unitful.ħ, Unitful.me, Unitful.ϵ0

@with_kw mutable struct ParamsPlasmon_G_wave
    m0 = 0.4me
    Ef = 0.5u"eV"
    kF = sqrt(2m0 * Ef) / ħ
    Δ = 0.25ħ^2 / m0 / kF^2
    T = 1u"K"
    η = 1u"μeV"
    ϵ = 1
    hbarm0_natural = ħ^2 / 2m0 / Ef * kF^2 |> upreferred
    Δ_natural = Δ / Ef * kF^4 |> upreferred
    vq_natural = Unitful.q^2 / (ϵ * ϵ0) * kF / Ef |> upreferred
end

p = ParamsPlasmon_G_wave()
filename = datadir("bulk_" * savename(p))
##

ϵ(k, σ, p) = (k[1]^2 + k[2]^2 + k[3]^2) * ħ^2 / 2p.m0 + σ * p.Δ * k[1] * k[2] * (k[1]^2 - k[2]^2)


ϵ_natural(k, σ, hbarm0_natural, Δ_natural) = (k[1]^2 + k[2]^2 + k[3]^2) * hbarm0_natural + σ * Δ_natural * k[1] * k[3] * (k[1]^2 - 3k[2]^2)


##
function f_iter(k, ω, q, σ, ħm, Δ)
    ϵk = ϵ_natural(k, σ, ħm, Δ)
    if ϵk > 1
        return 0im
    else
        ϵkq = ϵ_natural(k .+ q, σ, ħm, Δ)
        return 1 / (ω + ϵk - ϵkq + 1im * 1e-4) + 1 / (-ω + ϵk - ϵkq - 1im * 1e-4)
    end
end
##
function χ0_natural(ω, q, σ, Δ, ħm)
    hcubature(x -> f_iter(x, ω, SVector{3,Float64}(q), σ, Δ, ħm), [-1.3, -1.3, -1.3], [1.3, 1.3, 1.3], initdiv=5, norm=x -> sqrt(real(x)^2), rtol=1e-3, maxevals=Int(1e8))[1]
end



χ0_natural(ω, q, σ, p::ParamsPlasmon_G_wave) = χ0_natural(ω, q, σ, p.hbarm0_natural, p.Δ_natural)

# testq = 0.05p.kF .* [1, 0, 0]
# testω = p.Ef / 10 / ħ
# testq_natural = testq ./ p.kF .|> upreferred
# testω_natural = testω / (p.Ef / ħ) |> upreferred

# @btime χ0_natural(testω_natural, testq_natural, 1, p.hbarm0_natural, p.Δ_natural)


##

vq(q, p) = Unitful.q^2 / (p.ϵ * ϵ0 * norm(q)^2)
vq_natural(q, p) = p.vq_natural / norm(q)^2

function dielectric(ω, q, p)
    1 - vq_natural(q, p) * (χ0_natural(ω, q, 1, p) + χ0_natural(ω, q, -1, p))
end
function S(ω, q, p)
    χ0_natural(ω, q, 1, p) + χ0_natural(ω, q, -1, p)
end
function D(ω, q, p)
    χ0_natural(ω, q, 1, p) - χ0_natural(ω, q, -1, p)
end
function P(ω, q, p)
    χ0_natural(ω, q, 1, p) * χ0_natural(ω, q, -1, p)
end


χnn(ω, q, p) = S(ω, q, p) / dielectric(ω, q, p)

function χSzSz(ω, q, p)
    u = χ0_natural(ω, q, 1, p)
    d = χ0_natural(ω, q, -1, p)
    (u + d - 4vq_natural(q, p) * u * d) / (1 - vq_natural(q, p) * (u + d))
end
χnSz(ω, q, p) = D(ω, q, p) / dielectric(ω, q, p)



##
q = 0.05
Nθ = 100
Nω = 100
ωrange = range(0.001, stop=0.15, length=Nω)

θrange = range(0, stop=2pi, length=Nθ)
res = zeros(typeof(1.0im), Nω, Nθ)
q0 = SVector(q / sqrt(2), 0, q / sqrt(2))
@showprogress Threads.@threads for (i, ω) in collect(enumerate(ωrange))
    for (j, θ) in enumerate(θrange)
        qloc = [
            cos(θ) -sin(θ) 0
            sin(θ) cos(θ) 0
            0 0 1
        ] * q0
        res[i, j] = χSzSz(ω, qloc, p)
        # res[i, j] = χnn(ω, fq(q), p) / N0(p) |> upreferred
    end
end
print(filename)
safesave(filename * ".jld2", Dict(
    "θrange" => θrange, "ωrange" => ωrange, "res" => res
))
