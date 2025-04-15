using Parameters
using Unitful
using LinearAlgebra
import Unitful.ħ, Unitful.me, Unitful.ϵ0
# @with_kw mutable struct ParamsPlasmon
#     m0 = 0.4me
#     mstar = 1.2m0
#     mDOS = m0 * mstar / sqrt(mstar^2 - m0^2)
#     ϵ = 1
#     # Ef = sqrt(mstar^2 - m0^2) / mstar * (mstar^2 / (mstar^2 - m0^2))^(1 / 3) * 0.5u"eV"
#     Ef = ((mstar^2 - m0^2) / mstar^2)^(1 / 6) * 0.5u"eV"

# end

include("params.jl")
mDOS(p) = p.m0 * p.mstar / sqrt(p.mstar^2 - p.m0^2)
# convention is: 
ϵ(k, σ, p) = ħ^2 * k[1]^2 * (1 / 2p.m0 + σ / 2p.mstar) + ħ^2 * k[2]^2 * (1 / 2p.m0 - σ / 2p.mstar)
# such that
mx(σ, p) = p.m0 * p.mstar / (p.mstar + σ * p.m0)
my(σ, p) = p.m0 * p.mstar / (p.mstar - σ * p.m0)
N0(p) = mDOS(p) / (2pi * ħ^2)
n0(p) = 2p.Ef * N0(p)
qprime_squared(q, σ, p) = 1 / sqrt(p.mstar^2 - p.m0^2) * (
    q[1]^2 * (p.mstar + σ * p.m0) + q[2]^2 * (p.mstar - σ * p.m0)
)
qprime(q, σ, p) = sqrt(qprime_squared(q, σ, p))


kF(p) = sqrt(2mDOS(p) * p.Ef) / ħ
vF(p) = ħ * kF(p) / mDOS(p)

splitting(p) = abs(ħ^2 * kF(p)^2 / 2mx(1, p) - ħ^2 * kF(p)^2 / 2mx(2, p)) |> u"eV"


critical_angle(p) = 0.5acos(p.mstar / 7p.m0) |> rad2deg
qbar(q, σ, p) = qprime(q, σ, p) / kF(p)
νplus(ω, q, σ, p) = ω / (qprime(q, σ, p) * vF(p)) + qbar(q, σ, p) / 2
νmin(ω, q, σ, p) = ω / (qprime(q, σ, p) * vF(p)) - qbar(q, σ, p) / 2

ωplus(q, σ, p) = ħ * qprime_squared(q, σ, p) / 2mDOS(p) + qprime(q, σ, p) * vF(p)
ωmin(q, σ, p) = max(0u"Hz", ħ * qprime_squared(q, σ, p) / 2mDOS(p) - qprime(q, σ, p) * vF(p))

function χ0(ω, q, σ, p)
    real_part = 1 + 1 / qbar(q, σ, p) * (
        sign(νmin(ω, q, σ, p)) * begin
            νmin(ω, q, σ, p)^2 - 1 > 0 ? sqrt(νmin(ω, q, σ, p)^2 - 1) : 0
        end
        -
        sign(νplus(ω, q, σ, p)) * begin
            νplus(ω, q, σ, p)^2 - 1 > 0 ? sqrt(νplus(ω, q, σ, p)^2 - 1) : 0
        end)

    imag_part = 1 / qbar(q, σ, p) * (
        begin
            1 - νmin(ω, q, σ, p)^2 > 0 ? sqrt(1 - νmin(ω, q, σ, p)^2) : 0
        end
        -
        begin
            1 - νplus(ω, q, σ, p)^2 > 0 ? sqrt(1 - νplus(ω, q, σ, p)^2) : 0
        end
    )
    -N0(p) * (real_part + 1im * imag_part) |> u"meV^-1 * nm^-2"
end




##
function dielectric(ω, q, p)
    1 - 2π * Unitful.q^2 / (p.ϵ * ϵ0 * norm(q)) * (χ0(ω, q, 1, p) + χ0(ω, q, -1, p))
end
function S(ω, q, p)
    χ0(ω, q, 1, p) + χ0(ω, q, -1, p)
end
function D(ω, q, p)
    χ0(ω, q, 1, p) - χ0(ω, q, -1, p)
end
function P(ω, q, p)
    χ0(ω, q, 1, p) * χ0(ω, q, -1, p)
end
vq(q, p) = 2π * Unitful.q^2 / (p.ϵ * ϵ0 * norm(q))


χnn(ω, q, p) = S(ω, q, p) / dielectric(ω, q, p)

χSzSz(ω, q, p) = (S(ω, q, p) - 4vq(q, p) * P(ω, q, p)) / dielectric(ω, q, p)
χnSz(ω, q, p) = D(ω, q, p) / dielectric(ω, q, p)


##

Vq(q, p) = vq(q, p) * N0(p) / (1 + 2vq(q, p) * N0(p)) |> upreferred
# ωs(q, p) = min(qprime(q, -1, p) * sqrt(4 / 3) * vF(p), qprime(q, 1, p) * sqrt(4 / 3) * vF(p))

ωs(q, p) = min(
    qprime(q, -1, p) * vF(p) * sqrt(
        1 / (1 - Vq(q, p)^2) +
        qprime(q, -1, p)^2 / (4kF(p)^2 * Vq(q, p)^2)
    ),
    qprime(q, 1, p) * vF(p) * sqrt(
        1 / (1 - Vq(q, p)^2) +
        qprime(q, 1, p)^2 / (4kF(p)^2 * Vq(q, p)^2)
    ))

function ωs_corrected(q, p)
    if (ωplus(q, -1, p) < ωs(q, p) < ωplus(q, +1, p)) || (ωplus(q, +1, p) < ωs(q, p) < ωplus(q, -1, p))
        ωs(q, p)
    else
        NaN
    end
end

function ωs_plus_phase(q, p)
    if (ωplus(q, -1, p) < ωs(q, p) < ωplus(q, +1, p))
        ωs(q, p), 1
    elseif (ωplus(q, +1, p) < ωs(q, p) < ωplus(q, -1, p))
        ωs(q, p), -1
    else
        ωs(q, p), NaN
    end
end

ωp(q, p) = sqrt(
    π * Unitful.q^2 * n0(p) / (p.ϵ * ϵ0 * mDOS(p))
)

function phase(q, p)
    -vq(q, p) * N0(p) / (1 + vq(q, p) * N0(p)) |> upreferred
end

η(θ, σ, p) = sqrt(mDOS(p)) * sqrt(
    p.m0^-1 + σ * p.mstar^-1 * cos(2θ))
# function vs(θ, p)
#     if abs(cos(θ)) > abs(sin(θ))
#         if η(θ, -1, p) / η(θ, 1, p) < √(3) / 2
#             return 2 / √(3) * abs(η(θ, -1, p)) * vF(p)
#         else
#             return NaN
#         end
#     elseif abs(cos(θ)) < abs(sin(θ))
#         if η(θ, 1, p) / η(θ, -1, p) < √(3) / 2
#             return 2 / √(3) * abs(η(θ, 1, p)) * vF(p)
#         else
#             return NaN
#         end
#     end
#     return return NaN
# end

function vs(θ, p)
    if abs(cos(θ)) > abs(sin(θ))
        if η(θ, -1, p) / η(θ, 1, p) < √(3) / 2
            return 2 / √(3) * η(θ, -1, p) * vF(p)
        end
    elseif abs(cos(θ)) < abs(sin(θ))
        if η(θ, 1, p) / η(θ, -1, p) < √(3) / 2
            return 2 / √(3) * abs(η(θ, 1, p)) * vF(p)
        end
    end
    return NaN
end

function vs(q, θ, p)
    if abs(cos(θ)) > abs(sin(θ))
        if η(θ, -1, p) / η(θ, 1, p) < √(3) / 2
            return η(θ, -1, p) * vF(p) * sqrt(4 / 3 + (η(θ, -1, p) * q / kF(p))^2)
        end
    elseif abs(cos(θ)) < abs(sin(θ))
        if η(θ, 1, p) / η(θ, -1, p) < √(3) / 2
            return 2 / √(3) * abs(η(θ, 1, p)) * vF(p)
        end
    end
    return NaN
end
# vs(θ, p) = 2 / √(3) * η(θ, -1, p) * vF(p)
Q(θ, p) = -4η(θ, -1, p)^2 + 3η(θ, 1, p)^2 > 0 ? 3√(-4η(θ, -1, p)^2 + 3η(θ, 1, p)^2) / η(θ, -1, p) : NaN

## 
μ(q, θ, p) = ħ * vs(θ, p) * q * 2 * Unitful.μB / p.Ef * sign(abs(cos(θ)) - abs(sin(θ)))