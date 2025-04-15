using Unitful
using LinearAlgebra
import Unitful.ħ, Unitful.me, Unitful.ϵ0
using Roots
include("params.jl")
mDOS(p) = p.m0 * (p.mstar^2 / (p.mstar^2 - p.m0^2))^(1 / 3)
kF(p) = sqrt(2mDOS(p) * p.Ef) / ħ
# convention is: 
ϵ(k, σ, p) = ħ^2 * k[1]^2 * (1 / 2p.m0 + σ / 2p.mstar) + ħ^2 * k[2]^2 * (1 / 2p.m0 - σ / 2p.mstar) + ħ^2 * k[3]^2 / 2p.m0

# such that
mx(σ, p) = p.m0 * p.mstar / (p.mstar + σ * p.m0)
my(σ, p) = p.m0 * p.mstar / (p.mstar - σ * p.m0)
mz(p) = p.m0


eccentricity(σ, p) = sqrt(1 - my(σ, p) / mx(σ, p))

P(p) = (sqrt(mx(1, p)) + sqrt(mx(-1, p))) / sqrt(mDOS(p))
Q(p) = 4 / pi / P(p) * sqrt(mx(-1, p) / mDOS(p))
N0(p) = mDOS(p) * kF(p) / (2 * pi^2 * ħ^2)


n0(p) = 2p.Ef * N0(p) / 3
N0D(p) = -4 / 3 * n0(p) / p.Ef^2 |> u"meV^-2 * nm^-3"
Δ(x, p) = 2 * Unitful.μB * x * N0D(p) / N0(p)
dΔdB(p) = 2 * Unitful.μB * N0D(p) / N0(p)
qprime_squared(q, σ, p) = (
    q[1]^2 * (mDOS(p) / mx(σ, p)) + q[2]^2 * (mDOS(p) / my(σ, p)) + q[3]^2 * (mDOS(p) / mz(p))
)
qprime(q, σ, p) = sqrt(qprime_squared(q, σ, p))

δq(q, p) = qprime(q, -1, p) / qprime(q, 1, p)

vF(p) = ħ * kF(p) / mDOS(p)

splitting(p) = abs(ħ^2 * kF(p)^2 / 2mx(1, p) - ħ^2 * kF(p)^2 / 2mx(2, p)) |> u"eV"

plasmon_gap(p) = sqrt(2n0(p) * Unitful.q^2 / mDOS(p) / Unitful.ϵ0) |> u"THz"

plasmon_gap_in_Ef(p) = ħ * plasmon_gap(p) / p.Ef |> upreferred


qbar(q, σ, p) = qprime(q, σ, p) / kF(p)

ωs(q, p) = qprime(q, -1, p) * vF(p) * 1.04438

νplus(ω, q, σ, p) = ω / (qprime(q, σ, p) * vF(p)) + qbar(q, σ, p) / 2
νmin(ω, q, σ, p) = ω / (qprime(q, σ, p) * vF(p)) - qbar(q, σ, p) / 2

ωplus(q, σ, p) = ħ * qprime_squared(q, σ, p) / 2mDOS(p) + qprime(q, σ, p) * vF(p)
ωminabs(q, σ, p) = abs(ħ * qprime_squared(q, σ, p) / 2mDOS(p) - qprime(q, σ, p) * vF(p))
ωmin(q, σ, p) = max(0u"Hz", ħ * qprime_squared(q, σ, p) / 2mDOS(p) - qprime(q, σ, p) * vF(p))

function find_pole(q, p)
    f = x -> dielectric(x * p.Ef / ħ, q, p) |> real
    ω0 = min(ωplus(q, -1, p), ωplus(q, 1, p)) * 0.8 / (p.Ef / ħ)
    find_zero(f, ω0) * (p.Ef / ħ)
end

flogabs(x) = log(abs((x + 1) / (x - 1)))
function χ0(ω, q, σ, p)
    real_part = 1 / 2 - (1 - νmin(ω, q, σ, p)^2) / 4qbar(q, σ, p) * flogabs(νmin(ω, q, σ, p)) + (1 - νplus(ω, q, σ, p)^2) / 4qbar(q, σ, p) * flogabs(νplus(ω, q, σ, p))

    imag_part = π / 4qbar(q, σ, p) * (
        begin
            1 - νmin(ω, q, σ, p)^2 > 0 ? (1 - νmin(ω, q, σ, p)^2) : 0
        end
        -
        begin
            1 - νplus(ω, q, σ, p)^2 > 0 ? (1 - νplus(ω, q, σ, p)^2) : 0
        end
    )
    -N0(p) * (real_part + 1im * imag_part) |> u"meV^-1 * nm^-3"
end

##
vq(q, p) = Unitful.q^2 / (p.ϵ * ϵ0 * norm(q)^2)

function dielectric(ω, q, p)
    1 - vq(q, p) * (χ0(ω, q, 1, p) + χ0(ω, q, -1, p))
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


χnn(ω, q, p) = S(ω, q, p) / dielectric(ω, q, p)

χSzSz(ω, q, p) = (S(ω, q, p) - 4vq(q, p) * P(ω, q, p)) / dielectric(ω, q, p)
χnSz(ω, q, p) = D(ω, q, p) / dielectric(ω, q, p)

χupup(ω, q, p) = χnSz(ω, q, p) / 2 + (χnn(ω, q, p) + χSzSz(ω, q, p)) / 4
χdowndown(ω, q, p) = -χnSz(ω, q, p) / 2 + (χnn(ω, q, p) + χSzSz(ω, q, p)) / 4
χupdown(ω, q, p) = (χnn(ω, q, p) - χSzSz(ω, q, p)) / 4

##

