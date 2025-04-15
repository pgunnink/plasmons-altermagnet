using Parameters
using Unitful
using LinearAlgebra
import Unitful.ħ, Unitful.me, Unitful.ϵ0
@with_kw mutable struct ParamsPlasmon
    m0 = me
    mstar = 3me
    mDOS = m0 * (mstar^2 / (mstar^2 - m0^2))^(1 / 3)
    Ef = 2u"eV"
    ϵ = 1
end
kF(p) = sqrt(2p.m0 * p.Ef) / ħ
# convention is: 
ϵ(k, σ, p) = ħ^2 * k[1]^2 / 2p.m0 + ħ^2 * k[2]^2 / 2p.m0 + ħ^2 * k[3]^2 / 2p.m0
# such that


N0(p) = p.m0 * kF(p) / (2 * pi^2 * ħ^2)
n0(p) = 2p.Ef * N0(p)
qprime_squared(q, σ, p) = (
    q[1]^2 + q[2]^2 + q[3]^2
)
qprime(q, σ, p) = sqrt(qprime_squared(q, σ, p))


vF(p) = ħ * kF(p) / p.m0



qbar(q, σ, p) = qprime(q, σ, p) / kF(p)
νplus(ω, q, σ, p) = ω / (qprime(q, σ, p) * vF(p)) + qbar(q, σ, p) / 2
νmin(ω, q, σ, p) = ω / (qprime(q, σ, p) * vF(p)) - qbar(q, σ, p) / 2

ωplus(q, σ, p) = ħ * qprime_squared(q, σ, p) / 2p.m0 + qprime(q, σ, p) * vF(p)
ωmin(q, σ, p) = max(0u"Hz", ħ * qprime_squared(q, σ, p) / 2p.m0 - qprime(q, σ, p) * vF(p))


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
vq(q, p) = 4π * Unitful.q^2 / (p.ϵ * ϵ0 * norm(q)^2)

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
