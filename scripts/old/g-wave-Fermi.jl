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
    Δ = 0.4ħ^2 / m0 / kF^2
    T = 1u"K"
    η = 1u"μeV"
    ϵ = 1
    hbarm0_natural = ħ^2 / 2m0 / Ef * kF^2 |> upreferred
    Δ_natural = Δ / Ef * kF^4 |> upreferred
    vq_natural = Unitful.q^2 / (ϵ * ϵ0) * kF / Ef |> upreferred
end

p = ParamsPlasmon_G_wave()
##

ϵ(k, σ, p) = (k[1]^2 + k[2]^2 + k[3]^2) * ħ^2 / 2p.m0 + σ * p.Δ * k[1] * k[2] * (k[1]^2 - k[2]^2)


ϵ_natural(k, σ, hbarm0_natural, Δ_natural) = (k[1]^2 + k[2]^2 + k[3]^2) * hbarm0_natural + σ * Δ_natural * k[1] * k[2] * (k[1]^2 - k[2]^2)

ϵ_natural(k, σ, hbarm0_natural, Δ_natural) = (k[1]^2 + k[2]^2 + k[3]^2) * hbarm0_natural + σ * Δ_natural * k[1] * k[3] * (k[1]^2 - 3k[2]^2)

##
fig, ax = plt.subplots()
krange = range(0, stop=1.5, length=100)
θ = pi / 6
split = 0.7
kz = 0.5
ax.plot(krange, [ϵ_natural([cos(θ) * k, 0, sin(θ) * k], 1, 1, split) for k in krange])
ax.plot(krange, [ϵ_natural([cos(θ) * k, 0, sin(θ) * k], -1, 1, split) for k in krange])

ax.plot(-krange, [ϵ_natural(-[cos(θ) * k, 0, sin(θ) * k], 1, 1, split) for k in krange])
ax.plot(-krange, [ϵ_natural(-[cos(θ) * k, 0, sin(θ) * k], -1, 1, split) for k in krange])
ax.set_ylim(0, 2.5)
ax.axhline(1)
fig
##
fig, ax = plt.subplots()
krange = range(0, stop=1.2, length=100)
θ = 0pi
split = 0.5
kz = 0.5
ax.plot(krange, [ϵ_natural([cos(θ) * k, sin(θ) * k, kz], 1, 1, split) for k in krange])
ax.plot(krange, [ϵ_natural([cos(θ) * k, sin(θ) * k, kz], -1, 1, split) for k in krange])

ax.plot(-krange, [ϵ_natural([-cos(θ) * k, sin(θ) * k, kz], 1, 1, split) for k in krange], "C0")
ax.plot(-krange, [ϵ_natural([-cos(θ) * k, sin(θ) * k, kz], -1, 1, split) for k in krange])
ax.set_ylim(0, 2.5)
ax.axhline(1)
fig

##
Nx = 50
Ny = 60

resup = zeros(Float64, Nx, Ny)
resdown = zeros(Float64, Nx, Ny)
krangex = range(-1.2, stop=1.2, length=Nx)
krangey = range(-1.2, stop=1.2, length=Ny)
for (i, kx) in enumerate(krangex)
    for (j, ky) in enumerate(krangey)
        resup[i, j] = ϵ_natural([kx, ky, 0ky], 1, p.hbarm0_natural, p.Δ_natural)
        resdown[i, j] = ϵ_natural([kx, ky, 0ky], -1, p.hbarm0_natural, p.Δ_natural)

    end
end

resup[resup.>1] .= NaN
resdown[resdown.>1] .= NaN
##
fig, ax = plt.subplots()
ax.contour(krangex, krangey, resup')
# ax.contour(krangex, krangey, resdown')

fig