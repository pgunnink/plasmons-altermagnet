using DrWatson
@quickactivate
using Revise
using PyPlot
using LinearAlgebra
using Parameters
using StaticArrays
using HCubature
using Unitful
using ProgressMeter

import Unitful.ħ, Unitful.me, Unitful.ϵ0

@with_kw mutable struct ParamsPlasmonNumerics
    m0 = 0.4me
    mstar = m0 * 1.1
    Ef = 0.2u"eV"
    a = 5u"Å"
    mDOS = m0 * (mstar^2 / (mstar^2 - m0^2))^(1 / 3)
    ϵ = 1
    kiso = ħ^2 / 2m0 / Ef / a^2 |> upreferred
    kx = ħ^2 / 2mstar / Ef / a^2 |> upreferred
    vq_natural = Unitful.q^2 / (ϵ * ϵ0) / a / Ef |> upreferred
end

p = ParamsPlasmonNumerics()
##

ϵ_natural(k, σ, kiso, kx) = (k[1]^2 + k[2]^2 + k[3]^2) * kiso + σ * kx * k[1]^2 - σ * kx * k[2]^2


ϵ_natural(k, σ, kiso, kx) = (3 - cos(k[1]) - cos(k[2]) - cos(k[3])) * kiso + σ * kx * (1 - cos(k[1])) - σ * kx * (1 - cos(k[2]))


ϵ_natural(k, σ, p) = ϵ_natural(k, σ, p.kiso, p.kx)
ϵ_split(k, p) = abs(ϵ_natural(k, 1, p) - ϵ_natural(k, -1, p))


##
function f_iter(k, ω, q, σ, kiso, kx)
    ϵk = ϵ_natural(k, σ, kiso, kx)
    if ϵk > 1
        return 0im
    else
        ϵkq = ϵ_natural(k .+ q, σ, kiso, kx)
        return 1 / (ω + ϵk - ϵkq + 1im * 1e-4) + 1 / (-ω + ϵk - ϵkq - 1im * 1e-4)
    end
end

function f_iter_real(k, ω, q, σ, kiso, kx)
    ϵk = ϵ_natural(k, σ, kiso, kx)
    if ϵk > 1
        return 0.0
    else
        ϵkq = ϵ_natural(k .+ q, σ, kiso, kx)
        return real(1 / (ω + ϵk - ϵkq + 1im * 1e-4) + 1 / (-ω + ϵk - ϵkq - 1im * 1e-4))
    end
end
##
function χ0_natural(ω, q, σ, kiso, kx)
    hcubature(x -> f_iter(x, ω, SVector{3,Float64}(q), σ, kiso, kx), [-1.5, -1.5, -1.5], [1.5, 1.5, 1.5], initdiv=5, norm=x -> sqrt(real(x)^2), rtol=1e-3, maxevals=Int(1e7))[1]
end



χ0_natural(ω, q, σ, p::ParamsPlasmonNumerics) = χ0_natural(ω, q, σ, p.kiso, p.kx)

##
function χ0_natural_diff(ω, q, kiso, kx)
    (2pi)^(-3) .* hcubature(x -> f_iter_real(x, ω, SVector{3,Float64}(q), 1, kiso, kx) + f_iter_real(x, ω, SVector{3,Float64}(q), -1, kiso, kx), [-pi, -pi, -pi], [pi, pi, pi], initdiv=5, rtol=1e-3, maxevals=Int(1e8))[1]
end



χ0_natural_diff(ω, q, p::ParamsPlasmonNumerics) = χ0_natural_diff(ω, q, p.kiso, p.kx)


##
θ = 0
q = 0.01
q0 = SVector(q / sqrt(2), 0, q / sqrt(2))
q0 = SVector(q / sqrt(3), q / sqrt(3), q / sqrt(3))

qloc = [
    cos(θ) -sin(θ) 0
    sin(θ) cos(θ) 0
    0 0 1
] * q0

θ = pi / 2
ϕ = pi / 2
qloc = SVector(q * sin(θ) * cos(ϕ), q * sin(θ) * sin(ϕ), q * cos(θ))

# qloc = SVector(q / sqrt(3), q / sqrt(3), q / sqrt(3))

@info ϵ_split(qloc ./ q, p)

fig, ax = plt.subplots()
qrange = -2pi:0.01:2pi
ax.plot(qrange, [ϵ_natural(SVector(q * sin(θ) * cos(ϕ), q * sin(θ) * sin(ϕ), q * cos(θ)), 1, p) for q in qrange])
ax.plot(qrange, [ϵ_natural(SVector(q * sin(θ) * cos(ϕ), q * sin(θ) * sin(ϕ), q * cos(θ)), -1, p) for q in qrange])
@info ϵ_natural(SVector(1.2 * sin(θ) * cos(ϕ), 1.2 * sin(θ) * sin(ϕ), 1.2 * cos(θ)), -1, p), ϵ_natural(SVector(1.2 * sin(θ) * cos(ϕ), 1.2 * sin(θ) * sin(ϕ), 1.2 * cos(θ)), 1, p)
ax.axhline(1)
fig
##
N = 150
ωrange = range(0.02, stop=0.0283, length=N)


resdiff = @showprogress [χ0_natural_diff(ω, qloc, p) for ω in ωrange]
##

vq_natural(q, p) = p.vq_natural / norm(q)^2

dres = (1.0 .- vq_natural(q, p) .* resdiff)
fig, ax = plt.subplots(1, 1)
ax.plot(ωrange, dres, color="C0")
ax.set_xlabel(L"\omega")
ax.set_ylabel(L"\epsilon")
ax.axhline(0)

fig



##
resup = @showprogress [χ0_natural(ω, qloc, 1, p) for ω in ωrange]
resdown = @showprogress [χ0_natural(ω, qloc, -1, p) for ω in ωrange]

##
fig, ax = plt.subplots()

ax.plot(ωrange, .-resup .|> imag, color="C0")
ax.plot(ωrange, .-resdown .|> imag, color="C1")
ax.plot(ωrange, .-resup .|> real, "--", color="C0")
ax.plot(ωrange, .-resdown .|> real, "--", color="C1")

fig

##


vq_natural(q, p) = p.vq_natural / norm(q)^2

χres = (resup .+ resdown .- 4vq_natural(q, p) .* resup .* resdown) ./ (1.0 .- vq_natural(q, p) .* (resup .+ resdown))
dres = (1.0 .- vq_natural(q, p) .* (resup .+ resdown))
fig, axs = plt.subplots(2, 1)
axs[1].plot(ωrange, dres .|> real, color="C0")
axs[1].plot(ωrange, dres .|> imag, "--", color="C0")

axs[1].axhline(0)

axs[2].plot(ωrange, χres .|> imag, color="C0")
# ax.plot(ωrange, χres .|> real, "--", color="C0")

fig

##

χ0_natural(ω, q, σ, p::ParamsPlasmon_G_wave) = χ0_natural(ω, q, σ, p.hbarm0_natural, p.Δ_natural)

# testq = 0.05p.kF .* [1, 0, 0]
# testω = p.Ef / 10 / ħ
# testq_natural = testq ./ p.kF .|> upreferred
# testω_natural = testω / (p.Ef / ħ) |> upreferred

# @btime χ0_natural(testω_natural, testq_natural, 1, p.hbarm0_natural, p.Δ_natural)


##

vq(q, p) = Unitful.q^2 / (p.ϵ * ϵ0 * norm(q)^2)

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

χSzSz(ω, q, p) = (S(ω, q, p) - 4vq_natural(q, p) * P(ω, q, p)) / dielectric(ω, q, p)
χnSz(ω, q, p) = D(ω, q, p) / dielectric(ω, q, p)

##
N = 100
ωrange = range(0.01, stop=0.1, length=N)
q = 0.05
q0 = SVector(q / sqrt(2), 0, q / sqrt(2))
θ = pi / 2
qloc = [
    cos(θ) -sin(θ) 0
    sin(θ) cos(θ) 0
    0 0 1
] * q0

res = @showprogress [dielectric(ω, qloc, p) for ω in ωrange]

##
fig, ax = plt.subplots()
ax.plot(ωrange, res .|> imag, color="C0")
ax.plot(ωrange, res .|> real, color="C1")
fig
