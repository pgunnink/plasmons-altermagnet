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
    Ef = 0.2u"eV"
    t1 = 1u"eV" / Ef |> upreferred
    t2 = 1u"eV" / Ef |> upreferred
    t4 = 1u"eV" / Ef |> upreferred
    μ = 1u"eV" / Ef |> upreferred
    a = 5u"Å"
    ϵ = 1
    vq_natural = Unitful.q^2 / (ϵ * ϵ0) / a / Ef |> upreferred
end

p = ParamsPlasmonNumerics()
##

fx(x, y) = sin(x) + sin(x / 2) * cos(sqrt(3) * y / 2)
fy(x, y) = sqrt(3) * cos(x / 2) * sin(sqrt(3) * y / 2)

ϵ_natural(k, σ, p) = p.t1 * (cos(k[1]) + 2cos(k[1] / 2) * cos(sqrt(3) * k[2] / 2)) + p.t2 * cos(k[3]) - p.μ + σ * p.t4 * sin(k[3]) * fy(k[1], k[2]) * (fy(k[1], k[2])^2 - 3fx(k[1], k[2])^2)
ϵ_split(k, p) = abs(ϵ_natural(k, 1, p) - ϵ_natural(k, -1, p))


##
function f_iter(k, ω, q, σ, p)
    ϵk = ϵ_natural(k, σ, p)
    if ϵk > 1
        return 0im
    else
        ϵkq = ϵ_natural(k .+ q, σ, p)
        return 1 / (ω + ϵk - ϵkq + 1im * 1e-4) + 1 / (-ω + ϵk - ϵkq - 1im * 1e-4)
    end
end

function f_iter_real(k, ω, q, σ, p)
    ϵk = ϵ_natural(k, σ, p)
    if ϵk > 1
        return 0.0
    else
        ϵkq = ϵ_natural(k .+ q, σ, p)
        return real(1 / (ω + ϵk - ϵkq + 1im * 1e-4) + 1 / (-ω + ϵk - ϵkq - 1im * 1e-4))
    end
end
##
function χ0_natural(ω, q, σ, p)
    hcubature(x -> f_iter(x, ω, SVector{3,Float64}(q), σ, p), [-1.5, -1.5, -1.5], [1.5, 1.5, 1.5], initdiv=5, norm=x -> sqrt(real(x)^2), rtol=1e-3, maxevals=Int(1e7))[1]
end




##
function χ0_natural_diff(ω, q, p)
    (2pi)^(-3) .* hcubature(x -> f_iter_real(x, ω, SVector{3,Float64}(q), 1, p) + f_iter_real(x, ω, SVector{3,Float64}(q), -1, p), [-pi, -pi, -pi], [pi, pi, pi], initdiv=5, rtol=1e-3, maxevals=Int(1e7))[1]
end





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

θ = pi / 3
ϕ = 0pi / 3
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
