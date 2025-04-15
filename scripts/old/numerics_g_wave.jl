using DrWatson
@quickactivate
using Revise
using Roots
using Parameters
using PyPlot
using StaticArrays
using BenchmarkTools
using HCubature
using Unitful
using ProgressMeter
using LinearAlgebra
import Unitful.ħ, Unitful.me, Unitful.ϵ0
includet(scriptsdir("style.jl"))

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
    hcubature(x -> f_iter(x, ω, SVector{3,Float64}(q), σ, Δ, ħm), [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], initdiv=5, norm=x -> sqrt(real(x)^2), rtol=1e-3, maxevals=Int(1e7))[1]
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

χSzSz(ω, q, p) = (S(ω, q, p) - 4vq_natural(q, p) * P(ω, q, p)) / dielectric(ω, q, p)
χnSz(ω, q, p) = D(ω, q, p) / dielectric(ω, q, p)



## first cut along x
q = SVector(0.05, 0.05, 0.0)
N = 100
ωrange = range(0.01, stop=0.15, length=N)
res = @showprogress [dielectric(ω, q, p) for ω in ωrange]



##
fig, ax = plt.subplots()
ax.plot(ωrange, res .|> real, label="Real")
ax.plot(ωrange, res .|> imag, label="Imag")
ax.axhline(0, color="black")
plt.legend()
fig
##
# poles = find_zeros(x -> dielectric(x, q, p) |> real, 0.001, 0.2, rtol=1e-2, xatol=0.005)

##
q = 0.05
Nθ = 80
Nω = 80
ωrange = range(0.001, stop=0.15, length=Nω)

θrange = range(0, stop=2pi, length=Nθ)
res = zeros(typeof(1.0im), Nω, Nθ)

@showprogress for (i, ω) in enumerate(ωrange)
    for (j, θ) in enumerate(θrange)
        res[i, j] = χSzSz(ω, q .* [cos(θ), sin(θ), 0], p)
        # res[i, j] = χnn(ω, fq(q), p) / N0(p) |> upreferred
    end
end


##
toplot = -res .|> imag .|> ustrip

fnorm = plt.matplotlib.colors.Normalize(0, maximum(toplot))
# fnorm = plt.matplotlib.colors.LogNorm(5, maximum(toplot))
fig, ax = plt.subplots(figsize=(fig_width, fig_height), subplot_kw=Dict("projection" => "polar"))


function set_rticks_custom(ax)
    ax.set_rlabel_position(0)

    ax.set_rticks([0.05, 0.1, 0.15], [L"0.05", L"0.1", L"0.15"])
    ax.tick_params(axis="y", labelsize=7)
    ax.tick_params(axis="y", which="minor", bottom=false)
    ax.text(-0.2, 0.08, L"\hbar\omega/\epsilon_F", size=7)
    for label in ax.get_yticklabels()
        label.set_horizontalalignment("center")
    end
end
set_rticks_custom(ax)
# ax.set_xticks([0, 90, 180, 270] / 360 * 2pi)
cb = ax.pcolormesh(θrange, ωrange, toplot, cmap="Blues", norm=fnorm)
cb = fig.colorbar(cb, ax=ax, pad=0.12, shrink=0.8)
cb.ax.set_title(L"-\chi_{S_zS_z}/N_0", fontsize=6)
# ax.plot(θrange, [ωplus(q .* [cos(θ), sin(θ), 0], 1, p) * ħ / p.Ef |> upreferred for θ in θrange], color="white", linestyle="dashed", lw=0.5)
# ax.plot(θrange, [ωplus(q .* [cos(θ), sin(θ), 0], -1, p) * ħ / p.Ef |> upreferred for θ in θrange], color="white", linestyle="dashed", lw=0.5)
ax.set_ylim(extrema(ωrange)...)
ax.grid(false)
NN = 100

# ax.plot(range(-pi / 4, stop=pi / 4, length=NN), ones(NN) * maximum(ωrange) * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)
# ax.plot(range(3pi / 4, stop=5pi / 4, length=NN), ones(NN) * maximum(ωrange) * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)

# ax.plot(range(pi / 4, stop=3pi / 4, length=NN), ones(NN) * maximum(ωrange) * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)
# ax.plot(range(5pi / 4, stop=7pi / 4, length=NN), ones(NN) * maximum(ωrange) * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)

fig.savefig(joinpath(save_dir, "g-wave-polar-plot-chiSzSz.pdf"))

fig