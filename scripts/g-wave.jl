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
include(scriptsdir("style.jl"))
import Unitful.ħ, Unitful.me, Unitful.ϵ0

# @with_kw mutable struct ParamsPlasmon_G_wave
#     m0 = 0.4me
#     Ef = 0.5u"eV"
#     kF = sqrt(2m0 * Ef) / ħ
#     Δ = 0.35ħ^2 / m0 / kF^2
#     ϵ = 1
#     hbarm0_natural = ħ^2 / 2m0 / Ef * kF^2 |> upreferred
#     Δ_natural = Δ / Ef * kF^4 |> upreferred
#     vq_natural = Unitful.q^2 / (ϵ * ϵ0) * kF / Ef |> upreferred
# end
# vq_natural = Unitful.q^2 / (ϵ * ϵ0) * kF / Ef |> upreferred

##



ϵ_natural(k, σ, Js, Δ) = (k[1]^2 + k[2]^2 + k[3]^2) * Js + σ * Δ * k[1] * k[3] * (k[1]^2 - 3k[2]^2)
kmax(Js, Δ, Ef) = 2sqrt(2Js - sqrt(4Js^2 - 3sqrt(3) * Ef * Δ)) / sqrt(Δ) / 3^(3 / 4)



##
function f_iter(k, ω, q, σ, Js, Δ, Ef, kmax2)
    ϵk = ϵ_natural(k, σ, Js, Δ)
    if ϵk > Ef || (k[1]^2 + k[2]^2 + k[3]^2) > kmax2
        return 0im
    else
        ϵkq = ϵ_natural(k .+ q, σ, Js, Δ)
        return 1 / (ω + ϵk - ϵkq + 1im * 1e-3) + 1 / (-ω + ϵk - ϵkq - 1im * 1e-3)
    end
end


##
function χ0_natural(ω, q, σ, Js, Δ, Ef)
    klim = kmax(Js, Δ, Ef)
    (2pi)^(-3) .* hcubature(x -> f_iter(x, ω, q, σ, Js, Δ, Ef, klim^2), [-klim, -klim, -klim], [klim, klim, klim], initdiv=5, norm=x -> sqrt(real(x)^2), rtol=1e-3, maxevals=Int(1e8))
end
##

function epsilon_natural(ω, q, Js, Δ, vq, Ef)
    klim = kmax(Js, Δ, Ef)

    (2pi)^(-3) .* hcubature(x -> 1.0 - vq / norm(q)^2 * (f_iter(x, ω, q, 1, Js, Δ, Ef, klim^2) + f_iter(x, ω, q, -1, Js, Δ, Ef, klim^2)), [-klim, -klim, -klim], [klim, klim, klim], rtol=1e-3, maxevals=Int(1e8), norm=x -> sqrt(real(x)^2))
end



##
Ef = 1
Js = 1
Δ = 4 * Js^2 / (Ef * 3sqrt(3))

vq = 80

q = 0.05
θ = pi / 3
ϕ = 0pi / 3
qloc = SVector(q * sin(θ) * cos(ϕ), q * sin(θ) * sin(ϕ), q * cos(θ))
# qloc = SVector(q / sqrt(3), q / sqrt(3), q / sqrt(3))

fig, ax = plt.subplots()
qrange = -1.2kmax(Js, Δ, Ef):0.01:1.2kmax(Js, Δ, Ef)
vector_spher(θ, ϕ) = [sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)]
ax.plot(qrange, [ϵ_natural(vector_spher(θ, ϕ) * q, 1, Js, Δ) for q in qrange])
ax.plot(qrange, [ϵ_natural(vector_spher(θ, ϕ) * q, -1, Js, Δ) for q in qrange])
ax.axhline(Ef)
fig

##

N = 200
ωrange = range(0.001, stop=0.25, length=N)

resup = zeros(ComplexF64, N)
resdown = zeros(ComplexF64, N)

errup = zeros(Float64, N)
errdown = zeros(Float64, N)

@showprogress Threads.@threads for i in 1:N
    resup[i], errup[i] = χ0_natural(ωrange[i], qloc, 1, Js, Δ, Ef)
    resdown[i], errdown[i] = χ0_natural(ωrange[i], qloc, -1, Js, Δ, Ef)
end
##
fig, axs = plt.subplots(2, 1, sharex=true)
ax = axs[1]
ax.axhline(0, color="black")

ax.plot(ωrange, .-resup .|> real, color="C0", label=L"\mathrm{Re}[\chi_\uparrow]")
ax.plot(ωrange, .-resup .|> imag, color="C0", linestyle="dashed", label=L"\mathrm{Im}[\chi_\uparrow]")
ax.plot(ωrange, .-resdown .|> real, color="C1", label=L"\mathrm{Re}[\chi_\downarrow]")
ax.plot(ωrange, .-resdown .|> imag, color="C1", linestyle="dashed", label=L"\mathrm{Im}[\chi_\downarrow]")
ax.set_ylabel(L"-\chi_{\sigma}\,[\mathrm{arb. units}]")
ax.legend(loc="center right", frameon=true, fancybox=false, edgecolor="white")

ax = axs[2]
ax.axhline(0, color="black")
die = 1.0 .- vq / norm(q)^2 .* (resup + resdown)
ax.errorbar(ωrange, die .|> real, color="C0")
ax.set_xlabel(L"\omega")
ax.set_ylabel(L"\mathrm{Re}[\epsilon]")


ax.axhline(0, color="black")
fig.subplots_adjust(left=0.2, right=0.95, bottom=0.2)
fig.savefig(joinpath(save_dir, "g-wave.pdf"))

fig


# ##
# fig, axs = plt.subplots(3, 1)
# ax = axs[1]
# ax.errorbar(ωrange, .-resup .|> real, color="C0", yerr=errup)
# ax.plot(ωrange, .-resup .|> imag, color="C0", linestyle="dashed")
# ax.errorbar(ωrange, .-resdown .|> real, color="C1", yerr=errdown)
# ax.plot(ωrange, .-resdown .|> imag, color="C1", linestyle="dashed")

# ax.set_xlabel(L"\omega")
# ax.set_ylabel(L"\epsilon")

# ax = axs[2]
# ax.axhline(0, color="black")
# die = 1.0 .- vq / norm(q)^2 .* (resup + resdown)
# ax.errorbar(ωrange, die .|> real, color="C0", yerr=errup)

# ax = axs[3]
# chizz = (resup .+ resdown - 4vq .* resup .* resdown) ./ die
# ax.errorbar(ωrange, die .|> imag, color="C0", yerr=errup)


# ax.axhline(0)

# fig
