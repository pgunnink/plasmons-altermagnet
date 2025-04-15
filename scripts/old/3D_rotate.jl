using PyPlot
using ProgressMeter
using ForwardDiff
using Unitful
using Revise
using LsqFit
includet("functions_3D.jl")
includet("functions_analysis.jl")
includet("style.jl")
p = ParamsPlasmon()
@info splitting(p)
@info δq([0.01, 0, 0] .* kF(p), p)
##
q = 0.1kF(p)
fq = θ -> q .* [cos(θ), sin(θ), 0]
thetaloc = 0pi / 2 + pi / 8
@info 10u"T" * determine_μ(fq(thetaloc)) * Unitful.μB / (plasmon_pole(fq(thetaloc)) * p.Ef) * 100 |> upreferred
##

Nθ = 200
θrange = range(0, stop=2pi, length=Nθ)
res = zeros(Float64, Nθ)
res_Q = zeros(Float64, Nθ)

for (i, θ) in enumerate(θrange)
    res[i] = plasmon_pole(fq(θ))
    res_Q[i] = Q_factor(fq(θ), res[i])
end
##
fig, ax = plt.subplots(figsize=(fig_width, fig_height))
ax2 = ax.twinx()
@info δq(fq(0), p)
Nθ = 50
for i in 1:4
    # first divide in 4 segments
    θrange = range(2pi / 4 * (i - 1), stop=2pi / 4 * i, length=Nθ) .- 2pi / 8

    res = [plasmon_pole(fq(θ)) for θ in θrange]
    res_Q = [Q_factor(fq(θ)) for θ in θrange]
    if determine_μ(fq(θrange[Nθ/2|>Int])) > 0
        color = "C0"
    else
        color = "C1"
    end
    ax.plot(θrange, res .|> abs, color=color)
    ax2.plot(θrange, res_Q, color=color, "--")
end

# ax.fill_between(θrange, res .- (res_Q ./ res) .^ -1, res .+ (res_Q ./ res) .^ -1, color="C1", alpha=0.5)
ax.set_xlabel(L"\theta")
ax.set_ylabel(L"\omega_p")
ax2.set_ylabel(L"Q")
fig.savefig(joinpath(save_dir, "omega.pdf"))

fig
##
Nθ = 50
fig, ax = plt.subplots(figsize=(fig_width, fig_height), subplot_kw=Dict("projection" => "polar"))

for i in 1:4
    # first divide in 4 segments
    θrange = range(2pi / 4 * (i - 1), stop=2pi / 4 * i, length=Nθ) .- 2pi / 8

    res = [plasmon_pole(fq(θ)) for θ in θrange]

    if determine_μ(fq(θrange[Nθ/2|>Int])) > 0
        color = "C0"
    else
        color = "C1"
    end
    ax.plot(θrange, res .|> abs, color=color)
end
# ax.set_rlabel_position(0)
# ax.set_rticks([0.0, 0.05, 0.1, 0.15])
# ax.set_rticks([0.0, 0.2, 0.4])
# ax.set_ylabel(L"\mu", rotation=0, size=11)
ax.text(pi / 2, 0.02, L"\omega_p")
fig.savefig(joinpath(save_dir, "polar-plot-sign.pdf"))

fig
##
fig, ax = plt.subplots(figsize=(fig_width, fig_height), subplot_kw=Dict("projection" => "polar"))

idx_up = res_μ .> 1
ax.plot(θrange[idx_up], res[idx_up], color="C0", linestyle="-")
idx_down = res_μ .< 1
ax.plot(θrange[idx_down], res[idx_down], color="C1", linestyle="-")

# ax.fill_between(θrange, res .- (res_Q ./ res) .^ -1, res .+ (res_Q ./ res) .^ -1, color="C1", alpha=0.5)
ax.set_rlabel_position(0)
# ax.set_rticks([0.0, 0.05, 0.1, 0.15])
ax.set_rticks([0.0, 0.05, 0.1])
# ax.set_ylabel(L"\mu", rotation=0, size=11)
ax.text(pi / 2, 0.08, L"\omega_p")
fig.savefig(joinpath(save_dir, "polar-plot-omega.pdf"))

fig
##
fig, ax = plt.subplots(figsize=(fig_width, fig_height), subplot_kw=Dict("projection" => "polar"))
ax.plot(θrange, res_Q, color="C0", linestyle="-")
# ax.fill_between(θrange, res .- (res_Q ./ res) .^ -1, res .+ (res_Q ./ res) .^ -1, color="C1", alpha=0.5)
ax.set_rlabel_position(0)
# ax.set_rticks([0.0, 0.05, 0.1, 0.15])
# ax.set_rticks([0.0, 0.05, 0.1])
# ax.set_ylabel(L"\mu", rotation=0, size=11)
ax.text(pi / 2, 0.08, L"Q_p")
fig.savefig(joinpath(save_dir, "polar-plot-Q-factor.pdf"))

fig

##
Nθ = 50
fig, ax = plt.subplots(figsize=(fig_width, fig_height), subplot_kw=Dict("projection" => "polar"))

for i in 1:8
    # first divide in 8 segments
    θrange = range(2pi / 8 * (i - 1), stop=2pi / 8 * i, length=Nθ)
    res = [determine_μ(fq(θ)) for θ in θrange]
    if res[Nθ/2|>Int] > 0
        color = "C0"
    else
        color = "C1"
    end
    ax.plot(θrange, res .|> abs, color=color)
end
ax.set_rlabel_position(0)
# ax.set_rticks([0.0, 0.05, 0.1, 0.15])
# ax.set_rticks([0.0, 0.2, 0.4])
# ax.set_ylabel(L"\mu", rotation=0, size=11)
ax.text(pi / 2, 0.3, L"\mu_p")
fig.savefig(joinpath(save_dir, "polar-plot-magneton.pdf"))

fig
##
