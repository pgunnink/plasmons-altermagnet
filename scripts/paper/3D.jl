
using DrWatson
using PyPlot
using ProgressMeter
using ForwardDiff
using Unitful
using Revise
using LsqFit
includet(scriptsdir("functions_3D.jl"))
includet(scriptsdir("functions_analysis.jl"))

includet(scriptsdir("style.jl"))


p = ParamsPlasmon()
q = 0.05kF(p)

##
N = 100
krange = range(0kF(p), stop=3kF(p), length=N)
fig, ax = plt.subplots()
ax.plot(krange ./ kF(p), [ϵ([k, 0k, 0k], 1, p) |> u"eV" |> ustrip for k in krange], color="red")
ax.plot(krange ./ kF(p) .|> upreferred, [ϵ([k, 0k, 0k], -1, p) |> u"eV" |> ustrip for k in krange], color="blue")


ax.plot(.-krange ./ kF(p), [ϵ([0k, k, 0k], 1, p) |> u"eV" |> ustrip for k in .-krange], color="red", label=L"\sigma=\uparrow")
ax.plot(.-krange ./ kF(p) .|> upreferred, [ϵ([0k, k, 0k], -1, p) |> u"eV" |> ustrip for k in .-krange], color="blue", label=L"\sigma=\downarrow")
ax.set_xlabel(L"k/k_F")
ax.set_ylabel(L"\epsilon_{\mathbf k}^\sigma\,(\mathrm{eV})")
ax.axhline(p.Ef |> u"eV" |> ustrip, color="black", linestyle="dashed")
ax.set_ylim(-0.1, 0.3ϵ([0krange[end], krange[end], 0krange[end]], -1, p) |> u"eV" |> ustrip)
ax.text(krange[end] / kF(p) |> ustrip, -0.7, L"\mathbf k \parallel \hat x", ha="center")
ax.text(-krange[end] / kF(p) |> ustrip, -0.7, L"\mathbf k \parallel \hat y", ha="center")
ax.legend(loc="upper right")
ax.margins(x=0)
fig.subplots_adjust(bottom=0.2)
fig.savefig(joinpath(save_dir, "dispersion.pdf"))

fig
##
@info splitting(p)
Nq = 3 * 100
Nω = 3 * 80
qrange = range(0.001kF(p), stop=0.1kF(p), length=Nq)
ωrange = range(0.001p.Ef / ħ, stop=0.3p.Ef / ħ, length=Nω)


fq = q -> [q, 0q, 0q]

res = zeros(typeof(1.0im), Nω, Nq)
for (i, ω) in enumerate(ωrange)
    for (j, q) in enumerate(qrange)
        res[i, j] = χSzSz(ω, fq(q), p) / N0(p) |> upreferred
        # res[i, j] = χnn(ω, fq(q), p) / N0(p) |> upreferred
    end
end
##
fig, ax = plt.subplots(figsize=(fig_width, 0.8fig_height))
toplot = -res .|> imag .|> ustrip
toplot[toplot.<1e-5] .= 1e-5
fnorm = plt.matplotlib.colors.LogNorm(1e-1, maximum(toplot))
fnorm = plt.matplotlib.colors.Normalize(0, maximum(toplot))

cb = ax.pcolormesh(qrange ./ kF(p) .|> upreferred, ωrange .* ħ ./ p.Ef .|> upreferred, toplot, norm=fnorm, cmap="Blues", shading="nearest")
ax.plot(qrange ./ kF(p) .|> upreferred, [ωplus(fq(q), -1, p) * ħ / p.Ef |> upreferred for q in qrange], color="black", linestyle="dashed", lw=1, alpha=0.5)
ax.plot(qrange ./ kF(p) .|> upreferred, [ωplus(fq(q), 1, p) * ħ / p.Ef |> upreferred for q in qrange], color="black", linestyle="dashed", lw=1, alpha=0.5)

ax.text(0.052, 0.26, L"\omega_{+\uparrow}")
ax.text(0.089, 0.05, L"\omega_{+\downarrow}")
# ax.plot(qrange ./ kF(p) .|> upreferred, [ωs(fq(q), p) * ħ / p.Ef |> upreferred for q in qrange], color="black", linestyle="--", lw=0.5)
cb = fig.colorbar(cb, ax=ax)
cb.ax.set_title(L"-\mathrm{Im}[\chi_{S_zS_z}(\mathbf q,\omega)]/N_0")
ax.set_xlabel(L"q/k_F")
ax.set_ylabel(L"\hbar\omega/\epsilon_F")
ax.set_ylim(extrema(ωrange .* ħ ./ p.Ef .|> upreferred)...)
# ax.plot(qrange ./ kF(p) .|> upreferred, vF(p) .* qrange .* ħ ./ p.Ef, color="black", linestyle="dashed", lw=0.5)
ax.axvline(q / kF(p), color="black", alpha=0.3, linestyle="dotted")
fig.subplots_adjust(left=0.15, bottom=0.25, top=0.8, right=0.98)

fig.savefig(joinpath(save_dir, "SzSz_3D.pdf"))
fig
##

Nθ = 300
Nω = 300
ωrange = range(0.001p.Ef / ħ, stop=0.15p.Ef / ħ, length=Nω)

θrange = range(0, stop=2pi, length=Nθ)
res = zeros(typeof(1.0im), Nω, Nθ)
for (i, ω) in enumerate(ωrange)
    for (j, θ) in enumerate(θrange)
        res[i, j] = χSzSz(ω, q .* [cos(θ), sin(θ), 0], p) / N0(p) |> upreferred
        # res[i, j] = χnn(ω, fq(q), p) / N0(p) |> upreferred
    end
end
##
fnorm = plt.matplotlib.colors.Normalize(0, maximum(.-res .|> imag))

fig, ax = plt.subplots(figsize=(fig_width, fig_height), subplot_kw=Dict("projection" => "polar"))


set_rticks_custom(ax)
# ax.set_xticks([0, 90, 180, 270] / 360 * 2pi)
cb = ax.pcolormesh(θrange, ωrange .* ħ ./ p.Ef .|> upreferred, .-res .|> imag, cmap="Blues", norm=fnorm)
cb = fig.colorbar(cb, ax=ax, pad=0.12, shrink=0.7, label=L"-mathrm{Im}[\chi_{S_zS_z}(\mathbf q^*,\omega)]/N_0")
pos = cb.ax.get_position()
# cb.ax.set_position([pos.x0, pos.y0 + 0.1, pos.width, pos.height])  # Move up & shrink

# cb.ax.set_title(L"-\chi_{S_zS_z}(\mathbf q^*,\omega)/N_0", fontsize=6)
# ax.plot(θrange, [ωplus(q .* [cos(θ), sin(θ), 0], 1, p) * ħ / p.Ef |> upreferred for θ in θrange], color="white", linestyle="dashed", lw=0.5)
# ax.plot(θrange, [ωplus(q .* [cos(θ), sin(θ), 0], -1, p) * ħ / p.Ef |> upreferred for θ in θrange], color="white", linestyle="dashed", lw=0.5)
ax.set_ylim(extrema(ωrange * ħ / p.Ef)...)



ellipsdown = 0.01 ./ sqrt.(1 .- (eccentricity(-1, p) .* cos.(θrange)) .^ 2)
ellipsup = 0.01 ./ sqrt.(1 .- (eccentricity(-1, p) .* sin.(θrange)) .^ 2)

ax.plot(θrange, ellipsdown, color="red", alpha=0.5)
ax.fill_between(θrange, 0, ellipsdown, alpha=0.3, color="red")

ax.plot(θrange, ellipsup, color="blue", alpha=0.5)
ax.fill_between(θrange, 0, ellipsup, alpha=0.3, color="blue")


NN = 100

ax.plot(range(-pi / 4, stop=pi / 4, length=NN), ones(NN) * maximum(ωrange * ħ / p.Ef) * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)
ax.plot(range(3pi / 4, stop=5pi / 4, length=NN), ones(NN) * maximum(ωrange * ħ / p.Ef) * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)

ax.plot(range(pi / 4, stop=3pi / 4, length=NN), ones(NN) * maximum(ωrange * ħ / p.Ef) * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)
ax.plot(range(5pi / 4, stop=7pi / 4, length=NN), ones(NN) * maximum(ωrange * ħ / p.Ef) * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)


ax.text(0, 0.165, L"x", ha="center", va="center")
ax.text(pi / 2, 0.165, L"y", ha="center", va="center")

fig.savefig(joinpath(save_dir, "polar-plot-chiSzSz.pdf"))

fig

##
fq = q -> [q, 0q, 0q]

Nω = 300
ωrange = range(0.001p.Ef / ħ, stop=0.3p.Ef / ħ, length=Nω)
rescale = 1
res = rescale .* [dielectric(ω, fq(q), p) for ω in ωrange]

fig, ax = plt.subplots(figsize=(fig_width, 0.8fig_height))

ax.axvspan(ωrange[1] .* ħ ./ p.Ef |> upreferred, ωplus(fq(q), -1, p) * ħ / p.Ef |> upreferred, facecolor="blue",
    # hatch="//\\\\", 
    alpha=0.15, edgecolor="black")
ax.axvspan(ωplus(fq(q), -1, p) * ħ / p.Ef |> upreferred, ωplus(fq(q), 1, p) * ħ / p.Ef |> upreferred, facecolor="red",
    # hatch="//", 
    alpha=0.15, edgecolor="black")

ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, res .|> real, label=L"\mathrm{Re}\, \epsilon")

# find the poles
edges = (0.5ωplus(fq(q), -1, p) * ħ / p.Ef |> upreferred, ωplus(fq(q), 1, p) * ħ / p.Ef |> upreferred)
poles = find_zeros(x -> dielectric(x * p.Ef / ħ, fq(q), p) |> real |> upreferred, edges...)

ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, res .|> imag, label=L"\mathrm{Im}\, \epsilon")

ax.axvline(ωplus(fq(q), 1, p) * ħ / p.Ef, linestyle="--", color="black", zorder=-10)
ax.axvline(ωplus(fq(q), -1, p) * ħ / p.Ef, linestyle="--", color="black", zorder=-10)
# ax.axvline(ωs(fq(q), p) * ħ / p.Ef, linestyle="-", color="grey")

# for p in poles
#     ax.scatter(p, 0, color="black", facecolors="none", s=10, alpha=0.8, marker="*")
# end
ax.text(poles[1] - 0.01, -500, L"1")
ax.text(poles[2] + 0.005, -500, L"2")
ax.text(poles[3] - 0.005, -500, L"3")


ax.text(ωplus(fq(q), -1, p) * ħ / p.Ef |> upreferred, 2500, L"\omega_{+\downarrow}", ha="center")
ax.text(ωplus(fq(q), 1, p) * ħ / p.Ef |> upreferred, 2500, L"\omega_{+\uparrow}", ha="center")


ax.axhline(0, color="black")
ax.set_xlabel(L"\hbar\omega/\epsilon_F")
ax.set_ylabel(L"\epsilon(q,\omega)")
ax.legend()

ax.margins(x=0)
fig.subplots_adjust(left=0.2, bottom=0.22)
fig.savefig(joinpath(save_dir, "linecut-epsilon_3D.pdf"))
fig

##
Nθ = 101
fig, ax = plt.subplots(figsize=(fig_width, fig_height), subplot_kw=Dict("projection" => "polar"))
fq = θ -> q .* [cos(θ), sin(θ), 0]
@info determine_μ(fq(0), p)
set_rticks_custom(ax)

lines = []
for i in 1:4
    # first divide in 4 segments
    θrange = range(2pi / 4 * (i - 1), stop=2pi / 4 * i, length=Nθ) .- 2pi / 8

    res = [plasmon_pole(fq(θ), p) for θ in θrange]
    res_μ = [determine_μ(fq(θ), p) for θ in θrange]
    x = θrange[.~isnan.(res)]
    y = res[.~isnan.(res)]
    z = res_μ[.~isnan.(res)]

    fnorm = plt.matplotlib.colors.CenteredNorm(0, halfrange=maximum(abs.(z)))
    points = reshape([x'; y']', :, 1, 2)
    segments = cat(points[1:end-1, :, :], points[2:end, :, :], dims=2)

    lc = plt.matplotlib.collections.LineCollection(segments, cmap="seismic", norm=fnorm, linewidth=2)
    lc.set_array(z)
    ax.add_collection(lc)
    # if determine_μ(fq(θrange[Nθ/2|>Int]), p) > 0
    #     color = "C0"
    # else
    #     color = "C1"
    # end
    # ax.plot(θrange, res .|> abs, color=color)
    # ax
    push!(lines, lc)
end
ax.set_rlim(0.001, 0.15)
# ax.set_rticks([0.0, 0.1])
cb = fig.colorbar(lines[1], ax=ax, shrink=0.7, pad=0.12, label=L"\mu_d(\mathbf q^*,\omega_d)/\mu_B")
# cb.ax.set_title(L"\mu_d(\mathbf q^*,\omega_d)/\mu_B")

pos = cb.ax.get_position()
# cb.ax.set_position([pos.x0, pos.y0 + 0.1, pos.width, pos.height])  # Move up & shrink

ax.plot(range(-pi / 4, stop=pi / 4, length=NN), ones(NN) * 0.15 * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)
ax.plot(range(3pi / 4, stop=5pi / 4, length=NN), ones(NN) * 0.15 * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)

ax.plot(range(pi / 4, stop=3pi / 4, length=NN), ones(NN) * 0.15 * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)
ax.plot(range(5pi / 4, stop=7pi / 4, length=NN), ones(NN) * 0.15 * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)
# ax.text(pi / 2, 0.02, L"\omega_p")
ax.text(0, 0.165, L"x", ha="center", va="center")
ax.text(pi / 2, 0.165, L"y", ha="center", va="center")
ax.set_xlabel(L"\hbar\omega_d/E_F")
fig.savefig(joinpath(save_dir, "polar-plot-sign.pdf"))

fig


