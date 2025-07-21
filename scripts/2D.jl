using DrWatson
using PyPlot
using ProgressMeter
using Unitful

includet(scriptsdir("functions.jl"))
includet(scriptsdir("style.jl"))


function get_conversions(p)
    # m0 = 1
    # mstar = p.mstar / p.m0
    mup = sqrt(1 / p.m0 - 1 / p.mstar)
    mdown = sqrt(1 / p.m0 + 1 / p.mstar)
    delta_to_angle = @. x -> acos(
        (mdown - mup * x) / sqrt(
            -4mdown * mup * x + (mdown^2 + mup^2) * (1 + x^2)
        )
    ) |> rad2deg
    angle_to_delta = x -> (mup * cos.(deg2rad.(x)) .+ mdown * sin.(deg2rad.(x))) ./ (mup * sin.(deg2rad.(x)) .+ mdown * cos.(deg2rad.(x)))
    return angle_to_delta, delta_to_angle
end
##
p = ParamsPlasmon()
q = 0.05kF(p)
@info splitting(p)
Nq = 2 * 100
Nω = 2 * 80
qrange = range(0.001kF(p), stop=0.3kF(p), length=Nq)
ωrange = range(0.001p.Ef / ħ, stop=1p.Ef / ħ, length=Nω)
##
fq = q -> [0q, q]

res = zeros(typeof(1.0im), Nω, Nq)
for (i, ω) in enumerate(ωrange)
    for (j, q) in enumerate(qrange)
        res[i, j] = χSzSz(ω, fq(q), p) / N0(p) |> upreferred
    end
end

fig, ax = plt.subplots(figsize=(fig_width, fig_height))
toplot = -res .|> imag .|> ustrip
# fnorm = plt.matplotlib.colors.LogNorm(1e-2, maximum(toplot))
fnorm = plt.matplotlib.colors.Normalize(0, maximum(toplot))

cb = ax.pcolormesh(qrange ./ kF(p) .|> upreferred, ωrange .* ħ ./ p.Ef .|> upreferred, toplot, norm=fnorm, cmap="Blues")
ax.plot(qrange ./ kF(p) .|> upreferred, [ωplus(fq(q), -1, p) * ħ / p.Ef |> upreferred for q in qrange], color="black", lw=1, alpha=0.5, linestyle="--")
ax.plot(qrange ./ kF(p) .|> upreferred, [ωplus(fq(q), 1, p) * ħ / p.Ef |> upreferred for q in qrange], color="black", lw=1, alpha=0.5, linestyle="--")
# ax.plot(qrange ./ kF(p) .|> upreferred, [ωs(fq(q), p) * ħ / p.Ef |> upreferred for q in qrange], color="black", linestyle="--", lw=0.5)
fig.colorbar(cb, ax=ax, label=L"-\chi_{S_zS_z}/N_0")
ax.set_xlabel(L"q/k_F")
ax.set_ylabel(L"\hbar\omega/E_F")
ax.set_ylim(extrema(ωrange .* ħ ./ p.Ef .|> upreferred)...)
fig.savefig(joinpath(save_dir, "SzSz_2D.pdf"))

fig
##
θ = deg2rad(0)
fq = q -> q .* [cos(θ), sin(θ)]

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
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, res .|> imag, label=L"\mathrm{Im}\, \epsilon")

ax.axvline(ωplus(fq(q), 1, p) * ħ / p.Ef, linestyle="--", color="black", zorder=-10)
ax.axvline(ωplus(fq(q), -1, p) * ħ / p.Ef, linestyle="--", color="black", zorder=-10)



ax.text(ωplus(fq(q), -1, p) * ħ / p.Ef |> upreferred, 2500, L"\omega_{+\downarrow}", ha="center")
ax.text(ωplus(fq(q), 1, p) * ħ / p.Ef |> upreferred, 2500, L"\omega_{+\uparrow}", ha="center")


ax.axhline(0, color="black")
ax.set_xlabel(L"\hbar\omega/\epsilon_F")
ax.set_ylabel(L"\epsilon(q,\omega)")
ax.legend()

ax.margins(x=0)
fig.subplots_adjust(left=0.2, bottom=0.22)
fig.savefig(joinpath(save_dir, "linecut-epsilon_2D.pdf"))
fig

##
Nθ = 500
Nω = 500
ωrange = range(0.001p.Ef / ħ, stop=0.15p.Ef / ħ, length=Nω)

θrange = range(0, stop=2pi, length=Nθ)
res = zeros(typeof(1.0im), Nω, Nθ)
for (i, ω) in enumerate(ωrange)
    for (j, θ) in enumerate(θrange)
        res[i, j] = χSzSz(ω, q .* [cos(θ), sin(θ)], p) / N0(p) |> upreferred
        # res[i, j] = χnn(ω, fq(q), p) / N0(p) |> upreferred
    end
end
##
fnorm = plt.matplotlib.colors.Normalize(0, maximum(toplot))
fig, ax = plt.subplots(figsize=(fig_width, fig_height), subplot_kw=Dict("projection" => "polar"))
ax.set_rlabel_position(0 * 45 / 2)
ax.set_rticks([0.0, 0.05, 0.1,])
# ax.set_xticks([0, 90, 180, 270] / 360 * 2pi)
cb = ax.pcolormesh(θrange, ωrange .* ħ ./ p.Ef .|> upreferred, .-res .|> imag, cmap="Blues", norm=fnorm)
cb = fig.colorbar(cb, ax=ax, pad=0.12, shrink=0.7, panchor=(1.0, 0.8), label=L"-\chi_{S_zS_z}(\mathbf q^*,\omega)/N_0")

# pos = cb.ax.get_position()
# cb.ax.set_position([pos.x0, pos.y0 + 0.1, pos.width, pos.height])  # Move up & shrink
# cb.ax.set_title(L"-\chi_{S_zS_z}(\mathbf q^*,\omega)/N_0", fontsize=6)

# ax.plot(θrange, [ωplus(q .* [cos(θ), sin(θ), 0], 1, p) * ħ / p.Ef |> upreferred for θ in θrange], color="white", linestyle="dashed", lw=0.5)
# ax.plot(θrange, [ωplus(q .* [cos(θ), sin(θ), 0], -1, p) * ħ / p.Ef |> upreferred for θ in θrange], color="white", linestyle="dashed", lw=0.5)
set_rticks_custom(ax)
ax.set_ylim(extrema(ωrange * ħ / p.Ef)...)
NN = 100

ax.plot(range(-pi / 4, stop=pi / 4, length=NN), ones(NN) * maximum(ωrange * ħ / p.Ef) * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)
ax.plot(range(3pi / 4, stop=5pi / 4, length=NN), ones(NN) * maximum(ωrange * ħ / p.Ef) * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)

ax.plot(range(pi / 4, stop=3pi / 4, length=NN), ones(NN) * maximum(ωrange * ħ / p.Ef) * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)
ax.plot(range(5pi / 4, stop=7pi / 4, length=NN), ones(NN) * maximum(ωrange * ħ / p.Ef) * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)

ax.text(0, 0.165, L"x", ha="center", va="center")
ax.text(pi / 2, 0.165, L"y", ha="center", va="center")

fig.savefig(joinpath(save_dir, "2D_polar-plot-chiSzSz.pdf"))

θrange = range(0, stop=2pi, length=100)
ax.plot(θrange, [ħ * q * vs(q, θ, p) / p.Ef |> upreferred |> ustrip for θ in θrange], alpha=0.5, color="red")



fig.savefig(joinpath(save_dir, "2D_polar-plot-chiSzSz-with-vs.pdf"))


fig

##

N = 100
θrange = range(0, stop=pi / 4, length=N)

fig, axs = plt.subplots(2, 1, sharex=true)
ax.margins(x=0)
angle_to_delta, delta_to_angle = get_conversions(p)

twinx = ax.secondary_xaxis("top", functions=(angle_to_delta, delta_to_angle))
# twinx.set_xticks([0.2, 0.3])
twinx.set_xlabel(L"\eta_{\mathrm{maj}}(\theta)/\eta_{\mathrm{min}}(\theta)")



ω_N = [angle_to_delta(rad2deg(θ)) < sqrt(3) / 2 ? vs(θ, p) / vF(p) : NaN for θ in θrange] .|> upreferred
Q_N = [Q(θ, p) for θ in θrange] .|> upreferred
axs[1].plot(θrange .|> rad2deg, ω_N)
axs[2].plot(θrange .|> rad2deg, Q_N, color="C1")
axs[1].margins(x=0)


angle_to_delta, delta_to_angle = get_conversions(p)

fig.subplots_adjust(hspace=0.2, top=0.8, bottom=0.2, left=0.15)
fig.savefig(joinpath(save_dir, "2D_angles_vs_and_Q_no_eta.pdf"))

twinx = axs[1].secondary_xaxis("top", functions=(angle_to_delta, delta_to_angle))
# twinx.set_xticks([0.2, 0.3])
twinx.set_xlabel(L"\eta_{\mathrm{maj}}(\theta)/\eta_{\mathrm{min}}(\theta)")
axs[2].set_xlabel(L"\theta")
axs[1].set_ylabel(L"v_d/v_F")
axs[2].set_ylabel(L"Q")


twinx.set_xlim(extrema(angle_to_delta.(θrange .|> rad2deg))...)
twinx.set_xticks(0.4:0.1:1.0)
axs[1].xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter((x, _) -> L"%$(round(x) |> Int)^{\circ}"))

fig.savefig(joinpath(save_dir, "2D_angles_vs_and_Q.pdf"))

fig

##

Nθ = 100
fig, ax = plt.subplots(figsize=(fig_width, fig_height), subplot_kw=Dict("projection" => "polar"))
@info μ(q, 0, p) / Unitful.μB |> upreferred
lines = []
for i in 1:4
    # first divide in 4 segments
    θrange = range(2pi / 4 * (i - 1), stop=2pi / 4 * i, length=Nθ) .- 2pi / 8

    res = [ħ * q * vs(θ, p) / p.Ef for θ in θrange]
    res_μ = [μ(q, θ, p) / Unitful.μB |> upreferred for θ in θrange]
    x = θrange[.~isnan.(res)]
    y = res[.~isnan.(res)]
    z = res_μ[.~isnan.(res)]
    fnorm = plt.matplotlib.colors.CenteredNorm(0, halfrange=maximum(abs.(z)))
    @info maximum(abs.(z))
    points = reshape([x'; y']', :, 1, 2)
    segments = cat(points[1:end-1, :, :], points[2:end, :, :], dims=2)

    lc = plt.matplotlib.collections.LineCollection(segments, cmap="seismic", norm=fnorm)
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
ax.set_rlabel_position(-60)
# ax.set_rticks([0.0, 0.05, 0.1, 0.15])
# ax.set_rticks([0.0, 0.2, 0.4])
# ax.set_ylabel(L"\mu", rotation=0, size=11)
ax.set_rlim(0.0, 0.15)
cb = fig.colorbar(lines[1], ax=ax, shrink=0.7, pad=0.13, label=L"\mu_d(\mathbf q^*,\omega_d)/\mu_B")



set_rticks_custom(ax)
cb.ax.set_title(L"\mu/\mu_B")

ax.plot(range(-pi / 4, stop=pi / 4, length=NN), ones(NN) * 0.15 * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)
ax.plot(range(3pi / 4, stop=5pi / 4, length=NN), ones(NN) * 0.15 * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)

ax.plot(range(pi / 4, stop=3pi / 4, length=NN), ones(NN) * 0.15 * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)
ax.plot(range(5pi / 4, stop=7pi / 4, length=NN), ones(NN) * 0.15 * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)
set_rticks_custom

# ax.text(pi / 2, 0.02, L"\omega_p")
fig.savefig(joinpath(save_dir, "2D-polar-plot-sign.pdf"))

fig

