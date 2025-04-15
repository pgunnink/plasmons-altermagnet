using DrWatson
@quickactivate
using PyPlot
includet(scriptsdir("style.jl"))

filename = "/Users/pgunnink/backup_onedrive/Projects/plasmons-altermagnet/data/bulk_hbarm0_natural=1.0_vq_natural=82.9_Δ_natural=0.5_ϵ=1.jld2"
##
data = wload(filename)
res = data["res"]
θrange = data["θrange"]
ωrange = data["ωrange"]
toplot = -res .|> imag

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

ax.plot(range(0, stop=pi / 4, length=NN), ones(NN) * maximum(ωrange) * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)
ax.plot(range(pi / 2, stop=3pi / 4, length=NN), ones(NN) * maximum(ωrange) * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)

ax.plot(range(pi, stop=pi + pi / 4, length=NN), ones(NN) * maximum(ωrange) * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)

ax.plot(range(3pi / 2, stop=3pi / 2 + pi / 4, length=NN), ones(NN) * maximum(ωrange) * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)


ax.plot(range(pi / 4, stop=pi / 2, length=NN), ones(NN) * maximum(ωrange) * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)
ax.plot(range(pi - pi / 4, stop=pi, length=NN), ones(NN) * maximum(ωrange) * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)
ax.plot(range(3pi / 2 - pi / 4, stop=3pi / 2, length=NN), ones(NN) * maximum(ωrange) * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)
ax.plot(range(2pi - pi / 4, stop=2pi, length=NN), ones(NN) * maximum(ωrange) * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)


fig.savefig(joinpath(save_dir, "g-wave-polar-plot-chiSzSz.pdf"))
fig