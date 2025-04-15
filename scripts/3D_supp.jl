using DrWatson
@quickactivate

using PyPlot
using ProgressMeter
using ForwardDiff
using Unitful
using Revise
using LsqFit
includet(scriptsdir("functions_3D.jl"))
includet(scriptsdir("style.jl"))

##
p = ParamsPlasmon()
@info splitting(p)
Nq = 3 * 100
Nω = 3 * 80
qrange = range(0.001kF(p), stop=0.1kF(p), length=Nq)
ωrange = range(0.001p.Ef / ħ, stop=0.3p.Ef / ħ, length=Nω)


fq = q -> [q, 0q, 0q]

res = zeros(typeof(1.0im), Nω, Nq)
for (i, ω) in enumerate(ωrange)
    for (j, q) in enumerate(qrange)
        res[i, j] = χnSz(ω, fq(q), p) / N0(p) |> upreferred
        # res[i, j] = χnn(ω, fq(q), p) / N0(p) |> upreferred
    end
end

fig, ax = plt.subplots(figsize=(fig_width, fig_height))
toplot = res .|> imag .|> ustrip
# toplot[toplot.<1e-5] .= 1e-5
fnorm = plt.matplotlib.colors.LogNorm(1e-1, maximum(toplot))
fnorm = plt.matplotlib.colors.Normalize(0, maximum(toplot))

cb = ax.pcolormesh(qrange ./ kF(p) .|> upreferred, ωrange .* ħ ./ p.Ef .|> upreferred, toplot, norm=fnorm, cmap="Blues", shading="nearest")
ax.plot(qrange ./ kF(p) .|> upreferred, [ωplus(fq(q), -1, p) * ħ / p.Ef |> upreferred for q in qrange], color="black", linestyle="dashed", lw=1, alpha=0.5)
ax.plot(qrange ./ kF(p) .|> upreferred, [ωplus(fq(q), 1, p) * ħ / p.Ef |> upreferred for q in qrange], color="black", linestyle="dashed", lw=1, alpha=0.5)
# ax.plot(qrange ./ kF(p) .|> upreferred, [ωs(fq(q), p) * ħ / p.Ef |> upreferred for q in qrange], color="black", linestyle="--", lw=0.5)
fig.colorbar(cb, ax=ax, label=L"-mathrm{Im}[\chi_{nS_z}(\mathbf q,\omega)]/N_0")
ax.set_xlabel(L"q/k_F")
ax.set_ylabel(L"\hbar\omega/\epsilon_F")
ax.set_ylim(extrema(ωrange .* ħ ./ p.Ef .|> upreferred)...)
# ax.plot(qrange ./ kF(p) .|> upreferred, vF(p) .* qrange .* ħ ./ p.Ef, color="black", linestyle="dashed", lw=0.5)
fig.subplots_adjust(left=0.15, bottom=0.2)
fig.savefig(joinpath(save_dir, "nSz_3D.pdf"))
fig

##
p = ParamsPlasmon()
@info splitting(p)
Nq = 3 * 100
Nω = 3 * 80
qrange = range(0.001kF(p), stop=10 * 0.1kF(p), length=Nq)
ωrange = range(0.001p.Ef / ħ, stop=5 * 0.3p.Ef / ħ, length=Nω)


fq = q -> [q, 0q, 0q]

res = zeros(typeof(1.0im), Nω, Nq)
for (i, ω) in enumerate(ωrange)
    for (j, q) in enumerate(qrange)
        res[i, j] = χSzSz(ω, fq(q), p) / N0(p) |> upreferred
        # res[i, j] = χnn(ω, fq(q), p) / N0(p) |> upreferred
    end
end

fig, ax = plt.subplots(figsize=(fig_width, fig_height))
toplot = -res .|> imag .|> ustrip
# toplot[toplot.<1e-5] .= 1e-5
fnorm = plt.matplotlib.colors.LogNorm(1e-1, maximum(toplot))
fnorm = plt.matplotlib.colors.Normalize(0, maximum(toplot))

cb = ax.pcolormesh(qrange ./ kF(p) .|> upreferred, ωrange .* ħ ./ p.Ef .|> upreferred, toplot, norm=fnorm, cmap="Blues", shading="nearest")
ax.plot(qrange ./ kF(p) .|> upreferred, [ωplus(fq(q), -1, p) * ħ / p.Ef |> upreferred for q in qrange], color="black", linestyle="dashed", lw=1, alpha=0.1)
ax.plot(qrange ./ kF(p) .|> upreferred, [ωplus(fq(q), 1, p) * ħ / p.Ef |> upreferred for q in qrange], color="black", linestyle="dashed", lw=1, alpha=0.1)
# ax.plot(qrange ./ kF(p) .|> upreferred, [ωs(fq(q), p) * ħ / p.Ef |> upreferred for q in qrange], color="black", linestyle="--", lw=0.5)
fig.colorbar(cb, ax=ax, label=L"-mathrm{Im}[\chi_{S_zS_z}(\mathbf q,\omega)]/N_0")
ax.set_xlabel(L"q/k_F")
ax.set_ylabel(L"\hbar\omega/\epsilon_F")
ax.set_ylim(extrema(ωrange .* ħ ./ p.Ef .|> upreferred)...)
# ax.plot(qrange ./ kF(p) .|> upreferred, vF(p) .* qrange .* ħ ./ p.Ef, color="black", linestyle="dashed", lw=0.5)
fig.subplots_adjust(left=0.15, bottom=0.2)
fig.savefig(joinpath(save_dir, "SzSz_3D_large_range.pdf"))
fig


##
p = ParamsPlasmon()
@info splitting(p)
Nq = 3 * 100
Nω = 3 * 80
qrange = range(0.001kF(p), stop=0.1kF(p), length=Nq)
ωrange = range(0.001p.Ef / ħ, stop=0.3p.Ef / ħ, length=Nω)

for θ in [10, 20, 50, 70]
    fq = q -> q .* [cos(θ |> deg2rad), 0, sin(θ |> deg2rad)]

    res = zeros(typeof(1.0im), Nω, Nq)
    for (i, ω) in enumerate(ωrange)
        for (j, q) in enumerate(qrange)
            res[i, j] = χSzSz(ω, fq(q), p) / N0(p) |> upreferred
            # res[i, j] = χnn(ω, fq(q), p) / N0(p) |> upreferred
        end
    end

    fig, ax = plt.subplots(figsize=(fig_width, 0.7fig_height))
    toplot = -res .|> imag .|> ustrip
    toplot[toplot.<1e-5] .= 1e-5
    fnorm = plt.matplotlib.colors.LogNorm(1e-1, maximum(toplot))
    fnorm = plt.matplotlib.colors.Normalize(0, maximum(toplot))

    cb = ax.pcolormesh(qrange ./ kF(p) .|> upreferred, ωrange .* ħ ./ p.Ef .|> upreferred, toplot, norm=fnorm, cmap="Blues", shading="nearest")
    ax.plot(qrange ./ kF(p) .|> upreferred, [ωplus(fq(q), -1, p) * ħ / p.Ef |> upreferred for q in qrange], color="black", linestyle="dashed", lw=1, alpha=0.5)
    ax.plot(qrange ./ kF(p) .|> upreferred, [ωplus(fq(q), 1, p) * ħ / p.Ef |> upreferred for q in qrange], color="black", linestyle="dashed", lw=1, alpha=0.5)
    # ax.plot(qrange ./ kF(p) .|> upreferred, [ωs(fq(q), p) * ħ / p.Ef |> upreferred for q in qrange], color="black", linestyle="--", lw=0.5)
    fig.colorbar(cb, ax=ax, label=L"-mathrm{Im}[\chi_{S_zS_z}(\mathbf q^*,\omega)]/N_0")
    ax.set_xlabel(L"q/k_F")
    ax.set_ylabel(L"\hbar\omega/\epsilon_F")
    ax.set_ylim(extrema(ωrange .* ħ ./ p.Ef .|> upreferred)...)
    # ax.plot(qrange ./ kF(p) .|> upreferred, vF(p) .* qrange .* ħ ./ p.Ef, color="black", linestyle="dashed", lw=0.5)
    fig.subplots_adjust(left=0.15, bottom=0.25)
    fig.savefig(joinpath(save_dir, "SzSz_3D_theta_$(θ).pdf"))
end
# fig



##
p = ParamsPlasmon()
@info splitting(p)
Nq = 3 * 100
Nω = 3 * 80
qrange = range(0.001kF(p), stop=0.1kF(p), length=Nq)
ωrange = range(0.001p.Ef / ħ, stop=0.3p.Ef / ħ, length=Nω)


fq = q -> [q, 0q, 0q]


Nω = 300
ωrange = range(0.001p.Ef / ħ, stop=0.3p.Ef / ħ, length=Nω)
q = 0.05kF(p)
chinnresponse = [χnn(ω, fq(q), p) / N0(p) |> upreferred for ω in ωrange]

fig, ax = plt.subplots(figsize=(fig_width, fig_height))
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, 1e4chinnresponse .|> real, color="C0", label=L"\mathrm{Re}\, \chi_{nn}\times10^4")
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, 1e4chinnresponse .|> imag, color="C0", label=L"\mathrm{Im}\, \chi_{nn}\times10^4", linestyle="dashed")

chinSzresponse = [χnSz(ω, fq(q), p) / N0(p) |> upreferred for ω in ωrange]
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, 1e4chinSzresponse .|> real, color="C1", label=L"\mathrm{Re}\, \chi_{nS_z}\times10^{4}")
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, 1e4chinSzresponse .|> imag, color="C1", label=L"\mathrm{Im}\, \chi_{nS_z}\times10^{4}", linestyle="dashed")

chiSzSzresponse = [χSzSz(ω, fq(q), p) / N0(p) |> upreferred for ω in ωrange]
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, chiSzSzresponse .|> real, color="C2", label=L"\mathrm{Re}\, \chi_{S_zS_z}")
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, chiSzSzresponse .|> imag, color="C2", label=L"\mathrm{Im}\, \chi_{S_zS_z}", linestyle="dashed")

# ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, res .|> imag, label=L"\mathrm{Im}\, \epsilon")

ax.axvline(ωplus(fq(q), 1, p) * ħ / p.Ef, linestyle="--", color="black", zorder=-10)
ax.axvline(ωplus(fq(q), -1, p) * ħ / p.Ef, linestyle="--", color="black", zorder=-10)
ax.axhline(0, color="black")
ax.set_xlabel(L"\hbar\omega/\epsilon_F")
ax.set_ylabel(L"\chi_{\eta}(\mathbf q,\omega)/N_0")
ax.legend(loc="upper right", frameon=true,)
ax.margins(x=0)
fig.subplots_adjust(left=0.2, bottom=0.2)
fig.savefig(joinpath(save_dir, "linecut-chi_nn.pdf"))

fig
##
p = ParamsPlasmon()
@info splitting(p)
Nq = 3 * 100
Nω = 3 * 80
qrange = range(0.001kF(p), stop=0.1kF(p), length=Nq)
ωrange = range(0.001p.Ef / ħ, stop=0.3p.Ef / ħ, length=Nω)


fq = q -> [q, 0q, 0q]


Nω = 300
ωrange = range(0.001p.Ef / ħ, stop=0.3p.Ef / ħ, length=Nω)
q = 0.05kF(p)
chiupupresponse = [χupup(ω, fq(q), p) / N0(p) |> upreferred for ω in ωrange]


fig, ax = plt.subplots(figsize=(fig_width, fig_height))
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, 1e4chiupupresponse .|> real, color="C0", label=L"\mathrm{Re}\, \chi_{upup}\times10^4")
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, 1e4chiupupresponse .|> imag, color="C0", label=L"\mathrm{Im}\, \chi_{upup}\times10^4", linestyle="dashed")

chidowndownresponse = [χdowndown(ω, fq(q), p) / N0(p) |> upreferred for ω in ωrange]
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, 1e4chidowndownresponse .|> real, color="C1", label=L"\mathrm{Re}\, \chi_{nS_z}\times10^{4}")
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, 1e4chidowndownresponse .|> imag, color="C1", label=L"\mathrm{Im}\, \chi_{nS_z}\times10^{4}", linestyle="dashed")

chiupdownresponse = [χupdown(ω, fq(q), p) / N0(p) |> upreferred for ω in ωrange]
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, chiupdownresponse .|> real, color="C2", label=L"\mathrm{Re}\, \chi_{S_zS_z}")
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, chiupdownresponse .|> imag, color="C2", label=L"\mathrm{Im}\, \chi_{S_zS_z}", linestyle="dashed")

# ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, res .|> imag, label=L"\mathrm{Im}\, \epsilon")

ax.axvline(ωplus(fq(q), 1, p) * ħ / p.Ef, linestyle="--", color="black", zorder=-10)
ax.axvline(ωplus(fq(q), -1, p) * ħ / p.Ef, linestyle="--", color="black", zorder=-10)
ax.axhline(0, color="black")
ax.set_xlabel(L"\hbar\omega/\epsilon_F")
ax.set_ylabel(L"\chi_{\eta}(\mathbf q,\omega)/N_0")
ax.legend(loc="upper right", frameon=true,)
ax.margins(x=0)
fig.subplots_adjust(left=0.2, bottom=0.2)
fig.savefig(joinpath(save_dir, "linecut-individual-chi.pdf"))

fig
##

Nθ = 100
fig, ax = plt.subplots(figsize=(fig_width, fig_height),)
fq = θ -> q .* [cos(θ), sin(θ), 0]
@info determine_μ(fq(0), p)
lines = []
for i in 1:2
    # first divide in 4 segments
    θrange = range(2pi / 4 * (i - 1), stop=2pi / 4 * i, length=Nθ)

    res = [plasmon_pole(fq(θ), p) for θ in θrange]
    res_μ = [determine_μ(fq(θ), p) for θ in θrange]
    x = θrange[.~isnan.(res)]
    y = res[.~isnan.(res)]
    z = res_μ[.~isnan.(res)]
    ax.plot(θrange .|> rad2deg, res_μ, color="C0")
end
# ax.set_rlabel_position(0)
# ax.set_rticks([0.0, 0.05, 0.1, 0.15])
# ax.set_rticks([0.0, 0.2, 0.4])
# ax.set_ylabel(L"\mu", rotation=0, size=11)
# ax.set_ylim(0.0, 0.1)
ax.set_xlabel(L"\theta")
ax.set_xticks([0, 45, 90, 135, 180])
ax.set_ylabel(L"\mu/\mu_B")
fig.subplots_adjust(left=0.2, bottom=0.2)
fig.savefig(joinpath(save_dir, "non-polar-plot-sign.pdf"))

fig

