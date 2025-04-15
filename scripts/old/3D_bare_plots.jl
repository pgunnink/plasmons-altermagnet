using PyPlot
using ProgressMeter
using Unitful
using Revise
using LsqFit
includet("bare_3D.jl")
includet("style.jl")

##
p = ParamsPlasmon()
Nq = 2 * 100
Nω = 2 * 80
qrange = range(0.01kF(p), stop=0.1kF(p), length=Nq)
ωrange = range(0.005p.Ef / ħ, stop=0.3p.Ef / ħ, length=Nω)

fq = q -> [q, 0q, 0q]

res = zeros(typeof(1.0im * u"meV^-1 * nm^-3"), Nω, Nq)
for (i, ω) in enumerate(ωrange)
    for (j, q) in enumerate(qrange)
        res[i, j] = χSzSz(ω, fq(q), p)
        # res[i, j] = χnSz(ω, fq(q), p)
    end
end

fig, ax = plt.subplots(figsize=(fig_width, fig_height))
toplot = -res .|> imag .|> ustrip
fnorm = plt.matplotlib.colors.LogNorm(1e-2, maximum(toplot))
fnorm = plt.matplotlib.colors.Normalize(0, maximum(toplot))

cb = ax.pcolormesh(qrange ./ kF(p) .|> upreferred, ωrange .* ħ ./ p.Ef .|> upreferred, toplot, norm=fnorm, cmap="Blues")
ax.plot(qrange ./ kF(p) .|> upreferred, [ωplus(fq(q), -1, p) * ħ / p.Ef |> upreferred for q in qrange], color="black", lw=0.5)
ax.plot(qrange ./ kF(p) .|> upreferred, [ωplus(fq(q), 1, p) * ħ / p.Ef |> upreferred for q in qrange], color="black", lw=0.5)
# ax.plot(qrange ./ kF(p) .|> upreferred, [ωs(fq(q), p) * ħ / p.Ef |> upreferred for q in qrange], color="black", linestyle="--", lw=0.5)
fig.colorbar(cb, ax=ax, label=L"-\chi_{S_zS_z}/N_0")
ax.set_xlabel(L"q/k_F")
ax.set_ylabel(L"\hbar\omega/E_F")
ax.set_ylim(extrema(ωrange .* ħ ./ p.Ef .|> upreferred)...)
# ax.plot(qrange ./ kF(p) .|> upreferred, vF(p) .* qrange .* ħ ./ p.Ef, color="black", linestyle="dashed", lw=0.5)
fig.savefig(joinpath(save_dir, "SzSz_3D.pdf"))
fig


##
ωs = x -> vF(p) * sqrt(mDOS(p) / mx(-1, p)) * norm(x)
p = ParamsPlasmon(mstar=1.1me)
Nω = 200
ωrange = range(0.0001p.Ef / ħ, stop=0.004p.Ef / ħ, length=Nω)

q = 0.001kF(p)
fig, ax = plt.subplots(figsize=(fig_width, fig_height))
resup = [-((χ0(ω, fq(q), 1, p))) for ω in ωrange] .|> upreferred


ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, resup .|> ustrip .|> real, label=L"\chi_\uparrow")
resdown = [-((χ0(ω, fq(q), -1, p))) for ω in ωrange] .|> upreferred
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, resdown .|> ustrip .|> real, label=L"\chi_\downarrow")
resdown_approx = [-((n0(p) * qprime_squared(fq(q), -1, p) / (mDOS(p) * ω^2)) * (1 + 3 / 5 * qprime_squared(fq(q), -1, p) * vF(p)^2 / ω^2)) for ω in ωrange] .|> upreferred
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, resdown_approx .|> ustrip .|> real, label=L"\approx \chi_\downarrow")

ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, (resdown_approx .+ N0(p)) .|> ustrip .|> real, label=L"\approx \chi")

ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, (resdown .+ resup) .|> ustrip .|> real, label=L"\chi")

ax.axhline(vq(fq(q), p) |> upreferred |> ustrip)
ax.axvline(ωplus(fq(q), 1, p) * ħ / p.Ef, linestyle="--")
ax.axvline(ωplus(fq(q), -1, p) * ħ / p.Ef, linestyle="--")
# ax.axvline(ωs(fq(q)) * ħ / p.Ef, linestyle="-", color="grey")
ax.axhline(0, color="black")
plt.legend()
ax.set_xlabel(L"\hbar\omega/E_F")
ax.set_ylabel(L"\epsilon")
ax.set_ylim(extrema(resdown .|> ustrip .|> real)...)

fig


##
p = ParamsPlasmon()
Nω = 200
q = 1e-2kF(p)

ωrange = range(0.2ωplus(fq(q), -1, p), stop=1.5ωplus(fq(q), -1, p), length=Nω)
res = [χSzSz(ω, fq(q), p) / N0(p) |> upreferred for ω in ωrange]

fig, ax = plt.subplots(figsize=(fig_width, fig_height))
ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, res .|> imag)

ax.plot(ωrange .* ħ ./ p.Ef .|> upreferred, res .|> real)
ax.axhline(vq(fq(q), p) |> upreferred |> ustrip)
ax.axvline(ωplus(fq(q), 1, p) * ħ / p.Ef, linestyle="--")
# ax.axvline(ωmin(fq(q), -1, p) * ħ / p.Ef |> upreferred, linestyle="--")
# ax.axvline(ωs(fq(q), p) * ħ / p.Ef, linestyle="-", color="grey")
ax.axhline(0, color="black")
ax.set_xlabel(L"\hbar\omega/E_F")
ax.set_ylabel(L"\epsilon")

fig
