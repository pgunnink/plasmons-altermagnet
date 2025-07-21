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
cb = ax.pcolormesh(θrange, ωrange .* ħ ./ p.Ef .|> upreferred, .-res .|> imag, cmap="Blues", norm=fnorm, rasterized=true)

ax.set_ylim(extrema(ωrange * ħ / p.Ef)...)



ellipsdown = 0.01 ./ sqrt.(1 .- (eccentricity(-1, p) .* cos.(θrange)) .^ 2)
ellipsup = 0.01 ./ sqrt.(1 .- (eccentricity(-1, p) .* sin.(θrange)) .^ 2)

# ax.plot(θrange, ellipsdown, color="red", alpha=0.5)
# ax.fill_between(θrange, 0, ellipsdown, alpha=0.3, color="red")

# ax.plot(θrange, ellipsup, color="blue", alpha=0.5)
# ax.fill_between(θrange, 0, ellipsup, alpha=0.3, color="blue")


NN = 100

ax.plot(range(-pi / 4, stop=pi / 4, length=NN), ones(NN) * maximum(ωrange * ħ / p.Ef) * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)
ax.plot(range(3pi / 4, stop=5pi / 4, length=NN), ones(NN) * maximum(ωrange * ħ / p.Ef) * 0.98, lw=2, color="red", solid_capstyle="butt", alpha=0.5)

ax.plot(range(pi / 4, stop=3pi / 4, length=NN), ones(NN) * maximum(ωrange * ħ / p.Ef) * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)
ax.plot(range(5pi / 4, stop=7pi / 4, length=NN), ones(NN) * maximum(ωrange * ħ / p.Ef) * 0.98, lw=2, color="blue", solid_capstyle="butt", alpha=0.5)


ax.text(0, 0.165, L"x", ha="center", va="center")
ax.text(pi / 2, 0.165, L"y", ha="center", va="center")

fig.savefig(joinpath(save_dir, "emmy-polar-plot-chiSzSz.pdf"))

fig