using DrWatson
using PyPlot
using ProgressMeter
using Unitful
includet(scriptsdir("style.jl"))


includet(scriptsdir("functions_3D.jl"))
includet(scriptsdir("functions_analysis.jl"))
p = ParamsPlasmon()

@info splitting(p)


## as a function of rotation
N = 1000
θrange = range(0, stop=2pi, length=N)
q = 0.1kF(p)
frot(θ) = q .* [cos(θ), sin(θ), 0]
x = []
y = []
z = []
for θ in θrange
    qloc = frot(θ)
    edges = (min(0.5ωplus(qloc, -1, p) * ħ / p.Ef |> upreferred, 0.5ωplus(qloc, 1, p) * ħ / p.Ef |> upreferred), 1.5max(ωplus(qloc, -1, p) * ħ / p.Ef |> upreferred, ωplus(qloc, 1, p) * ħ / p.Ef |> upreferred))


    poles = find_zeros(x -> dielectric(x * p.Ef / ħ, qloc, p) |> real |> upreferred, edges...)
    for (i, p) in enumerate(poles)
        push!(x, θ)
        push!(y, p)
        push!(z, "C$(i-1)")
    end
end


##
fig, ax = plt.subplots(figsize=(fig_width, 1.3fig_height), subplot_kw=Dict("projection" => "polar"))

ax.scatter(x, y, c=z, s=1)

ax.plot(θrange, [ωplus(frot(θ), -1, p) * ħ / p.Ef |> upreferred for θ in θrange], color="red", linestyle="--")
ax.plot(θrange, [ωplus(frot(θ), 1, p) * ħ / p.Ef |> upreferred for θ in θrange], color="blue", linestyle="--")

fig.subplots_adjust(top=0.7)


ax.legend(handles=[
        plt.matplotlib.lines.Line2D([], [], color="C0", markersize=2, marker="o", label="1st zero (spin-minority acoustic plasmon)"),
        plt.matplotlib.lines.Line2D([], [], color="C1", markersize=2, marker="o", label="2nd zero (spin demon)"),
        plt.matplotlib.lines.Line2D([], [], color="C2", markersize=2, marker="o", label="3rd zero (spin-majority acoustic plasmon)")
    ],
    bbox_to_anchor=(0.5, 1), loc="upper center",
    bbox_transform=fig.transFigure)

# set_rticks_custom(ax)
ax.grid(color="black", alpha=0.2)
ax.set_rlabel_position(45)
ax.set_rlim(0, 0.3)

ax.text(pi / 4, 0.3, L"\hbar\omega/\epsilon_F", size=7)
ax.set_rticks([0.1, 0.2], [0.1, 0.2])
fig.savefig(joinpath(save_dir, "all_zeros.pdf"))

fig
