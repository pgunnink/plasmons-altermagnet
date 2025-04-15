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
N = 100
θrange = range(0, stop=pi / 4, length=N)
q = 0.1kF(p)
frot(θ) = q .* [cos(θ), sin(θ), 0]

ωp_N = [plasmon_pole(frot(θ), p) * p.Ef / ħ / (q * vF(p)) for θ in θrange] .|> upreferred
Q_N = [Q_factor(frot(θ), p) for θ in θrange]

ωp_a = [analytical_plasmon_pole(δq(frot(θ), p)) * (qprime(frot(θ), -1, p) / q) for θ in θrange]
Q_a = [Q(δq(frot(θ), p)) for θ in θrange]

fig, (ax, ax2) = plt.subplots(2, 1, sharex=true)
ax.plot(θrange .|> rad2deg, ωp_N, color="C0", linestyle="-", label="Numerical")
ax.plot(θrange .|> rad2deg, ωp_a, color="C0", linestyle="--", label="Analytical")
# ax.set_xlim(0, )
ax.legend()
ax.margins(x=0)
angle_to_delta, delta_to_angle = get_conversions(p)



# ax2 = ax.twinx()
ax2.plot(θrange .|> rad2deg, Q_N, color="C1", linestyle="-")
ax2.plot(θrange .|> rad2deg, Q_a, color="C1", linestyle="--")
ax2.set_xlabel(L"\theta")
ax.set_ylabel(L"v_s/v_F")
ax2.set_ylabel(L"Q")
# twinx.set_xticks(twinx.get_xticks())


fig.subplots_adjust(hspace=0.2, top=0.8, bottom=0.2, left=0.15)

fig.savefig(joinpath(save_dir, "angles_vs_and_Q.pdf"))


twinx = ax.secondary_xaxis("top", functions=(angle_to_delta, delta_to_angle))
# twinx.set_xticks([0.2, 0.3])
twinx.set_xlabel(L"\eta_{\mathrm{min}}(\theta)/\eta_{\mathrm{maj}}(\theta)")

twinx.set_xlim(extrema(angle_to_delta.(θrange .|> rad2deg))...)
twinx.set_xticks(0.4:0.1:1.0)
ax2.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter((x, _) -> L"%$x^{\circ}"))

# ax2.margins(x=0)
# ax.set_xlim(0, 5)
fig.savefig(joinpath(save_dir, "angles_vs_and_Q_with_eta.pdf"))

fig