using Roots
using IntervalRootFinding
include("style.jl")

fapprox = δq -> 1 + 2exp(-4) * ((1 + δq) / (1 - δq))^δq
##

f = (x, ratio) -> real(8 - 2log((1 + x) / (x - 1)) - 2ratio * x * log((1 + x * ratio) / (x * ratio - 1)))
vs_range = 0.5:0.01:2.0
fig, ax = plt.subplots()
ax.plot(vs_range, [f(x + 1e-8im, 0.8) for x in vs_range])
ax.axhline(0)
fig
##

find_zeros(x -> f(x + 1e-4im, 1), 1.0, 2.0)

##

qratio_range = 0.1:0.001:0.7
res = [find_zeros(x -> f(x + 1e-8im, q), 1.0, 2.0)[1] < q^-1 ? find_zeros(x -> f(x + 1e-8im, q), 1.0, 2.0)[1] : NaN for q in qratio_range]
res2 = [find_zeros(x -> f(x + 1e-8im, q), 1.0, 2.0)[1] > q^-1 ? find_zeros(x -> f(x + 1e-8im, q), 1.0, 2.0)[1] : NaN for q in qratio_range]
@info res[1], fapprox(qratio_range[1])
fig, (ax, ax2) = plt.subplots(2, 1, sharex=true)
ax.plot(qratio_range, res, label="Numerical")
# ax.plot(qratio_range, res2, linestyle="--")
ax.plot(qratio_range[.~isnan.(res)], fapprox.(qratio_range[.~isnan.(res)]), "--", label="Approximation")
# ax2 = ax.twinx()

# g = x -> x / (x^2 - 1) - 1 / 2 * log((1 + x) / (x - 1))
function g(x, y)
    δ = y + 1e-6im
    real(x / (x^2 - 1) - 1 / 2 * log((1 + x) / (x - 1)) +
         x * δ^2 / (x^2 * δ^2 - 1) - 1 / 2 * δ * log((1 + δ * x) / (δ * x - 1)))
end
function g(deltaq)
    g(find_zeros(x -> f(x + 1e-8im, deltaq), 1.0, 2.0)[1], deltaq)
end
function Q(x, y)
    2g(x, y) / pi / y
end

function Q(deltaq)
    Q(find_zeros(x -> f(x + 1e-8im, deltaq), 1.0, 2.0)[1], deltaq)
end

ax2.plot(qratio_range, [Q(x, y) for (x, y) in zip(res, qratio_range)], color="C2")
@info g(0.19)
@info Q(0.19)
# ax.plot(qratio_range, qratio_range.^-1)

ax2.set_ylabel(L"Q")
ax2.set_xlabel(L"\delta_q")
ax.set_ylabel(L"\tilde v_d")
ax.legend()
# ax.set_ylim(extrema(res)...)
ax.margins(x=0)
fig.subplots_adjust(bottom=0.2, right=0.8)
fig.savefig(joinpath(save_dir, "vs_poles.pdf"))

fig


##
m0 = 1
mstar = 1.05
mup = sqrt(1 / m0 - 1 / mstar)
mdown = sqrt(1 / m0 + 1 / mstar)

delta_to_angle = @. x -> acos(
    (mdown - mup * x) / sqrt(
        -4mdown * mup * x + (mdown^2 + mup^2) * (1 + x^2)
    )
) |> rad2deg
angle_to_delta = x -> (mup * cos.(deg2rad.(x)) .+ mdown * sin.(deg2rad.(x))) ./ (mup * sin.(deg2rad.(x)) .+ mdown * cos.(deg2rad.(x)))

delta_to_mstar = @. x -> (1 + x^2) / (1 - x^2)
mstar_to_delta = @. x -> sqrt((1 - 1 / x) / (1 + 1 / x))

qratio_range = angle_to_delta(0):0.001:0.7
res = [find_zeros(x -> f(x + 1e-8im, q), 1.0, 2.0)[1] < q^-1 ? find_zeros(x -> f(x + 1e-8im, q), 1.0, 2.0)[1] : NaN for q in qratio_range]
res2 = [find_zeros(x -> f(x + 1e-8im, q), 1.0, 2.0)[1] > q^-1 ? find_zeros(x -> f(x + 1e-8im, q), 1.0, 2.0)[1] : NaN for q in qratio_range]
@info res[1], fapprox(qratio_range[1])
fig, ax = plt.subplots(1, 1, sharex=true)
ax.plot(qratio_range, res, label="Numerical")
# ax.plot(qratio_range, res2, linestyle="--")
# ax.plot(qratio_range[.~isnan.(res)], fapprox.(qratio_range[.~isnan.(res)]), "--", label="Approximation")
# ax2 = ax.twinx()

# g = x -> x / (x^2 - 1) - 1 / 2 * log((1 + x) / (x - 1))
function g(x, y)
    δ = y + 1e-6im
    real(x / (x^2 - 1) - 1 / 2 * log((1 + x) / (x - 1)) +
         x * δ^2 / (x^2 * δ^2 - 1) - 1 / 2 * δ * log((1 + δ * x) / (δ * x - 1)))
end
function g(deltaq)
    g(find_zeros(x -> f(x + 1e-8im, deltaq), 1.0, 2.0)[1], deltaq)
end
function Q(x, y)
    2g(x, y) / pi / y
end

function Q(deltaq)
    Q(find_zeros(x -> f(x + 1e-8im, deltaq), 1.0, 2.0)[1], deltaq)
end
ax2 = ax.twinx()
ax2.plot(qratio_range, [Q(x, y) for (x, y) in zip(res, qratio_range)], color="C2")
@info g(0.19)
@info Q(0.19)
# ax.plot(qratio_range, qratio_range.^-1)

ax2.set_ylabel(L"Q")
ax2.set_xlabel(L"\delta_q")
ax.set_ylabel(L"v_d/v_F")
# ax.legend()
ax.set_xlabel(L"\delta_q")
twinx = ax.secondary_xaxis("top", functions=(delta_to_angle, angle_to_delta))
twinx.set_xlabel(L"\theta")

# twinx2 = ax.secondary_xaxis(1.2, functions=(identity, identity))

# twinx2 = ax.secondary_xaxis(1.3, functions=(delta_to_mstar, mstar_to_delta))
# twinx2.set_xlabel(L"m^*")

# twinx2.spines["top"].set_position(("axes", 1.2))
# twinx.set_xticks(0:10:30)

# ax.set_ylim(extrema(res)...)
ax.margins(x=0)
fig.subplots_adjust(bottom=0.2, right=0.8, top=0.7, left=0.2)
fig.savefig(joinpath(save_dir, "vs_and_Q_main.pdf"))

fig

##
