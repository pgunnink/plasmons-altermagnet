include("functions_3D.jl")
using Roots

## conversion
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
## numerical functions
function plasmon_pole(q, p)
    edges = (ωplus(q, -1, p) * ħ / p.Ef |> upreferred, ωplus(q, 1, p) * ħ / p.Ef |> upreferred)
    poles = find_zeros(x -> dielectric(x * p.Ef / ħ, q, p) |> real |> upreferred, edges...)
    if length(poles) > 1
        return poles[1]
    else
        return NaN
    end
end

##
function Q_factor(q, ω, p)
    dω = 1e-4
    dω_pole = (dielectric((ω - dω) * p.Ef / ħ, q, p) - dielectric((ω + dω) * p.Ef / ħ, q, p)) / 2dω |> real |> upreferred
    γ = imag(dielectric((ω) * p.Ef / ħ, q, p)) / dω_pole
    -ω / γ |> upreferred
end

Q_factor(q, p) = Q_factor(q, plasmon_pole(q, p), p)
##
function determine_μ(q, p)
    edges = (ωplus(q, -1, p) * ħ / p.Ef |> upreferred, ωplus(q, 1, p) * ħ / p.Ef |> upreferred)
    probe_T = 0.1u"T"
    local_Δ = Δ(probe_T, p)
    custom_dielectric = (δ, ω, q, p) -> 1 - vq(q, p) * ((1 + δ) * χ0(ω, q, 1, p) + (1 - δ) * χ0(ω, q, -1, p))
    poles_plus = find_zeros(x -> custom_dielectric(local_Δ, x * p.Ef / ħ, q, p) |> real |> upreferred, edges...)
    poles_min = find_zeros(x -> custom_dielectric(-local_Δ, x * p.Ef / ħ, q, p) |> real |> upreferred, edges...)
    if length(poles_plus) > 1 && length(poles_min) > 1
        (poles_plus[1] - poles_min[1]) / 2local_Δ * p.Ef * dΔdB(p) / Unitful.μB |> upreferred
    else
        return NaN
    end
end

## analytical functions




function analytical_plasmon_pole(deltaq)
    f = x -> real(4 - log((1 + x) / (x - 1)) - deltaq * x * log((1 + x * deltaq) / (x * deltaq - 1)))
    res = find_zeros(x -> f(x + 1e-6im), 1.0, 2.0)[1]
    res < deltaq^-1 ? res : NaN
end

function c0(x, y)
    δ = y + 1e-6im
    real(x / (x^2 - 1) - 1 / 2 * log((1 + x) / (x - 1)) +
         x * δ^2 / (x^2 * δ^2 - 1) - 1 / 2 * δ * log((1 + δ * x) / (δ * x - 1)))
end
function c0(deltaq)
    c0(analytical_plasmon_pole(deltaq), deltaq)
end
function Q(x, y)
    2c0(x, y) / pi / y
end

function Q(deltaq)
    Q(analytical_plasmon_pole(deltaq), deltaq)
end
