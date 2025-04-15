using PyPlot
using ProgressMeter
using Unitful
using Revise
includet("style.jl")


includet("functions_magnons.jl")
p = ParamsMagnons()
# @info "splitting ratio " p.α * kT(p)^2 / p.vs * kT(p) |> upreferred
@info "gap" p.Δ |> u"meV"

@info ωk([0.2pi / p.a, 0 / p.a, 0 / p.a], 1, p)
@info ωk([0.2pi / p.a, 0 / p.a, 0 / p.a], 1, p) - ωk([0.2pi / p.a, 0 / p.a, 0 / p.a], -1, p)
@info "temperature" Unitful.k * p.T |> u"meV"
##

##
Nω = 10

ωrange = range(0.1u"GHz", stop=5u"GHz", length=Nω)

q = 1e-4u"nm^-1" .* [1, 0, 0]
@info q[1] * 8p.J * p.a / ħ |> u"GHz"

res1 = @showprogress [-χ0(x, q, 1, p)[1] for x in ωrange]
res2 = @showprogress [-χ0(x, q, -1, p)[1] for x in ωrange]

##
fig, ax = plt.subplots()
ax.plot(ωrange .|> ustrip, res1 .|> real .|> ustrip .|> real, label="up")
ax.plot(ωrange .|> ustrip, res2 .|> real .|> ustrip .|> real, label="down")
ax.plot(ωrange .|> ustrip, (res1 .+ res2) .|> real .|> ustrip .|> real, label="sum")
ax.plot(ωrange .|> ustrip, (res1 .+ res2) .|> imag .|> ustrip, "--")
ax.legend()
fig

##
fig, ax = plt.subplots()
ax.plot(ωrange .|> ustrip, 1 ./ (res1 .+ res2) .|> imag .|> ustrip)
fig