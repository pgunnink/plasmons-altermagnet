using DrWatson
@quickactivate
using Unitful
using PyPlot


includet(scriptsdir("functions_3D.jl"))
includet(scriptsdir("style.jl"))

p = ParamsPlasmon()

fig, ax = plt.subplots()
Nk = 100
krange = range(0u"cm^-1", stop=pi / 5u"Å" / 5, length=Nk)

ax.plot(krange .|> ustrip, [ϵ([k, 0k, 0k], 1, p) for k in krange] .|> u"eV" .|> ustrip)
ax.plot(krange .|> ustrip, [ϵ([k, 0k, 0k], -1, p) for k in krange] .|> u"eV" .|> ustrip)


@info splitting(p)
fig