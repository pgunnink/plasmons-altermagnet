using Unitful
import Unitful.ħ, Unitful.me, Unitful.ϵ0

bandwidth = 0.5u"eV"
a = 5u"Å"

path_along_S = sqrt(2 * (pi / a)^2)
kF_measured = 0.2path_along_S

meff = ħ^2 * kF_measured^2 / 2bandwidth |> upreferred
α = 0.4ħ^2 / meff
@info ħ^2 * kF_measured^2 / 2meff + α * kF_measured^2 |> u"eV"

@info meff / me
##
