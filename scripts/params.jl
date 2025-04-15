using Parameters
using Unitful

@with_kw mutable struct ParamsPlasmon
    m0 = 0.4me
    mstar = m0 * 1.25
    # mDOS = m0 * (mstar^2 / (mstar^2 - m0^2))^(1 / 3)
    Ef = 0.5u"eV"
    ϵ = 1
end