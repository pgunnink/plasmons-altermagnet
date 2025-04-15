using DrWatson
@quickactivate
using Printf
function generate_latex_command(name, cmd)
    "\\newcommand{\\$(name)}{$cmd}"
end

function export_params(p, processing_dict)
    for x in keys(processing_dict)
        println(processing_dict[x](getfield(p, x)))
    end
end
function convert_e_to_num_latex(x)
    exp_string = @sprintf "%.0e" x
    return exp_string[2:end]
end

includet(scriptsdir("functions_3D.jl"))
include(scriptsdir("functions_analysis.jl"))

function export_params_3D(p=ParamsPlasmon())
    processing_dict = Dict(
        :m0 => x -> generate_latex_command("mzero", "$(round(p.m0 / me ;digits=2))"),
        :mstar => x -> generate_latex_command("mstar", "$(round(p.mstar / p.m0 ;digits=2))"),
        :Ef => x -> generate_latex_command("Fermilevel", "\\SI{$(p.Ef |> u"eV" |> ustrip)}{eV}"),)
    export_params(p, processing_dict)
end
export_splitting(p=ParamsPlasmon()) = println(generate_latex_command("splitting", "\\SI{$(round(splitting(p) |> u"eV" |> ustrip;digits=1))}{eV}"))

export_magnetic_moment(p=ParamsPlasmon(), q=0.05) = println(
    generate_latex_command("muB", round(determine_μ([q * kF(p), 0kF(p), 0kF(p)], p); digits=3))
)


export_magnetic_moment_shift(p=ParamsPlasmon(), q=0.05, B=1u"T") = println(
    generate_latex_command("muBShift",
        "\\SI{$(
        round(determine_μ([q * kF(p), 0kF(p), 0kF(p)], p) * Unitful.μB *  B |> u"μeV" |> ustrip ; digits=1)
        )}{\\mu eV}")
)

export_energy_demon(p=ParamsPlasmon(), q=0.05) = println(
    generate_latex_command("demonenergy",
        "\\SI{$(
        round(plasmon_pole([q * kF(p), 0kF(p), 0kF(p)], p) * vF(p) * q * kF(p) * ħ  |> u"meV" |> ustrip ; digits=1)
        )}{meV}")
)

export_magnetic_moment_shift_relative(p=ParamsPlasmon(), q=0.05, B=1u"T") = println(
    generate_latex_command("muBShiftrelative",
        round(determine_μ([q * kF(p), 0kF(p), 0kF(p)], p) * Unitful.μB * B / (plasmon_pole([q * kF(p), 0kF(p), 0kF(p)], p) * vF(p) * q * kF(p) * ħ) * 100 |> upreferred; digits=1)
    ))


function export_max_magnetic_moment(p=ParamsPlasmon(), q=0.05)
    θrange = range(0, stop=2pi / 4, length=101) .- 2pi / 8
    mus = [abs(determine_μ(q * kF(p) * [cos(θ), sin(θ), 0], p)) for θ in θrange]
    maxmu = maximum(filter(!isnan, mus))
    println(
        generate_latex_command("muBmax",
            round(maxmu; digits=1)))

end

export_params_3D()
export_splitting()
export_magnetic_moment()
export_magnetic_moment_shift()
export_energy_demon()
export_magnetic_moment_shift_relative()
export_max_magnetic_moment()