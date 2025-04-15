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

includet(scriptsdir("functions.jl"))



function export_params_2D(p=ParamsPlasmon())
    processing_dict = Dict(
        :Ef => x -> generate_latex_command("FermileveltwoD", "\\SI{$(round(p.Ef |> u"eV" |> ustrip;digits=1))}{eV}"),)
    export_params(p, processing_dict)
end

export_2D_ratio(p=ParamsPlasmon()) = println(generate_latex_command("twoDratio", "$(round(((p.mstar^2 - p.m0^2) / p.mstar^2)^(1 / 6);digits=2))"))

export_critical_angle(p=ParamsPlasmon()) = println(
    generate_latex_command("criticalangle", "\\SI{$(round(critical_angle(p); digits=1))}{\\degree}")
)

export_params_2D()
export_2D_ratio()
export_critical_angle()