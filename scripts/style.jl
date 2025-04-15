using PyPlot
using DrWatson
plt.style.use(joinpath(@__DIR__, "prl.mplstyle"))

save_dir = joinpath(papersdir(), "figures/")
column_width = 3 + 3 / 8
fig_width = column_width
fig_height = fig_width / 1.618


function set_rticks_custom(ax)
    ax.set_rlabel_position(0)
    ax.grid(color="black", alpha=0.2)

    ax.set_rticks([0.05, 0.1], [L"0.05", L"0.1"])
    ax.tick_params(axis="y", labelsize=7)
    ax.tick_params(axis="y", which="minor", bottom=false)
    ax.text(-0.2, 0.085, L"\hbar\omega/\epsilon_F", size=7)
    for label in ax.get_yticklabels()
        label.set_horizontalalignment("center")
    end
end