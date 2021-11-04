using CSV
using DataFrames
using Interpolations
using PyCall
using SiriusB
using Statistics

pro = pyimport("proplot")
# rcParams
pro.use_style("ggplot")
pro.rc["image.origin"] = "lower"
pro.rc["image.cmap"] = "inferno"
pro.rc["grid"] = false


contrast_curves = [
    DataFrame(CSV.File(datadir("epoch_2020feb04", "residuals", "2020feb04_contrast-curve_median.csv"))),
    DataFrame(CSV.File(datadir("epoch_2020nov21", "residuals", "2020nov21_contrast-curve_median.csv"))),
    DataFrame(CSV.File(datadir("epoch_2020nov28", "residuals", "2020nov28_contrast-curve_annular-pca-2.csv"))),
]

# convert distances from pixels to AU
distances = [curve.distance .* SiriusB.auscale for curve in contrast_curves]

function contrast_to_absmag(y)
    Δmag = contrast_to_dmag(y)
    appmag = SiriusB.appmag + Δmag
    return appmag - distance_modulus(SiriusB.distance)
end

## plot each epoch as contrast/absolute mag
py"""
import proplot as pro
fig, axs = pro.subplots(width="7.5in", height="4in")

# curves
axs.plot($(distances[1]), $(contrast_curves[1].contrast), c="C0", label="2020-02-04")
axs.plot($(distances[1]), $(contrast_curves[1].contrast_corr), c="C0", ls="--")

axs.plot($(distances[2]), $(contrast_curves[2].contrast), c="C1", label="2020-11-21")
axs.plot($(distances[2]), $(contrast_curves[2].contrast_corr), c="C1", ls="--")

axs.plot($(distances[3]), $(contrast_curves[3].contrast), c="C5", label="2020-11-28")
axs.plot($(distances[3]), $(contrast_curves[3].contrast_corr), c="C5", ls="--")

# formatting
axs.dualx(lambda x: x * $(SiriusB.parallax), label="separation [arcsec]")
axs.format(
    yscale="log",
    xlabel="projected separation [AU]",
    ylabel="contrast ",
    grid=True,
    yformatter="log",
    xlim=(None, 2),
)

# he's mad, he's gona and made a custom legend!
from matplotlib.lines import Line2D
leg1 = axs.legend(ncol=1, loc="ur")
elements = [
    Line2D([0], [0], color="k", alpha=0.6, label="Gaussian"),
    Line2D([0], [0], color="k", alpha=0.6, ls="--", label="Student-t"),
]
axs.legend(handles=elements, queue=True, loc="uc")
axs.add_artist(leg1)

# stability limit
ylims=axs.get_ylim()
axs.vlines(1.5, *ylims, color="k", alpha=0.4, ls="-.")
mid = (ylims[0] + ylims[1]) / 40
# print(mid)
axs.text(1.45, mid, "dynamical stability limit", color="k", alpha=0.5, va="center", ha="left", rotation="vertical")
axs.set_ylim(ylims)

ax2 = axs.alty(reverse=True, label="abs. mag [$\mathrm{L^\prime}$]")
ax2.plot($(distances[2]), $(contrast_to_absmag.(contrast_curves[2].contrast)), alpha=0)

fig.savefig($(figuredir("contrast_curves.pdf")))
pro.close(fig)
"""

for curve in contrast_curves
    curve.absmags = contrast_to_absmag.(curve.contrast)
    curve.masses_225 = map(curve.absmags) do absmag
        table_interpolate(ATMO2020, :Age => 0.225, :MKO_Lp => absmag, :Mass_MJ)
    end
    curve.masses_126 = map(curve.absmags) do absmag
        table_interpolate(ATMO2020, :Age => 0.126, :MKO_Lp => absmag, :Mass_MJ)
    end
    curve.absmags_corr = contrast_to_absmag.(curve.contrast_corr)
    curve.masses_corr_225 = map(curve.absmags_corr) do absmag
        table_interpolate(ATMO2020, :Age => 0.225, :MKO_Lp => absmag, :Mass_MJ)
    end
    curve.masses_corr_126 = map(curve.absmags_corr) do absmag
        table_interpolate(ATMO2020, :Age => 0.126, :MKO_Lp => absmag, :Mass_MJ)
    end
    # these try catch curves are because Sonora grid bottoms out
    curve.masses_225_SonoraSolar = map(curve.absmags) do absmag
        table_interpolate(SonoraSolar, :age_Gyr => 0.225, :Keck_Lp => absmag, :M_Mjup)
    end
    curve.masses_126_SonoraSolar = map(curve.absmags) do absmag
        table_interpolate(SonoraSolar, :age_Gyr => 0.126, :Keck_Lp => absmag, :M_Mjup)
    end
    curve.masses_corr_225_SonoraSolar = map(curve.absmags_corr) do absmag
        table_interpolate(SonoraSolar, :age_Gyr => 0.225, :Keck_Lp => absmag, :M_Mjup)
    end
    curve.masses_corr_126_SonoraSolar = map(curve.absmags_corr) do absmag
        table_interpolate(SonoraSolar, :age_Gyr => 0.126, :Keck_Lp => absmag, :M_Mjup)
    end
    # metalrich
    curve.masses_225_SonoraMetalRich = map(curve.absmags) do absmag
        table_interpolate(SonoraMetalRich, :age_Gyr => 0.225, :Keck_Lp => absmag, :M_Mjup)
    end
    curve.masses_126_SonoraMetalRich = map(curve.absmags) do absmag
        table_interpolate(SonoraMetalRich, :age_Gyr => 0.126, :Keck_Lp => absmag, :M_Mjup)
    end
    curve.masses_corr_225_SonoraMetalRich = map(curve.absmags_corr) do absmag
        table_interpolate(SonoraMetalRich, :age_Gyr => 0.225, :Keck_Lp => absmag, :M_Mjup)
    end
    curve.masses_corr_126_SonoraMetalRich = map(curve.absmags_corr) do absmag
        table_interpolate(SonoraMetalRich, :age_Gyr => 0.126, :Keck_Lp => absmag, :M_Mjup)
    end
end


## plot best epoch as mass limits
py"""
import proplot as pro
fig, axs = pro.subplots(ncols=2, share=2, width="7.5in", height="4in")

# atmo2020 curves
axs[0].plot($(distances[2]), $(contrast_curves[2].masses_225), c="C2", label="225 Myr")
axs[0].plot($(distances[2]), $(contrast_curves[2].masses_corr_225), c="C2", ls="--")
axs[0].fill_between($(distances[2]), $(contrast_curves[2].masses_225), $(contrast_curves[2].masses_corr_225), color="C2", lw=0, alpha=0.1)

axs[0].plot($(distances[2]), $(contrast_curves[2].masses_126), c="C3", label="126 Myr")
axs[0].plot($(distances[2]), $(contrast_curves[2].masses_corr_126), c="C3", ls="--")
axs[0].fill_between($(distances[2]), $(contrast_curves[2].masses_126), $(contrast_curves[2].masses_corr_126), color="C3", lw=0, alpha=0.1)


# curves
axs[1].plot($(distances[2]), $(contrast_curves[2].masses_225_SonoraSolar), c="C4", label="225 Myr")
axs[1].plot($(distances[2]), $(contrast_curves[2].masses_126_SonoraSolar), c="C5", label="126 Myr")
axs[1].plot($(distances[2]), $(contrast_curves[2].masses_225_SonoraMetalRich), c="C4", ls=":")
axs[1].plot($(distances[2]), $(contrast_curves[2].masses_126_SonoraMetalRich), c="C5", ls=":")


# # model mass limit
[ax.hlines([1, 3, 5, 6, 7, 9, 10], *ax.get_xlim(), color="w", alpha=0.5, lw=0.75) for ax in axs]

# formatting
axs.format(
    abc=True,
    yscale=pro.LogScale(base=2),
    xlabel="projected separation [AU]",
    ylabel="companion mass [$M_J$] ",
    grid=True,
    yformatter="auto",
    xlim=(None, 1.5),
)
axs[0].format(title="ATMO2020")
axs[1].format(title="Sonora Bobcat")

# custom legend
axs.legend(ncol=1, loc="ur")
elements = [
    Line2D([0], [0], color="k", alpha=0.6, label="Gaussian"),
    Line2D([0], [0], color="k", alpha=0.6, ls="--", label="Student-t"),
]
axs[0].legend(handles=elements, queue=True, loc="t")
elements = [
    Line2D([0], [0], color="k", alpha=0.6, label="[M/H] = 0.0"),
    Line2D([0], [0], color="k", alpha=0.6, ls=":", label="[M/H] = +0.5"),
]
axs[1].legend(handles=elements, queue=True, loc="t")

fig.savefig($(figuredir("mass_curves.pdf")))
pro.close(fig)
"""


