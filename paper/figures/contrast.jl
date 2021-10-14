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
    yformatter="sci",
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

ax2 = axs.alty(reverse=True, label="abs. mag [$M^{Lp}$]")
ax2.plot($(distances[2]), $(contrast_to_absmag.(contrast_curves[2].contrast)), alpha=0)

fig.savefig($(figuredir("contrast_curves.pdf")))
pro.close(fig)
"""

for curve in contrast_curves
    curve.absmags = contrast_to_absmag.(curve.contrast)
    curve.masses_226 = map(curve.absmags) do absmag
        table_interpolate(ATMO2020, :Age => 0.226, :MKO_Lp => absmag, :Mass_MJ)
    end
    curve.masses_125 = map(curve.absmags) do absmag
        table_interpolate(ATMO2020, :Age => 0.125, :MKO_Lp => absmag, :Mass_MJ)
    end
    curve.absmags_corr = contrast_to_absmag.(curve.contrast_corr)
    curve.masses_corr_226 = map(curve.absmags_corr) do absmag
        table_interpolate(ATMO2020, :Age => 0.226, :MKO_Lp => absmag, :Mass_MJ)
    end
    curve.masses_corr_125 = map(curve.absmags_corr) do absmag
        table_interpolate(ATMO2020, :Age => 0.125, :MKO_Lp => absmag, :Mass_MJ)
    end
end

## plot best epoch as mass limits
    py"""
import proplot as pro
fig, axs = pro.subplots(width="7.5in", height="4in")

# curves
axs.plot($(distances[2]), $(contrast_curves[2].masses_226), c="C0", label="226 Myr")
axs.plot($(distances[2]), $(contrast_curves[2].masses_corr_226), c="C0", ls="--")
axs.fill_between($(distances[2]), $(contrast_curves[2].masses_226), $(contrast_curves[2].masses_corr_226), color="C0", lw=0, alpha=0.1)

axs.plot($(distances[2]), $(contrast_curves[2].masses_125), c="C1", label="125 Myr")
axs.plot($(distances[2]), $(contrast_curves[2].masses_corr_125), c="C1", ls="--")
axs.fill_between($(distances[2]), $(contrast_curves[2].masses_125), $(contrast_curves[2].masses_corr_125), color="C1", lw=0, alpha=0.1)

# model mass limit
xlims = axs.get_xlim()
axs.hlines($(minimum(ATMO2020.Mass_MJ)), *xlims, color="k", alpha=0.3, linestyle=":", label="model mass limit")
axs.hlines([1, 3, 5, 6, 7, 9, 10], *xlims, color="w", alpha=0.5, lw=0.75)

# formatting
axs.dualx(lambda x: x * $(SiriusB.parallax), label="separation [arcsec]")
axs.format(
    yscale=pro.LogScale(base=2),
    xlabel="projected separation [AU]",
    ylabel="companion mass [$M_J$] ",
    grid=True,
    yformatter="auto",
    xlim=(None, 2),
)

# custom legend
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
mid = $log2(ylims[0] + ylims[1]) - 0.5
# print(mid)
axs.text(1.45, mid, "dynamical stability limit", color="k", alpha=0.5, va="center", ha="left", rotation="vertical")
axs.set_ylim(ylims)


fig.savefig($(figuredir("mass_curves.pdf")))
pro.close(fig)
"""
