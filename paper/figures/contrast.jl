using CSV
using DataFrames
using PyCall
using SiriusB

pro = pyimport("proplot")
# rcParams
pro.use_style("ggplot")
pro.rc["image.origin"] = "lower"
pro.rc["image.cmap"] = "inferno"
pro.rc["grid"] = false

py"""
import sys
sys.path.append($(srcdir()))
"""
models = pyimport("mass_models")
model = models.MassModel(srcdir("mass_models", "data", "model.AMES-Cond-2000.M-0.0.MKO.Vega.txt"))

parallax = 376.6801e-3 # arcseconds
pxscale = 0.01 # arcseconds/px
auscale = pxscale / parallax # AU/px
age = 226; # Myr

contrast_curves = [
    DataFrame(CSV.File(datadir("epoch_2020feb04", "residuals", "2020feb04_contrast-curve_median.csv"))),
    DataFrame(CSV.File(datadir("epoch_2020nov21", "residuals", "2020nov21_contrast-curve_median.csv"))),
    DataFrame(CSV.File(datadir("epoch_2020nov28", "residuals", "2020nov28_contrast-curve_annular-pca-2.csv"))),
]

contrast_to_dmag(contrast) = -2.5 * log10(contrast)

function dmag_to_mass(dmag)
    model.contrast_to_mass(
        deltaMag=dmag,
        age_Myr=age,
        band="Lp",
        dist_pc=1/parallax,
        stellar_apparent_mag=9.1
    )
end

distances = [curve.distance .* auscale for curve in contrast_curves]

masses = map(contrast_curves) do curve
    dmag = contrast_to_dmag.(curve.contrast)
    dmag_to_mass.(dmag)
end

masses_corr = map(contrast_curves) do curve
    dmag = contrast_to_dmag.(curve.contrast_corr)
    dmag_to_mass.(dmag)
end

py"""
import proplot as pro
fig, axs = pro.subplots(width="7.5in", height="4in")

axs.plot($(distances[1]), $(masses[1]), c="C0", label="  2020-02-04\n     median")
axs.plot($(distances[1]), $(masses_corr[1]), c="C0", ls="--")

axs.plot($(distances[2]), $(masses[2]), c="C1", label="  2020-11-21\n     median")
axs.plot($(distances[2]), $(masses_corr[2]), c="C1", ls="--")

axs.plot($(distances[3]), $(masses[3]), c="C5", label="  2020-11-28\nannular PCA(2)")
axs.plot($(distances[3]), $(masses_corr[3]), c="C5", ls="--")

ylims=axs.get_ylim()
axs.vlines(1.5, *ylims, color="k", alpha=0.5, ls="-.")
mid = (ylims[0] + ylims[1]) / 2
axs.text(1.45, mid, "dynamic stability limit", color="k", alpha=0.5, va="center", ha="left", rotation="vertical")

xlims = axs.get_xlim()
axs.hlines(0.0005 * 1047.9, *xlims, color="k", alpha=0.3, linestyle=":", label="model mass limit")

axs.dualx(lambda x: x * $parallax, label="separation [arcsec]")
axs.legend(ncol=1)
axs.format(
    # yscale="log",
    xlabel="projected separation [AU]",
    ylabel="companion mass [M$_J$]",
    grid=True,
    ylim=(0, ylims[1]),
    xlim=(xlims[0], 2),
    yticks="auto"
)

ax2 = axs.alty(reverse=True, label="Δ mag")
ax2.plot($(distances[2]), $(contrast_to_dmag.(contrast_curves[2].contrast)), alpha=0)

fig.savefig($(figuredir("contrast_curves.pdf")))
"""
