using HCIToolbox
using PlotUtils: zscale
using PyCall
using SiriusB


fits = pyimport("astropy.io.fits")
pro = pyimport("proplot")

# rcParams
pro.rc["image.origin"] = "lower"
pro.rc["image.cmap"] = "inferno"
pro.rc["grid"] = false

cube= fits.getdata(datadir("epoch_2020nov21", "processed", "2020nov21_sirius-b_cube_calib_registered_crop.fits"))

flat = collapse(cube)

function norm_psf(frame)
    tmp = frame .- minimum(frame)
    return tmp ./ maximum(tmp)
    # return log10.(tmp ./ maximum(tmp))
end
psf = norm_psf(flat)

tick_locs = range(20, 179, length=5)
tick_labs = @. string(round(0.01 * (tick_locs - 99.5), digits=1))

py"""
import proplot as pro
from matplotlib.colors import LogNorm

fig, axs = pro.subplots(figwidth="3.5in")

axs[0].imshow($psf, norm=LogNorm(vmin=3e-3, vmax=5e-2))
axs[0].format(title="Epoch 2020-11-21")
axs.format(
    xticks=$tick_locs,
    xticklabels=$tick_labs,
    yticks=$tick_locs,
    yticklabels=$tick_labs,
    xlabel="x [arcsec]",
    ylabel="y [arcsec]",
)
fig.save($(joinpath(@__DIR__, "psf.pdf")))
"""