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
    tmp = frame .- minimum(frame) .+ 1
    return tmp ./ maximum(tmp)
end
psf = norm_psf(flat)

tick_locs = range(25, 174, length=5)
tick_labs = @. string(round(SiriusB.auscale * (tick_locs - 99.5), digits=1))

py"""
import proplot as pro
from matplotlib.colors import LogNorm

fig, axs = pro.subplots(figwidth="4in")

axs.imshow($psf, 
    norm=LogNorm(vmin=3e-3, vmax=5e-2),
    colorbar="r",
)
# axs.format(title="Epoch 2020-11-21")
axs.format(
    xticks=$tick_locs,
    xticklabels=$tick_labs,
    yticks=$tick_locs,
    yticklabels=$tick_labs,
    xlabel="x [AU]",
    ylabel="y [AU]",
)
fig.save($(figuredir("psf.pdf")))
"""