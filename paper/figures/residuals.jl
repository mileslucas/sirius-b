using ADI
using PyCall
using SiriusB


fits = pyimport("astropy.io.fits")
pro = pyimport("proplot")

# rcParams
pro.rc["image.origin"] = "lower"
pro.rc["image.cmap"] = "inferno"
pro.rc["grid"] = false

epochs = [
    "2020feb04",
    "2020nov21",
    "2020nov28"
]
res = [
    fits.getdata(datadir("epoch_2020feb04", "residuals", "2020feb04_sirius-b_residual_median.fits")),
    fits.getdata(datadir("epoch_2020nov21", "residuals", "2020nov21_sirius-b_residual_median.fits")),
    fits.getdata(datadir("epoch_2020nov28", "residuals", "2020nov28_sirius-b_residual_annular_pca-2.fits")),
]
angles = [
    fits.getdata(datadir("epoch_$epoch", "processed", "$(epoch)_sirius-b_pa.fits"))
    for epoch in epochs
]
fwhms = [
    7.985399074538929,
    7.644946278252633,
    8.216162087416189
]

flat_res = collapse.(res, angles)
sigs = detectionmap.(significance, flat_res, fwhms)
stims = stimmap.(res, angles)

tick_locs = range(25, 174, length=5)
tick_labs = @. string(round(SiriusB.auscale * (tick_locs - 99.5), digits=1))

py"""
import proplot as pro
from matplotlib import patches
fig, axs = pro.subplots(ncols=3, wspace="1em", figwidth="7.5in")

axs[0].imshow($(flat_res[1]))
axs[0].text(6, 6, "median", color="w")
axs[0].format(title="Epoch 2020-02-04")

axs[1].imshow($(flat_res[2]))
axs[1].text(6, 6, "median", color="w")
axs[1].format(title="Epoch 2020-11-21")

axs[2].imshow($(flat_res[3]))
axs[2].text(6, 6, "annular PCA(2)", color="w")
axs[2].format(title="Epoch 2020-11-28")

axs.format(
    xticks=$tick_locs,
    xticklabels=$tick_labs,
    yticks=$tick_locs,
    yticklabels=$tick_labs,
    xlabel="x [AU]",
    ylabel="y [AU]",
)

fig.save($(figuredir("residuals.pdf")))
"""

py"""
import proplot as pro
fig, axs = pro.subplots(ncols=3, nrows=2, wspace="1em", figwidth="7.5in")

axs[0, 0].imshow($(sigs[1]), vmin=0, vmax=5)#, colorbar="r", colorbar_kw=dict(space=0))
axs[0, 0].text(6, 6, "median", color="w")
axs[0, 0].format(title="Epoch 2020-02-04")
axs[0, 1].imshow($(sigs[2]), vmin=0, vmax=5)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[0, 1].text(6, 6, "median", color="w")
axs[0, 1].format(title="Epoch 2020-11-21")
m= axs[0, 2].imshow($(sigs[3]), vmin=0, vmax=5)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[0, 2].text(6, 6, "annular PCA(2)", color="w")
axs[0, 2].format(title="Epoch 2020-11-28")
axs[0, 2].colorbar(m, loc="r", label="significance")


axs[1, 0].imshow($(stims[1]), vmin=0, vmax=1)#, colorbar="r", colorbar_kw=dict(space=0))
axs[1, 0].text(6, 6, "median", color="w")
axs[1, 1].imshow($(stims[2]), vmin=0, vmax=1)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[1, 1].text(6, 6, "median", color="w")
m= axs[1, 2].imshow($(stims[3]), vmin=0, vmax=1)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[1, 2].text(6, 6, "annular PCA(2)", color="w")
axs[1, 2].colorbar(m, loc="r", label="STIM probability")

axs.format(
    xticks=$tick_locs,
    xticklabels=$tick_labs,
    yticks=$tick_locs,
    yticklabels=$tick_labs,
    xlabel="x [AU]",
    ylabel="y [AU]",
)

fig.save($(figuredir("detection_maps.pdf")))
pro.close(fig)

fig, axs = pro.subplots(ncols=3, wspace="1em", figwidth="7.5in")

axs[0].imshow($(sigs[1]), vmin=3)#, colorbar="r", colorbar_kw=dict(space=0))
axs[0].text(6, 6, "median", color="w")
axs[0].format(title="Epoch 2020-02-04")
axs[1].imshow($(sigs[2]), vmin=3)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[1].text(6, 6, "median", color="w")
axs[1].format(title="Epoch 2020-11-21")
m= axs[2].imshow($(sigs[3]), vmin=3)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[2].text(6, 6, "annular PCA(2)", color="w")
axs[2].format(title="Epoch 2020-11-28")
fig.colorbar(m, loc="r", label="significance")

axs.format(
    xticks=$tick_locs,
    xticklabels=$tick_labs,
    yticks=$tick_locs,
    yticklabels=$tick_labs,
    xlabel="x [AU]",
    ylabel="y [AU]",
)

fig.save($(figuredir("sig_threshold.pdf")))
pro.close(fig)
"""

py"""
import proplot as pro
# fig, axs = pro.subplots(ncols=3, wspace="1em", figwidth="7.5in")

# axs[0].imshow($(stims[1]), vmin=0, vmax=1)#, colorbar="r", colorbar_kw=dict(space=0))
# axs[0].text(6, 6, "median", color="w")
# axs[0].format(title="Epoch 2020-02-04")
# axs[1].imshow($(stims[2]), vmin=0, vmax=1)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
# axs[1].text(6, 6, "median", color="w")
# axs[1].format(title="Epoch 2020-11-21")
# m= axs[2].imshow($(stims[3]), vmin=0, vmax=1)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
# axs[2].text(6, 6, "annular PCA(2)", color="w")
# axs[2].format(title="Epoch 2020-11-28")
# fig.colorbar(m, loc="r", label="STIM probability")

# axs.format(
#     xticks=$tick_locs,
#     xticklabels=$tick_labs,
#     yticks=$tick_locs,
#     yticklabels=$tick_labs,
#     xlabel="x [AU]",
#     ylabel="y [AU]",
# )

# fig.save($(figuredir("stim.pdf")))

fig, axs = pro.subplots(ncols=3, wspace="1em", figwidth="7.5in")

axs[0].imshow($(stims[1]), vmin=0.5, vmax=1)#, colorbar="r", colorbar_kw=dict(space=0))
axs[0].text(6, 6, "median", color="w")
axs[0].format(title="Epoch 2020-02-04")
axs[1].imshow($(stims[2]), vmin=0.5, vmax=1)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[1].text(6, 6, "median", color="w")
axs[1].format(title="Epoch 2020-11-21")
m= axs[2].imshow($(stims[3]), vmin=0.5, vmax=1)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[2].text(6, 6, "annular PCA(2)", color="w")
axs[2].format(title="Epoch 2020-11-28")
fig.colorbar(m, loc="r", label="STIM probability")

axs.format(
    xticks=$tick_locs,
    xticklabels=$tick_labs,
    yticks=$tick_locs,
    yticklabels=$tick_labs,
    xlabel="x [AU]",
    ylabel="y [AU]",
)

fig.save($(figuredir("stim_threshold.pdf")))
pro.close(fig)
"""
