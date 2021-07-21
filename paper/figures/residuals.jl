using ADI
using PyCall
using SiriusB


fits = pyimport("astropy.io.fits")
pro = pyimport("proplot")

# rcParams
pro.rc["image.origin"] = "lower"
pro.rc["image.cmap"] = "inferno"
pro.rc["grid"] = false

cube1 = fits.getdata(datadir("epoch_2020feb04", "processed", "2020feb04_sirius-b_cube_calib_registered_crop.fits"))
angles1 = fits.getdata(datadir("epoch_2020feb04", "processed", "2020feb04_sirius-b_pa.fits"))

cube2 = fits.getdata(datadir("epoch_2020nov21", "processed", "2020nov21_sirius-b_cube_calib_registered_crop.fits"))
angles2 = fits.getdata(datadir("epoch_2020nov21", "processed", "2020nov21_sirius-b_pa.fits"))

cube3 = fits.getdata(datadir("epoch_2020nov28", "processed", "2020nov28_sirius-b_cube_calib_registered_crop.fits"))
angles3 = fits.getdata(datadir("epoch_2020nov28", "processed", "2020nov28_sirius-b_pa.fits"))


res1 = subtract(Classic(), mask_circle(cube1, 7.985399074538929))
res2 = subtract(Classic(), mask_circle(cube2, 7.644946278252633))
res3 = subtract(PCA(2), mask_circle(cube3, 8.216162087416189))

flat_res1 = collapse(res1, angles1)
flat_res2 = collapse(res2, angles2)
flat_res3 = collapse(res3, angles3)

sn_1 = detectionmap(flat_res1, 7.985399074538929)
sn_2 = detectionmap(flat_res2, 7.644946278252633)
sn_3 = detectionmap(flat_res3, 8.216162087416189)

stim_1 = stimmap(res1, angles1)
stim_2 = stimmap(res2, angles2)
stim_3 = stimmap(res3, angles3)


tick_locs = range(20, 179, length=5)
tick_labs = @. string(round(0.01 * (tick_locs - 99.5), digits=1))

py"""
import proplot as pro
fig, axs = pro.subplots(ncols=3, wspace="1em", figwidth="7in")

axs[0].imshow($flat_res1)#, colorbar="r", colorbar_kw=dict(space=0))
axs[0].format(title="Epoch 2020-02-04")
axs[1].imshow($flat_res2)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[1].format(title="Epoch 2020-11-21")
axs[2].imshow($flat_res3)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[2].format(title="Epoch 2020-11-28")

axs.format(
    xticks=$tick_locs,
    xticklabels=$tick_labs,
    yticks=$tick_locs,
    yticklabels=$tick_labs,
    xlabel="x [arcsec]",
    ylabel="y [arcsec]",
)

fig.save($(joinpath(@__DIR__, "residuals.pdf")))
"""

py"""
import proplot as pro
fig, axs = pro.subplots(ncols=3, wspace="1em", figwidth="7in")

axs[0].imshow($sn_1, vmin=0, vmax=5)#, colorbar="r", colorbar_kw=dict(space=0))
axs[0].format(title="Epoch 2020-02-04")
axs[1].imshow($sn_2, vmin=0, vmax=5)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[1].format(title="Epoch 2020-11-21")
m= axs[2].imshow($sn_3, vmin=0, vmax=5)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[2].format(title="Epoch 2020-11-28")
fig.colorbar(m, loc="r", label="S/N")

axs.format(
    xticks=$tick_locs,
    xticklabels=$tick_labs,
    yticks=$tick_locs,
    yticklabels=$tick_labs,
    xlabel="x [arcsec]",
    ylabel="y [arcsec]",
)

fig.save($(joinpath(@__DIR__, "snr.pdf")))
"""

py"""
import proplot as pro
fig, axs = pro.subplots(ncols=3, wspace="1em", figwidth="7in")

axs[0].imshow($stim_1, vmin=0, vmax=1)#, colorbar="r", colorbar_kw=dict(space=0))
axs[0].format(title="Epoch 2020-02-04")
axs[1].imshow($stim_2, vmin=0, vmax=1)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[1].format(title="Epoch 2020-11-21")
m= axs[2].imshow($stim_3, vmin=0, vmax=1)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[2].format(title="Epoch 2020-11-28")
fig.colorbar(m, loc="r", label="STIM probability")

axs.format(
    xticks=$tick_locs,
    xticklabels=$tick_labs,
    yticks=$tick_locs,
    yticklabels=$tick_labs,
    xlabel="x [arcsec]",
    ylabel="y [arcsec]",
)

fig.save($(joinpath(@__DIR__, "stim.pdf")))
"""

###

# py"""
# import proplot as pro
# fig, axs = pro.subplots(ncols=3, wspace="3.5em", figwidth="7in")

# axs[0].imshow($flat_res1, colorbar="r", colorbar_kw=dict(space=0))
# axs[0].format(title="Residual Frame")
# axs[1].imshow($sn_1, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
# axs[1].format(title="Gaussian S/N map")
# axs[2].imshow($stim_1, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
# axs[2].format(title="STIM map")

# axs.format(
#     abc=True,
#     xticks=$tick_locs,
#     xticklabels=$tick_labs,
#     yticks=$tick_locs,
#     yticklabels=$tick_labs,
#     xlabel="x [arcsec]",
#     ylabel="y [arcsec]",
#     suptitle="Epoch 2020-02-04"
# )

# fig.save($(joinpath(@__DIR__, "residuals_2020feb04.pdf")))
# """

# py"""
# import proplot as pro
# fig, axs = pro.subplots(ncols=3, wspace="3.5em", figwidth="7in")

# axs[0].imshow($flat_res2, colorbar="r", colorbar_kw=dict(space=0))
# axs[0].format(title="Residual Frame")
# axs[1].imshow($sn_2, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
# axs[1].format(title="Gaussian S/N map")
# axs[2].imshow($stim_2, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
# axs[2].format(title="STIM map")

# axs.format(
#     abc=True,
#     xticks=$tick_locs,
#     xticklabels=$tick_labs,
#     yticks=$tick_locs,
#     yticklabels=$tick_labs,
#     xlabel="x [arcsec]",
#     ylabel="y [arcsec]",
#     suptitle="Epoch 2020-11-21"
# )

# fig.save($(joinpath(@__DIR__, "residuals_2020nov21.pdf")))
# """

# py"""
# import proplot as pro
# fig, axs = pro.subplots(ncols=3, wspace="3.5em", figwidth="7in")

# axs[0].imshow($flat_res3, colorbar="r", colorbar_kw=dict(space=0))
# axs[0].format(title="Residual Frame")
# axs[1].imshow($sn_3, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
# axs[1].format(title="Gaussian S/N map")
# axs[2].imshow($stim_3, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
# axs[2].format(title="STIM map")

# axs.format(
#     abc=True,
#     xticks=$tick_locs,
#     xticklabels=$tick_labs,
#     yticks=$tick_locs,
#     yticklabels=$tick_labs,
#     xlabel="x [arcsec]",
#     ylabel="y [arcsec]",
#     suptitle="Epoch 2020-11-28"
# )

# fig.save($(joinpath(@__DIR__, "residuals_2020nov28.pdf")))
# """
