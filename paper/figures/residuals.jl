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
cubes = [
    fits.getdata(datadir("epoch_$epoch", "processed", "$(epoch)_sirius-b_cube_calib_registered_crop.fits"))
    for epoch in epochs
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

res = [
    subtract(Classic(), AnnulusView(cubes[1]; inner=fwhms[1])),
    subtract(Classic(), AnnulusView(cubes[2]; inner=fwhms[2])),
    subtract(Framewise(PCA(2), delta_rot=0.5), MultiAnnulusView(cubes[3], fwhms[3]; inner=fwhms[3]); angles=angles[3])
]
flat_res = collapse.(res, angles)
sigs = detectionmap.(significance, flat_res, fwhms)
stims = stimmap.(res, angles)

tick_locs = range(20, 179, length=5)
tick_labs = @. string(round(0.01 * (tick_locs - 99.5), digits=1))

py"""
import proplot as pro
fig, axs = pro.subplots(ncols=3, wspace="1em", figwidth="7.5in")

axs[0].imshow($(flat_res[1]))#, colorbar="r", colorbar_kw=dict(space=0))
axs[0].text(6, 6, "median", color="w")
axs[0].format(title="Epoch 2020-02-04")
axs[1].imshow($(flat_res[2]))#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[1].text(6, 6, "median", color="w")
axs[1].format(title="Epoch 2020-11-21")
axs[2].imshow($(flat_res[3]))#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[2].text(6, 6, "annular PCA(2)", color="w")
axs[2].format(title="Epoch 2020-11-28")

axs.format(
    xticks=$tick_locs,
    xticklabels=$tick_labs,
    yticks=$tick_locs,
    yticklabels=$tick_labs,
    xlabel="x [arcsec]",
    ylabel="y [arcsec]",
)

fig.save($(figuredir("residuals.pdf")))
"""

py"""
import proplot as pro
fig, axs = pro.subplots(ncols=3, wspace="1em", figwidth="7.5in")

axs[0].imshow($(sigs[1]), vmin=0, vmax=5)#, colorbar="r", colorbar_kw=dict(space=0))
axs[0].text(6, 6, "median", color="w")
axs[0].format(title="Epoch 2020-02-04")
axs[1].imshow($(sigs[2]), vmin=0, vmax=5)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[1].text(6, 6, "median", color="w")
axs[1].format(title="Epoch 2020-11-21")
m= axs[2].imshow($(sigs[3]), vmin=0, vmax=5)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[2].text(6, 6, "annular PCA(2)", color="w")
axs[2].format(title="Epoch 2020-11-28")
fig.colorbar(m, loc="r", label="significance")

axs.format(
    xticks=$tick_locs,
    xticklabels=$tick_labs,
    yticks=$tick_locs,
    yticklabels=$tick_labs,
    xlabel="x [arcsec]",
    ylabel="y [arcsec]",
)

fig.save($(figuredir("sig.pdf")))
"""

py"""
import proplot as pro
fig, axs = pro.subplots(ncols=3, wspace="1em", figwidth="7.5in")

axs[0].imshow($(stims[1]), vmin=0, vmax=1)#, colorbar="r", colorbar_kw=dict(space=0))
axs[0].text(6, 6, "median", color="w")
axs[0].format(title="Epoch 2020-02-04")
axs[1].imshow($(stims[2]), vmin=0, vmax=1)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[1].text(6, 6, "median", color="w")
axs[1].format(title="Epoch 2020-11-21")
m= axs[2].imshow($(stims[3]), vmin=0, vmax=1)#, colorbar="r", vmin=0, colorbar_kw=dict(space=0))
axs[2].text(6, 6, "annular PCA(2)", color="w")
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

fig.save($(figuredir("stim.pdf")))
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
