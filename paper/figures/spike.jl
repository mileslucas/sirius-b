using PlotUtils: zscale
using PyCall
using SiriusB

fits = pyimport("astropy.io.fits")
pro = pyimport("proplot")

# rcParams
pro.rc["image.origin"] = "lower"
pro.rc["image.cmap"] = "inferno"
pro.rc["grid"] = false

data_cube = fits.getdata(datadir("epoch_2020nov21", "processed", "2020nov21_sirius-b_cube_calib.fits"))
frame = Float32.(data_cube[169, :, :])

fig, axs = pro.subplots(figwidth="4in")
ax = first(axs)
vmin, vmax = zscale(frame)
ax.imshow(frame; vmin, vmax, colorbar="r")
ax.format(xlabel="x [px]", ylabel="y [px]")
fig.save(figuredir("spike.pdf"))
