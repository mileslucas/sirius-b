using PlotUtils: zscale
using PyCall
using SiriusB

fits = pyimport("astropy.io.fits")
pro = pyimport("proplot")

# rcParams
pro.rc["image.origin"] = "lower"
pro.rc["image.cmap"] = "inferno"

data_cube = fits.getdata(datadir("epoch_2020feb04", "processed", "2020feb04_sirius-b_cube_calib.fits"))
frame = data_cube[88, :, :]

fig, axs = pro.subplots(figwidth="3.5in")
ax = first(axs)
vmin, vmax = zscale(frame)
ax.imshow(frame; vmin, vmax)
ax.format(xlabel="x", ylabel="y", grid=false)
fig.save(joinpath(@__DIR__, "spike.pdf"))
