using ADI
using AstroTime
using CSV
using DataFrames
using LinearAlgebra
using Measurements
using PyCall
using PyPlot
using SiriusB

fits = pyimport("astropy.io.fits")
driver = pyimport("orbitize.driver")
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
av_cubes = [AnnulusView(cube; inner=fwhm) for (cube, fwhm) in zip(cubes, fwhms)]
res_cubes = [subtract(GreeDS(2), av_cube; angles=angle) for (av_cube, angle) in zip(av_cubes, angles)]
stims = [stimmap(res_cube, angle) for (res_cube, angle) in zip(res_cubes, angles)]
flat_res = [collapse(res_cube, angle) for (res_cube, angle) in zip(res_cubes, angles)]

parallax = 376.36801 # mas
pxscale = 10 # mas/px


errs = fwhms ./ 2

# guesses from GreeDS outputs of three epochs
pos = [
    110.821 ± errs[1]  95.6414 ± errs[1];
    111.024 ± errs[2]  94.5881 ± errs[2];
    110.98 ± errs[3]  95.9881 ± errs[3];
] # px coordinates

ctr = hcat(center(flat_res[1])...)

delta_px = pos .- ctr
delta = delta_px .* pxscale

###
# plot the epochs
# first, zoom in
crop_res = @. crop(stims, 50)

tick_locs = range(5, 43, length=3)
tick_labs = @. string(round(0.01 * (tick_locs - 24.5), digits=1))


pycenter = hcat(center(crop_res[1])...) .- 1
pos_py = Measurements.value.(delta_px .+ pycenter)


py"""
import proplot as pro
import matplotlib.pyplot as plt
fig, axs = pro.subplots(ncols=3, wspace="1em", figwidth="7.5in")

axs[0].imshow($(crop_res[1]))
axs[0].add_artist(plt.Circle( $(pos_py[1, :]), $(errs[1]), lw=2, color="w", fill=False))
axs[0].format(title="Epoch 2020-02-04")

axs[1].imshow($(crop_res[2]))
axs[1].add_artist(plt.Circle( $(pos_py[2, :]), $(errs[2]), lw=2, color="w", fill=False))
axs[1].format(title="Epoch 2020-11-21")

axs[2].imshow($(crop_res[3]))
axs[2].add_artist(plt.Circle( $(pos_py[3, :]), $(errs[3]), lw=2, color="w", fill=False))
axs[2].format(title="Epoch 2020-11-28")

axs.format(
    xticks=$tick_locs,
    xticklabels=$tick_labs,
    yticks=$tick_locs,
    yticklabels=$tick_labs,
    xlabel="x [arcsec]",
    ylabel="y [arcsec]",
)

fig.save($(figuredir("orbit_frames.pdf")))
"""



###

# translate into offset and angle

dist = mapslices(norm, delta, dims=2)
angs = mapslices(v -> atand(reverse(v)...) - 90, delta, dims=2)

epochs_utc = [
    TTEpoch(2020, 2, 4),
    TTEpoch(2020, 11, 21),
    TTEpoch(2020, 11, 28)
]

epochs_mjd = @. AstroTime.value(modified_julian(epochs_utc))

data = DataFrame(
    epoch = epochs_mjd,
    object = fill(1, length(epochs)),
    sep = Measurements.value.(dist) |> vec,
    sep_err = Measurements.uncertainty.(dist) |> vec,
    pa = Measurements.value.(angs) |> vec,
    pa_err = Measurements.uncertainty.(angs) |> vec
)
@info data

fname = datadir("sirius-b_astrometry_data.csv")
data |> CSV.write(fname)

# now start setting up orbitize
mydriver = driver.Driver(
    fname,
    "OFTI",
    1, # number of bodies
    3.081, # total system mass M_sun https://iopscience.iop.org/article/10.3847/1538-4357/aa6af8/pdf
    parallax, # parallax in mas https://iopscience.iop.org/article/10.3847/1538-4357/aa6af8/pdf
    mass_err = 0.034,
    plx_err = 1.4
)

s = mydriver.sampler
orbits = s.run_sampler(Int(1e6))

# orbit plot
s.results.plot_orbits(start_mjd=epochs_mjd[1], num_orbits_to_plot=1000)
savefig(figuredir("orbit.pdf"))

# corner plot
s.results.plot_corner(param_list=["sma1", "ecc1", "inc1", "aop1", "pan1","tau1"])
savefig(figuredir("orbit_corner.pdf"))
