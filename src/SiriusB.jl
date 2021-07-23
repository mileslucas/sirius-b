module SiriusB

using Dates
using FITSIO
using SkyCoords

export rootdir, datadir, srcdir, notebookdir, paperdir, figuredir,
       parallactic_angles, gaussian_fit, gaussian_fit_offset

# path configuration
rootdir(args...) = joinpath(@__DIR__(), "..", args...)
datadir(args...) = rootdir("data", args...)
srcdir(args...) = rootdir("src", args...)
notebookdir(args...) = rootdir("notebooks", args...)
paperdir(args...) = rootdir("paper", args...)
figuredir(args...) = paperdir("figures", args...)


"""
    parallactic_angles(filelist; latitude=19.82636)

Calculate the parallactic angless(in degrees) for the given FITS files in the middle of the exposure. Uses the given 
`latitude` which defaults to the Keck II site.
"""
function parallactic_angles(filelist; latitude=19.82636)
    map(filelist) do filename
        # get header
        header = read_header(filename)

        # compute the hour angle at the middle of the frame
        exp_start = Time(header["EXPSTART"])
        exp_stop = Time(header["EXPSTOP"])
        exp_mid = Dates.value((exp_stop - exp_start) / 2)
        # convert exp_mid from nanoseconds of time to degrees
        Δha = exp_mid / 1e9 / 24 / 10
        ha_mid = header["HA"] + Δha

        # precess the star coordinates to the appropriate epoch
        ra = deg2rad(header["RA"])
        dec = deg2rad(header["DEC"])
        coord = FK5Coords{2000}(ra, dec)
        obs_epoch = Date(header["DATE-OBS"])
        yr = year(obs_epoch) + (dayofyear(obs_epoch) - 1) / daysinyear(obs_epoch)
        coord_curr = convert(FK5Coords{yr}, coord)

        # derive true PA at the middle of the frame
        ra = coord_curr.ra
        dec = coord_curr.dec
        PA = -rad2deg(atan(-sind(ha_mid), cos(dec) * tand(latitude) - sin(dec) * cosd(ha_mid)))
        # apply instrumental effects
        PA += header["ROTPOSN"] - header["INSTANGL"] + (header["PARANG"] - header["PARANTEL"])

        return PA
    end
end

include("fitting.jl")
include("badpix.jl")

end