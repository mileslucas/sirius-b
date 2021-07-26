

"""
    parallactic_angles(filelist; latitude=19.82636)

Calculate the parallactic angles (in degrees) for the given FITS files in the middle of the exposure. Uses the given 
`latitude` which defaults to the Keck II site.
"""
parallactic_angles(filelist; kwargs...) = map(f -> parallactic_angle(read_header(f); kwargs...), filelist)

parallactic_angles(filelist::DataFrame; kwargs...) = map(r -> parallactic_angle(r; kwargs...), eachrow(filelist))

function parallactic_angle(row; latitude=19.82636)
    # compute the hour angle at the middle of the frame
    exp_start = Time(row["EXPSTART"])
    exp_stop = Time(row["EXPSTOP"])
    exp_mid = Dates.value((exp_stop - exp_start) / 2)
    # convert exp_mid from nanoseconds of time to degrees
    Δha = exp_mid / 1e9 / 24 / 10
    ha_mid = row["HA"] + Δha

    # precess the star coordinates to the appropriate epoch
    ra = deg2rad(row["RA"])
    dec = deg2rad(row["DEC"])
    coord = FK5Coords{2000}(ra, dec)
    obs_epoch = Date(row["DATE-OBS"])
    yr = year(obs_epoch) + (dayofyear(obs_epoch) - 1) / daysinyear(obs_epoch)
    coord_curr = convert(FK5Coords{yr}, coord)

    # derive true PA at the middle of the frame
    ra = coord_curr.ra
    dec = coord_curr.dec
    PA = -rad2deg(atan(-sind(ha_mid), cos(dec) * tand(latitude) - sin(dec) * cosd(ha_mid)))
    # apply instrumental effects
    PA += row["ROTPOSN"] - row["INSTANGL"] + (row["PARANG"] - row["PARANTEL"])

    return PA
end