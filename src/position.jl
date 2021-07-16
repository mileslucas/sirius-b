#=
This file was used for calculating the offset of Sirius B from Sirius A. We
needed to do this because it is very hard to lock onto Sirius B directly, 
in addition to the fact we were using Sirius A as the AO NGS and targeting
Sirius B in an off-axis manner.
=#

using AstroLib
using Measurements
using SkyCoords
using Unitful, UnitfulAstro

# http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=sirius+a&submit=SIMBAD+search
sirius = ICRSCoords("06:45:08.91728", "-16:42:48.0171")
parallax = 379.21e-3 #mas
distance = 1 / parallax # pc


## orbital parameters
## https://www.stelledoppie.it/index2.php?iddoppia=27936
# a = deg2rad((7.4957 ± 0.0025) / distance / 3600) # semi-major axis
a = (7.4957 ± 0.0025)
e = (0.59142 ± 0.00037)          # eccentricity
Ω = deg2rad(45.4 ± 0.071)        # logitude of ascending node
i = deg2rad(136.336 ± 0.04)      # inclination
ω = deg2rad(149.161 ± 0.075)     # argument of periapsis
T = (1994.5715 ± 0.0058)u"yr" # periastron epoch
P = (50.1284 ± 0.0043)u"yr"   # orbital period

## Solve Kepler's equations
# mean anomaly
M(t) = ustrip(u"NoUnits", 2π * (t - T) / P)
# eccentric anomaly
E(t) = kepler_solver(M(t), e)

## equation of a conic section
# true anomaly
ν(t) = 2atan(sqrt((1 + e) / (1 - e)) * tan(E(t)/2))
r(t) = a * (1 - e^2) / (1 + e * cos(ν(t)))

pa(t) = Ω + atan(tan(ν(t) + ω) * cos(i))
sep(t) = r(t) * cos(ν(t) + ω) / cos(pa(t) - Ω)

t = 2020.8876712328767123u"yr" # time for 11/20/2020
t = 2020.9094227697923797u"yr" # time for 11/27/2020
ρ_arcsec = sep(t)
θ = pa(t)

ρ = deg2rad(ρ_arcsec / 3600)
θ_deg = rad2deg(θ)


println("ρ= $ρ_arcsec\", θ= $(θ_deg)°")
println("offset: $(ρ_arcsec .* sincosd(θ_deg))")


# sirius_b = offset(sirius, Measurements.value(ρ), Measurements.value(θ))

# dra = sirius_b.ra - sirius.ra
# # ra = sirius_b.ra
# ra_hrs = 24 * dra / 2π
# ra_h, ra_r = divrem(sirius_b.ra * 24 / 2π, 1)
# ra_m, ra_r = divrem(ra_r * 60, 1)
# ra_s = ra_r * 60

# ddec = sirius_b.dec - sirius.dec
# # dec = sirius_b.dec
# dec_degs = rad2deg(ddec)
# dec_d, dec_r = divrem(rad2deg(sirius_b.dec), 1)
# dec_m, dec_r = divrem(dec_r * 60, 1)
# dec_s = dec_r * 60


# println("ρ= $ρ_arcsec\", θ= $(θ_deg)°")
# println("sirius b = $(ra_h)h$(ra_m)m$(ra_s)s $(dec_d):$(dec_m):$(dec_s)")

# println("Δ(ra, dec)= $(ra_hrs * 24 * 3600 / 2π )s $(dec_degs * 3600)\"")
# println("Δ(ra, dec)= $(rad2deg(ra_hrs) * 3600)\" $(dec_degs * 3600)\"")
