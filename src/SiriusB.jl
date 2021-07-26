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

function load_or_produce(func, filename; force=false)


include("angles.jl")
include("fitting.jl")
include("badpix.jl")

end