module SiriusB

using CSV
using DataFrames
using Dates
using FITSIO
import ImageTransformations: center
using SkyCoords
using PyCall

export rootdir, datadir, srcdir, notebookdir, paperdir, figuredir, load_or_produce,
       parallactic_angles, gaussian_fit, gaussian_fit_offset, center

# path configuration
rootdir(args...) = joinpath(@__DIR__(), "..", args...)
datadir(args...) = rootdir("data", args...)
srcdir(args...) = rootdir("src", args...)
notebookdir(args...) = rootdir("notebooks", args...)
paperdir(args...) = rootdir("paper", args...)
figuredir(args...) = paperdir("figures", args...)

"""
    load_or_produce(f, filename; force=false)

Load the data from `filename` if it exists, otherwise (or if `force` is true) produce the data using the function `f`. This is a simplified version of the system DrWatson.jl uses for organization and reproducibility. To facilitate loading and saving, only ".fits" and ".csv" files are supported.

# Examples

```julia
using ADI
load_or_produce("flat_res_pca-2.csv") do
    PCA(2)(cube, angles)
end
```
"""
function load_or_produce(func, filename; force=false)
    filename = normpath(filename)
    ext = last(splitext(filename))

    if !isfile(filename) || force
        @info "$(normpath(filename)) not found, producing..."
        res = func()
        SAVER[ext](filename, res)
        @info "file saved to $filename"
        return res
    else
        @info "file loaded from $filename"
        return LOADER[ext](filename)
    end        
end


function load_fits(filename; kwargs...)
    PYFITS = pyimport("astropy.io.fits")
    Array(PYFITS.getdata(filename; kwargs...))
end
load_csv(filename; kwargs...) = DataFrame(CSV.File(filename; kwargs...))
const LOADER = Dict(
    ".csv" => load_csv,
    ".fits" => load_fits
)

function save_fits(filename, data; kwargs...)
    mkpath(dirname(filename))
    PYFITS = pyimport("astropy.io.fits")
    PYFITS.writeto(filename, Array(data); overwrite=true, kwargs...)
end
function save_csv(filename, data; kwargs...)
    mkpath(dirname(filename))
    CSV.write(filename, data; kwargs...)
end
const SAVER = Dict(
    ".csv" => save_csv,
    ".fits" => save_fits
)

include("angles.jl")
include("fitting.jl")
include("badpix.jl")

end