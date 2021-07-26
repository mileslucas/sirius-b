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

fits = pyimport("astropy.io.fits")

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
    ext = last(splitext(filename))

    if !isfile(filename) || force
        @info "$filename not found, producing..."
        res = func()
        SAVER[ext](filename, res)
        @info "file saved to $filename"
        return res
    else
        @info "file loaded from $filename"
        return LOADER[ext](filename)
    end        
end

const LOADER = Dict(
    ".csv" => load_csv,
    ".fits" => load_fits
)
load_fits(filename; kwargs...) = Array(fits.getdata(filename; kwargs...))
load_csv(filename; kwargs...) = DataFrame(CSV.File(filename; kwargs...))

const SAVER = Dict(
    ".csv" => save_csv,
    ".fits" => save_fits
)
save_fits(filename, data; kwargs...) = fits.writeto(filename, Array(data); overwrite=true, kwargs...)
save_csv(filename, data; kwargs...) = CSV.write(data, filename; kwargs...)


include("angles.jl")
include("fitting.jl")
include("badpix.jl")

end