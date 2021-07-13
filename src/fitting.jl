using HCIToolbox
using LossFunctions
using Optim
using PSFModels: Gaussian

# generative model
function model(X::AbstractVector{T}) where T
    position = @view X[1:2] # x, y position
    fwhm     = X[3]         # fwhm
    amp      = exp(X[4])    # amplitude
    return amp * Gaussian{T}(position, fwhm)
end

# objective function
function loss(X::AbstractVector{T}, target::SubArray) where T
    # cheap way to enforce positivity
    all(>(0), X) || return T(Inf)
    # get generative model
    m = model(X)
    # l2-distance loss (χ² loss) (LossFunctions.jl)
    stamp = @view m[target.indices...]
    return value(L2DistLoss(), target, stamp, AggMode.Sum())
end

"""
    gaussian_fit(frame, guess=(256.5, 256.5))

Fit a 2D Gaussian PSF profile to `frame` centered around `guess` and return the best fitting parameters `[x, y, FWHM, log(Amp)]`.
"""
function gaussian_fit(frame::AbstractMatrix{T}, guess=(256.5, 256.5)) where T
    # crop target view for speed
    target = HCIToolbox.cropview(frame, 270)
    X0 = T[guess..., 7.5, log(maximum(target))]
    res = optimize(P -> loss(P, target), X0, LBFGS(); autodiff=:forward)
    # turn best-fit location into index
    return Optim.minimizer(res)
end

"""
    gaussian_fit_offset(frame, guess=(256.5, 256.5))

Fit a 2D Gaussian PSF profile to `frame` centered around `guess` and return the offset from the center. Combine this with `HCIToolbox.shift_frame` for registering frames.
"""
function gaussian_fit_offset(frame::AbstractMatrix{T}, guess=(256.5, 256.5)) where T
    Xbest = gaussian_fit(frame, guess)
    dx = 256.5 - Xbest[1]
    dy = 256.5 - Xbest[2]
    return dx, dy
end