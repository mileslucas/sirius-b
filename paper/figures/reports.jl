using ADI
using ProgressLogging
using PSFModels: Gaussian
using PyCall
using SiriusB

fits = pyimport("astropy.io.fits")
pro = pyimport("proplot")

pro.use_style("ggplot")
pro.rc["image.cmap"] = "inferno"
pro.rc["image.origin"] = "lower"
pro.rc["grid"] = false

parallax = 376.6801e-3 # arcseconds
pxscale = 0.01 # arcsec / px
auscale = pxscale / parallax # AU / px

function plot_results(res_cube, angles, curve; fwhm, epoch, label)
    # produce images
    flat_res = collapse(res_cube, angles)
    snr_map = detectionmap(flat_res, fwhm)
    sig_map = detectionmap(significance, flat_res, fwhm)
    stim_map = stimmap(res_cube, angles)

    arcsec = curve.distance .* pxscale

    ticks = range(20, 179, length=5)
    tick_labs = @. string(round(pxscale * (ticks - 99.5), digits=1))
    # start plotting
    # NB: proplot does tricky things with the Axes class that doesn't play nicely
    # with PyCall.jl because it's treated as a Vector in Julia, so use Python
    # syntax instead
    py"""
    import proplot as pro

    layout = [[1, 2, 5, 5],
              [3, 4, 5, 5]]
    fig, axs = pro.subplots(
        layout,
        width="7.5in",
        wspace=("3.25em", "2.5em", "0"),
        hspace="1.25em",
        share=0
    )
    axs[0].imshow($flat_res, colorbar="r", colorbar_kw=dict(space=0))
    axs[0].text(6, 6, "flat residual", color="w")
    axs[2].imshow($stim_map, vmin=0, vmax=1, colorbar="r", colorbar_kw=dict(space=0))
    axs[2].text(6, 6, "STIM prob.", color="w")
    axs[1].imshow($snr_map, vmin=0, vmax=5, colorbar="r", colorbar_kw=dict(space=0))
    axs[1].text(6, 6, "S/N", color="w")
    axs[3].imshow($sig_map, vmin=0, vmax=5, colorbar="r", colorbar_kw=dict(space=0))
    axs[3].text(6, 6, "significance", color="w")

    axs[0].format(xticks="null")
    axs[1].format(xticks="null", yticks="null")
    axs[3].format(yticks="null")
    axs[0].format(
        yticks=$ticks,
        yticklabels=$tick_labs,
        ylabel="y [arcsec]"
    )
    axs[2].format(
        yticks=$ticks,
        yticklabels=$tick_labs,
        ylabel="y [arcsec]"
    )
    axs[2:4].format(
        xticks=$ticks,
        xticklabels=$tick_labs,
        xlabel="x [arcsec]",
    )

    axs[4].plot($arcsec, $(curve.contrast), label="Guassian", lw=1)
    axs[4].plot($arcsec, $(curve.contrast_corr), color="C0", ls="--", label="Student-t", lw=1)
    axs[4].fill_between($arcsec, $(curve.contrast), $(curve.contrast_corr), color="C0", alpha=0.2)
    axs[4].legend(ncol=1)
    axs[4].format(
        ytickloc="right",
        yformatter="log",
        xlabel="separation [arcsec]",
        ylabel="5σ contrast",
        yscale="log",
        grid=True,
        ylim=(10**(-3.5), 1)
    )
    fig.savefig($(figuredir("reports", make_filename_friendly("$(epoch)_$label") * ".pdf")))
    """
end

function make_filename_friendly(filename)
    tmp = lowercase(filename)
    tmp = replace(tmp, " " => "_")
    tmp = replace(tmp, "(" => "-")
    tmp = replace(tmp,  ")" => "")
end

_epochs = Dict("2020feb04" => "2020-02-04", "2020nov21" => "2020-11-21", "2020nov28" => "2020-11-28")

function plot_curves(contrast_curves; epoch)
    py"""
    import proplot as pro
    fig, axs = pro.subplots(width="7.5in", height="4in")
    """
    i = 0
    for (label, curve) in contrast_curves
        dist = curve.distance .* auscale
        color = "C$i"
        py"""
        axs.plot($dist, $(curve.contrast), color=$color, label=$label)
        axs.plot($dist, $(curve.contrast_corr), color=$color, ls="--")
        """
        i += 1
    end
    py"""
    axs.text(0.25, 5e-4, $(_epochs[epoch]), fontsize=14)
    axs.legend(ncol=2)
    axs.dualx(lambda x: x * $parallax, label="separation [arcsec]")
    ax2 = axs.alty(
        yticks=[0.09, 0.31, 0.54, 0.77, 1],
        yticklabels=["8", "6", "4", "2", "0"],
        label="Δ mag"
    )

    axs.format(
        ylim=(10**(-3.5), 1),
        yscale="log",
        yformatter="log",
        xlabel="projected separation [AU]",
        ylabel="5σ contrast",
        grid=True
    )
    fig.savefig($(figuredir(make_filename_friendly("$(epoch)_contrast_curves") * ".pdf")))
    """
end


@progress name="epoch" for epoch in ["2020feb04", "2020nov21", "2020nov28"]
    cube = fits.getdata(datadir("epoch_$epoch", "processed", "$(epoch)_sirius-b_cube_calib_registered_crop.fits"))
    angles = fits.getdata(datadir("epoch_$epoch", "processed", "$(epoch)_sirius-b_pa.fits"))

    raw_median = collapse(cube)
    res = gaussian_fit(raw_median, (100.5, 100.5))
    fwhm = res[3]

    psf_model = Gaussian{eltype(cube)}(res[1:2], fwhm)
    @info "FWHM for epoch $epoch is $fwhm px"
    starphot = Metrics.estimate_starphot(cube, fwhm)
    psfphot = Metrics.estimate_starphot(psf_model, fwhm) * exp(res[4])
    @info "star photometry is ≈ $starphot (PSF amp was $psfphot)"
    
    # naturally mask out inner FWHM
    av_cube = AnnulusView(cube .- minimum(cube); inner=fwhm)
    # create annular view for potential reductions
    mav_cube = MultiAnnulusView(cube .- minimum(cube), fwhm; inner=fwhm)

    # orient algorithms with appropriate data structures
    targets = [
        (Classic(), av_cube, "median"),
        (PCA(2), av_cube, "PCA(2)"),
        (NMF(2), av_cube, "NMF(2)"),
        (GreeDS(2), av_cube, "GreeDS(2)"),
        (Framewise(Classic(), delta_rot=0.5), mav_cube, "annular median"),
        (Framewise(PCA(2), delta_rot=0.5), mav_cube, "annular PCA(2)"),
        (Framewise(NMF(2), delta_rot=0.5), mav_cube, "annular NMF(2)"),
    ]

    contrast_curves = Dict{String, Any}()

    @progress "algorithm" for (alg, target, label) in targets
        @info label
        res_cube = subtract(alg, target; angles, fwhm)
        cc = contrast_curve(alg, target, angles, psf_model; fwhm, angles, starphot, nbranch=6)
        push!(contrast_curves, label => cc)
        plot_results(res_cube, angles, cc; fwhm, label, epoch)
    end

    plot_curves(contrast_curves; epoch)
end
