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


function plot_mosaic(res_cubes, angles; epoch, label, ncomps)
    py"""
    import proplot as pro
    fig, axs = pro.subplots(nrows=2, ncols=5, share=3, width="7.5in", space=0)
    """
    i = 0
    for (n, res) in zip(res_cubes, ncomps)
        lab = "ncomp=$n"
        flat = collapse(res, angles)
        py"""
        axs[$i].imshow($flat)
        axs[$i].text(6, 6, $lab, color="w")
        """
        i += 1
    end

    ticks = range(40, 159, length=3)
    tick_labs = @. string(round(pxscale * (ticks - 99.5), digits=1))

    py"""
    axs.format(
        yticks=$ticks,
        yticklabels=$tick_labs,
        ylabel="y [arcsec]"
        xticks=$ticks,
        xticklabels=$tick_labs,
        xlabel="x [arcsec]",
    )
    fig.savefig($(figuredir(make_filename_friendly("$(epoch)_$(label)_mosaic") * ".pdf")))
    """
end

function plot_results(res_cubes, angles, contrast_curves; epoch, label, fwhm)
    py"""
    import proplot as pro
    fig, axs = pro.subplots([[1, 3], [2, 3]], share=0, width="7.5in")
    cycle = pro.Cycle("viridis", $(length(res_cubes)))
    lines = []
    """
    # plot the average STIM map and SLIM probability map
    ticks = range(20, 179, length=5)
    tick_labs = @. string(round(pxscale * (ticks - 99.5), digits=1))

    stimmap, slimmask = slimmap(res_cubes, angles; N=cil(Int, π/4 * fwhm^2))
    slimprob = stimmap .* slimmask
    py"""
    axs[0].imshow($stimmap, vmin=0, vmax=1)
    axs[0].text(6, 6, "mean STIM", color="w")
    m = axs[1].imshow($slimprob, vmin=0, vmax=1)
    axs[1].text(6, 6, "SLIM map", color="w")
    fig.colorbar(m, loc="l")

    axs[0].format(
        ylabel="y [arcsec]"
        yticks=$ticks
        yticklabels=$tick_labs
        xticks="null"
    )
    axs[1].format(
        ylabel="y [arcsec]"
        xlabel="x [arcsec]"
        yticks=$ticks
        yticklabels=$tick_labs
        xticks=$ticks
        yticklabels=$tick_labs
    )
    """
    # now plot contrast curves
    i = 0
    for (label, curve) in contrast_curves
        dist = curve.distance .* auscale
        color = "C$i"
        py"""
        l = axs[2].plot($dist, $(curve.contrast), cycle=cycle, color=$color)
        lines.append(l)
        """
        i += 1
    end
    py"""
    axs[2].text(0.25, 5e-4, $(_epochs[epoch]), fontsize=14)
    fig.colorbar(lines, loc="r", label="ncomp")
    axs.dualx(lambda x: x * $parallax, label="separation [arcsec]")
    ax2 = axs.alty(
        yticks=[0.09, 0.31, 0.54, 0.77, 1],
        yticklabels=["8", "6", "4", "2", "0"],
        label="Δ mag"
    )

    axs[2].format(
        ylim=(10**(-3.5), 1),
        yscale="log",
        yformatter="log",
        xlabel="projected separation [AU]",
        ylabel="5σ contrast",
        grid=True
    )
    fig.savefig($(figuredir(make_filename_friendly("$(epoch)_$(label)_contrast_curves") * ".pdf")))
    """
end

function make_filename_friendly(filename)
    tmp = lowercase(filename)
    tmp = replace(tmp, " " => "_")
    tmp = replace(tmp, "(" => "-")
    tmp = replace(tmp,  ")" => "")
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
    @info "star photometry is ≈ $starphot (PSF photometry was $psfphot)"
    
    # naturally mask out inner FWHM
    av_cube = AnnulusView(cube; inner=fwhm)

    ncomps = 1:10
    @progress name="algorithm" for alg in [PCA, NMF, GreeDS]
        @withprogress name="ncomp" begin
        i = 1
        N = length(ncomps)
        outputs = map(ncomps) do n
            res_cube = subtract(alg(n), av_cube; angles, fwhm)
            cc = contrast_curve(alg, target, angles, psf_model; fwhm, angles, starphot, nbranches=6)
            @logprogress i/N
            i += 1
            res_cube, cc
        end
        res_cubes = map(r -> r[1], outputs)
        contrast_curves = map(r -> r[1], outputs)
        label = _label(alg)
        plot_mosaic(res_cubes, angles; epoch, label, ncomps)
        plot_results(res_cubes, angles, contrast_curves; epoch, label, fwhm)
    end
end

_label(a::Type{<:PCA}) = "PCA"
_label(a::Type{<:NMF}) = "NMF"
_label(a::Type{<:GreeDS}) = "GreeDS"