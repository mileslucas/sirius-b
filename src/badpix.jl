using ImageFiltering
using Statistics
using StatsBase

function badpixremoval_clump(frame::AbstractMatrix{T}; max_niter=15, kwargs...) where T
    # create bad pixel map by sigma clipping outliers
    bad_pixels = badpixel_clip(frame; kwargs...)
    bad_pixels_cumul = copy(bad_pixels)
    total = sum(bad_pixels)

    # iterate replacing bad pixels with sigma clipping
    out = copy(frame)
    it = 0
    while total > 0 && it < max_niter
        it += 1
        out = sigma_filter!(out, bad_pixels)
        bad_pixels = badpixel_clip(out)
        bad_pixels_cumul .|= bad_pixels
        total = sum(bad_pixels)
    end
    return out, bad_pixels_cumul
end

function badpixel_clip(frame; sigma=4, box_size=9)
    med = mapwindow(median, frame, (box_size, box_size); border="reflect")
    std = mapwindow(mad, frame, (box_size, box_size); border="reflect")

    return @. abs(frame - med) > sigma * std
end

function sigma_filter!(frame, bad_pixels; box_size=9)
    bpm = copy(bad_pixels)

    min_neigh = sum(3:2:box_size)
    half_box = box_size ÷ 2
    while sum(bpm) > 0
        pad_bpm = BorderArray(bpm, Pad(:reflect, half_box, half_box))
        pad_im = BorderArray(frame, Pad(:reflect, half_box, half_box))
        for idx in findall(bpm)
            ax1 = idx.I[1] - half_box:idx.I[1] + half_box
            ax2 = idx.I[2] - half_box:idx.I[2] + half_box
            block_bpm = @view pad_bpm[ax1, ax2]

            # ensure enough "good" neighbors for imputation
            if sum(block_bpm) < min_neigh
                block_im = @view pad_im[ax1, ax2]
                frame[idx] = median(block_im)
                bpm[idx] = false
            end
        end
    end
    return frame
end