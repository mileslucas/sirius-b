function get_below(subtable, arg)
    _tab = Iterators.filter(r -> r[arg.first] ≤ arg.second, eachrow(subtable))
    y = maximum(r -> r[arg.first], _tab)
    return filter(r -> r[arg.first] ≈ y, subtable)
end

function get_above(subtable, arg)
    _tab = Iterators.filter(r -> r[arg.first] > arg.second, eachrow(subtable))
    y = minimum(r -> r[arg.first], _tab)
    return filter(r -> r[arg.first] ≈ y, subtable)
end


function table_interpolate(table, arg1, arg2, arg3)
    sub = copy(table)

    sub0 = get_below(sub, arg1)
    x0 = x1 = sub0[1, arg1.first]

    tab0 = get_below(sub0, arg2)
    y0 = tab0[1, arg2.first]
    z0 = tab0[1, arg3]

    tab1 = get_above(sub0, arg2)
    y1 = tab1[1, arg2.first]
    z1 = tab1[1, arg3]

    sub1 = get_above(sub, arg1)
    x2 = x3 = sub1[1, arg1.first]

    tab2 = get_below(sub1, arg2)
    y2 = tab2[1, arg2.first]
    z2 = tab2[1, arg3]

    tab3 = get_above(sub1, arg2)
    y3 = tab3[1, arg2.first]
    z3 = tab3[1, arg3]

    xs = [x0, x1, x2, x3]
    ys = [y0, y1, y2, y3]
    zs = [z0, z1, z2, z3]
    spl = Spline2D(xs, ys, zs, kx=1, ky=1)
    return spl(arg1.second, arg2.second)
end

#################################################

function read_ATMO2020_grid(dir)
    filenames = filter(f -> endswith(f, "vega.txt"), readdir(dir))

    header = ["Mass", "Age", "Teff", "Luminosity", "Radius", "Gravity", "MKO_Y", "MKO_J", "MKO_H", 
              "MKO_K", "MKO_Lp", "MKO_Mp", "W1", "W2", "W3", "W4", "IRAC_CH1", "IRAC_CH2"]
    
    full_table = mapreduce(vcat, filenames) do filename
        mat = readdlm(joinpath(dir, filename), skipstart=2)
        # a bunch of magnitudes listed with 0, treat as NaN
        @. mat[iszero(mat)] = NaN
        return DataFrame(mat, header)
    end
    full_table.Mass_MJ = full_table.Mass .* 1047.57
    return full_table
end

const ATMO2020 = read_ATMO2020_grid(srcdir("mass_models", "ATMO2020_NEQ_weak"))

#################################################


