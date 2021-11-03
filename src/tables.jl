function get_below(subtable, arg)
    _tab = Iterators.filter(r -> r[arg.first] ≤ arg.second, eachrow(subtable))
    local y
    try
        y = maximum(r -> r[arg.first], _tab)
    catch e
        @debug "out of bounds" subtable arg
        throw(e)
    end
    return filter(r -> r[arg.first] ≈ y, subtable)
end

function get_above(subtable, arg)
    _tab = Iterators.filter(r -> r[arg.first] > arg.second, eachrow(subtable))
    local y
    try
        y = minimum(r -> r[arg.first], _tab)
    catch e
        @debug "out of bounds" subtable arg
        throw(e)
    end
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


function read_Sonora_grid(dir; metallicity=0, CO_ratio=1.0)
    # prepare file names
    if CO_ratio ≈ 1
        dex_token = @sprintf "%+02.1f" metallicity
        CO_token = ""
    else
        dex_token = @sprintf "_m%+02.1f" metallicity
        CO_token = @sprintf "_co%02.1f" CO_ratio
    end
    evo_path = joinpath(dir, "evolution_tables", @sprintf("evo_tables%+02.1f", metallicity), @sprintf("nc%+02.1f_co%02.1f_age", metallicity, CO_ratio))
    mag_path = joinpath(dir, "photometry_tables", @sprintf("mag_table%s%s", dex_token, CO_token))

    # read evolution_table
    evo_header = ["age_Gyr", "M_Msun", "logL", "Teff", "logg", "R_Rsun"]
    evo_mat = open(evo_path, "r") do fh
        # skip first (header) row
        readline(fh)
        # iterate over rows
        rows = map(readlines(fh)) do row
            tokens = split(replace(row, "\n" => ""))
            if length(tokens) == 1
                nothing
            else
                parse.(Float64, tokens)
            end
        end
        rows = filter(!isnothing, rows)
        return mapreduce(r -> reduce(hcat, r), vcat, rows)
    end
    evo_table = DataFrame(evo_mat, evo_header)
    evo_table.M_Mjup = evo_table.M_Msun .* 1047.57
    
    # read photometry_table
    mag_header = ["Teff", "logg", "M_Mjup", "R_Rsun", "Y", "log_Kzz", "MKO_Y", "MKO_Z", "MKO_J", "MKO_H", "MKO_K", 
    "MKO_Lp", "MKO_Mp", "2MASS_J", "2MASS_H", "2MASS_Ks", "Keck_Ks", "Keck_Lp", "Keck_Ms", "SDSS_g",
    "SDSS_r", "SDSS_i", "SDSS_z", "IRAC_3.6", "IRAC_4.5", "IRAC_5.7", "IRAC_7.9", "WISE_W1", "WISE_W2",
    "WISE_W3", "WISE_W4"]
    mag_mat = readdlm(mag_path, skipstart=11)
    # remove values with R/Rsun has a *
    mag_mat = @. Float64(mag_mat[mag_mat[:, 4] isa Float64, :])
    # replace Kzz values
    @. mag_mat[mag_mat[:, 6] ≈ -99, 6] = -Inf
    mag_table = DataFrame(mag_mat, mag_header)
    
    # now we need to combine the magnitude and evolution tables
    # this is accomplished by taking (Teff, logg) pairs and interpolating age from the evo table
    for key in mag_header[7:end]
        evo_table[!, key] = map(eachrow(evo_table)) do row
            try
                table_interpolate(mag_table, :Teff => row.Teff, :logg => row.logg, Symbol(key))
            catch
                missing
            end
        end
    end
    dropmissing!(evo_table)
    return evo_table
end

const SonoraSolar = read_Sonora_grid(srcdir("mass_models", "Sonora_Bobcat"))
const SonoraMetalRich = read_Sonora_grid(srcdir("mass_models", "Sonora_Bobcat"); metallicity=0.5)
