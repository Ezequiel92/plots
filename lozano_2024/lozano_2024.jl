cd(@__DIR__)

using Pkg

Pkg.activate("../")
Pkg.instantiate()

using CSV,
    CairoMakie,
    ColorSchemes,
    Colors,
    DataFrames,
    DelimitedFiles,
    FileIO,
    GLM,
    Glob,
    HDF5,
    Images,
    JLD2,
    LaTeXStrings,
    LinearAlgebra,
    Logging,
    Measurements,
    NearestNeighbors,
    PrettyTables,
    ProgressMeter,
    QuadGK,
    Rotations,
    Statistics,
    StatsBase,
    Unitful,
    UnitfulAstro

push!(LOAD_PATH, "../../codes/GalaxyInspector/src/")
using GalaxyInspector

function lozano2024(
    simulation_paths::Vector{String},
    labels::Vector{String},
    base_out_path::String,
    r1::Unitful.Length,
    r2::Unitful.Length,
    r3::Unitful.Length,
    logging::Bool,
)::Nothing

    # Set the paths
    Au6_MOL_path = simulation_paths[1]
    Au6_BLT_path = simulation_paths[2]
    Au6_STD_path = simulation_paths[3]
    Au6_MOL_LR_path = simulation_paths[4]

    # Set the labels
    Au6_MOL_label = labels[1]
    Au6_BLT_label = labels[2]
    Au6_STD_label = labels[3]
    Au6_MOL_LR_label = labels[4]

    # Construct the necessary folders and files
    figures_path = mkpath(joinpath(base_out_path, "figures"))
    report_path = mkpath(joinpath(base_out_path, "reports"))
    info_file = open(joinpath(report_path, "info.txt"), "w")
    table_file = open(joinpath(report_path, "table.txt"), "w")
    log_file = open(joinpath(report_path, "logs.txt"), "w+")

    GalaxyInspector.setLogging!(logging; stream=log_file)

    # Select the last snapshot for Au6_MOL
    Au6_MOL_slice = GalaxyInspector.countSnapshot(Au6_MOL_path)

    ################################################################################################
    # Main simulation (Au6_MOL)
    ################################################################################################

    ####################################################################
    # Stellar circularity vs. stellar xy distance (scatter density map)
    ####################################################################

    scatterDensityMap(
        [Au6_MOL_path],
        Au6_MOL_slice,
        :stellar_xy_distance,
        :stellar_circularity,
        :stellar_mass,
        u"Msun";
        x_range=(0.0u"kpc", r1 / 2.0),
        y_range=(-2.0, 2.0),
        n_bins=200,
        output_path=joinpath(figures_path, Au6_MOL_label),
        filter_mode=:subhalo,
        theme=Theme(
            size=(800, 500),
            figure_padding=(15, 20, 10, 20),
            Axis=(aspect=AxisAspect(1.7),),
        ),
    )

    ##################################
    # Stellar surface density profile
    ##################################

    densityProfile(
        [Au6_MOL_path],
        Au6_MOL_slice,
        :stellar_mass;
        yscale=log10,
        radius=r1,
        n_bins=40,
        output_path=joinpath(figures_path, Au6_MOL_label),
        filter_mode=:subhalo,
        sim_labels=nothing,
        theme=Theme(
            palette=(color=[Makie.wong_colors()[2]],),
            size=(800, 470),
            figure_padding=(0, 15, 5, 15),
            Axis=(aspect=AxisAspect(1.7),),
        ),
    )

    ######################################################################################
    # SFR vs physical time (radio separated using an stellar age histogram at redshift 0)
    ######################################################################################

    da_ff = [
        dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r1), :zero),
        dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r3), :zero),
        dd -> GalaxyInspector.filterWithinSphere(dd, (r3, r1), :zero),
    ]

    x_plot_params = GalaxyInspector.plotParams(:physical_time)
    y_plot_params = GalaxyInspector.plotParams(:sfr)

    filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
        :subhalo,
        Dict(:stars => ["GAGE"]),
    )

    r1_label = string(round(Int, ustrip(u"kpc", r1)))
    r3_label = string(round(Int, ustrip(u"kpc", r3)))

    plotSnapshot(
        [Au6_MOL_path, Au6_MOL_path, Au6_MOL_path],
        request,
        [lines!];
        output_path=joinpath(figures_path, Au6_MOL_label),
        base_filename="sfr_vs_physical_time_radio_separated",
        output_format=".pdf",
        slice=Au6_MOL_slice,
        filter_function,
        da_functions=[GalaxyInspector.daStellarHistory],
        da_kwargs=[
            (; filter_function=da_ff[1]),
            (; filter_function=da_ff[2]),
            (; filter_function=da_ff[3]),
        ],
        transform_box=true,
        translation,
        rotation,
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=y_plot_params.exp_factor,
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        yaxis_scale_func=log10,
        theme=Theme(
            size=(800, 470),
            figure_padding=(1, 20, 5, 15),
            palette=(
                color=[:gray65, :orangered2, :navy],
                linestyle=[:solid],
            ),
            Axis=(
                aspect=AxisAspect(1.7),
                limits=(nothing, (1.0e-2, nothing)),
                xticks=0:2:14,
            ),
            Legend=(
                nbanks=1,
                halign=:right,
                valign=:bottom,
                padding=(0, 45, 5, 0),
                labelsize=25,
                rowgap=-4,
            ),
        ),
        sim_labels=[
            L"r \,\, \le \,\, %$(r1_label) \, \mathrm{kpc}",
            L"r \,\, \le \,\, %$(r3_label) \, \mathrm{kpc}",
            L"%$(r3_label) \, \mathrm{kpc} \,\, < \,\, r \,\, \le \,\, %$(r1_label) \, \mathrm{kpc}",
        ],
    )

    ##########################################
    # Circularity histogram (radio separated)
    ##########################################

    plot_params = GalaxyInspector.plotParams(:stellar_circularity)
    filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
        :subhalo,
        plot_params.request,
    )
    grid = GalaxyInspector.LinearGrid(-1.5, 1.5, 200)

    da_ff = [
        dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r1), :zero),
        dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r3), :zero),
        dd -> GalaxyInspector.filterWithinSphere(dd, (r3, r1), :zero),
    ]

    r1_label = string(round(Int, ustrip(u"kpc", r1)))
    r3_label = string(round(Int, ustrip(u"kpc", r3)))

    plotSnapshot(
        [Au6_MOL_path, Au6_MOL_path, Au6_MOL_path],
        request,
        [lines!];
        output_path=joinpath(figures_path, Au6_MOL_label),
        base_filename="circularity_histogram_radio_separated",
        output_format=".pdf",
        slice=Au6_MOL_slice,
        filter_function,
        da_functions=[GalaxyInspector.daLineHistogram],
        da_args=[(:stellar_circularity, grid, :stars)],
        da_kwargs=[
            (; filter_function=da_ff[1], norm=23760),
            (; filter_function=da_ff[2], norm=23760),
            (; filter_function=da_ff[3], norm=23760),
        ],
        transform_box=true,
        translation,
        rotation,
        x_unit=plot_params.unit,
        x_exp_factor=plot_params.exp_factor,
        xaxis_label=plot_params.axis_label,
        xaxis_var_name=plot_params.var_name,
        yaxis_var_name=L"\mathrm{Normalized \,\, counts}",
        theme=Theme(
            size=(800, 470),
            figure_padding=(5, 15, 5, 15),
            palette=(
                color=[:gray65, :orangered2, :navy],
                linestyle=[:solid],
            ),
            Axis=(aspect=AxisAspect(1.7),),
            Legend=(
                nbanks=1,
                halign=:left,
                valign=:top,
                padding=(55, 0, 0, 10),
                labelsize=25,
                rowgap=-4,
            ),
        ),
        sim_labels=[
            L"r \,\, \le \,\, %$(r1_label) \, \mathrm{kpc}",
            L"r \,\, \le \,\, %$(r3_label) \, \mathrm{kpc}",
            L"%$(r3_label) \, \mathrm{kpc} \,\, < \,\, r \,\, \le \,\, %$(r1_label) \, \mathrm{kpc}",
        ],
    )

    ###############################################################################
    # Stellar density maps (face-on and edge-on projections, with velocity fields)
    ###############################################################################

    temp_folder = joinpath(figures_path, Au6_MOL_label, "_stellar_density_maps")

    box_size = 50.0u"kpc"

    densityMapVelField(
        [Au6_MOL_path],
        Au6_MOL_slice;
        quantities=[:stellar_mass],
        types=[:particles],
        output_path=temp_folder,
        filter_mode=:all_subhalo,
        box_size,
        pixel_length=box_size / 400.0,
        theme=Theme(
            size=(880, 620),
            figure_padding=(1, 1, 20, 20),
            Colorbar=(
                ticks=4:1:11,
                label=L"\mathrm{log}_{10} \, \Sigma_\star \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
            ),
            Axis=(
                xminorticksvisible=false,
                xticksvisible=false,
                xminorgridvisible=false,
                xlabelvisible=false,
                xticklabelsvisible=false,
                yticks=[-20, -10, 0, 10, 20],
            ),
        ),
        colorbar=true,
        colorrange=(4.0, 11.0),
    )

    densityMapVelField(
        [Au6_MOL_path],
        Au6_MOL_slice;
        quantities=[:stellar_mass],
        types=[:particles],
        output_path=temp_folder,
        filter_mode=:all_subhalo,
        projection_planes=[:xz],
        box_size,
        pixel_length=box_size / 400.0,
        theme=Theme(
            size=(880, 380),
            figure_padding=(15, 154, 1, 1),
            Axis=(
                limits=(-ustrip(u"kpc", box_size) / 2, ustrip(u"kpc", box_size) / 2, -10.0, 10.0),
                aspect=DataAspect(),
                xticks=[-20, -10, 0, 10, 20],
            ),
        ),
        colorrange=(4.0, 11.0),
    )

    GalaxyInspector.hvcatImages(
        1,
        glob("$(Au6_MOL_label)_stellar_mass_x*_density_map_snap_*.png", temp_folder);
        output_path=joinpath(
            figures_path,
            Au6_MOL_label,
            "stellar_density_maps.png",
        ),
    )

    rm(temp_folder; recursive=true)

    #########################################################################
    # Gas components density map (face-on/edge-on projections at redshift 0)
    #########################################################################

    temp_folder = joinpath(figures_path, Au6_MOL_label, "_density_maps")

    grid = GalaxyInspector.CubicGrid(65.0u"kpc", 400)
    projections = [:xy, :xz]
    quantities = [:hydrogen_mass, :molecular_mass, :atomic_mass, :ionized_mass]

    colorbar_labels = [
        L"\log_{10} \, \Sigma_\mathrm{H} \,\, [\mathrm{M_\odot \, kpc^{-2}}]"
        L"\log_{10} \, \Sigma_\mathrm{H_2} \,\, [\mathrm{M_\odot \, kpc^{-2}}]"
        L"\log_{10} \, \Sigma_\mathrm{HI} \,\, [\mathrm{M_\odot \, kpc^{-2}}]"
        L"\log_{10} \, \Sigma_\mathrm{HII} \,\, [\mathrm{M_\odot \, kpc^{-2}}]"
    ]

    x_label = GalaxyInspector.getLabel("x", 0, u"kpc")
    y_label = GalaxyInspector.getLabel("y", 0, u"kpc")
    z_label = GalaxyInspector.getLabel("z", 0, u"kpc")
    n_rows = length(projections)
    n_cols = length(quantities)
    x_size = 1700
    y_size = (x_size / n_cols) * n_rows + 120
    colorrange = (3.0, 9.0)

    paths = joinpath.(
        temp_folder,
        [
            "hydrogen_mass_xy.jld2",
            "molecular_mass_xy.jld2",
            "atomic_mass_xy.jld2",
            "ionized_mass_xy.jld2",
            "hydrogen_mass_xz.jld2",
            "molecular_mass_xz.jld2",
            "atomic_mass_xz.jld2",
            "ionized_mass_xz.jld2",
        ]
    )

    for projection_plane in projections

        for quantity in quantities

            filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
                :all_subhalo,
                GalaxyInspector.plotParams(quantity).request,
            )

            GalaxyInspector.plotSnapshot(
                [Au6_MOL_path],
                request,
                [heatmap!];
                output_path=temp_folder,
                base_filename="$(quantity)_$(projection_plane)",
                slice=Au6_MOL_slice,
                filter_function,
                da_functions=[GalaxyInspector.daDensity2DProjection],
                da_args=[(grid, quantity, :cells)],
                da_kwargs=[(; projection_plane)],
                transform_box=true,
                translation,
                rotation,
                x_unit=u"kpc",
                y_unit=u"kpc",
                save_figures=false,
                backup_results=true,
            )

        end

    end

    with_theme(merge(theme_latexfonts(), GalaxyInspector.DEFAULT_THEME)) do

        f = Figure(size=(x_size, y_size), figure_padding=(0, 0, 0, 0))

        for (idx, path) in enumerate(paths)

            row = ceil(Int, idx / n_cols)
            col = mod1(idx, n_cols)

            xaxis_v = row == 2
            yaxis_v = col == 1

            filename = split(basename(path), '.')[1]

            ax = CairoMakie.Axis(
                f[row+1, col];
                xlabel=x_label,
                ylabel=(row == 1 ? y_label : z_label),
                xminorticksvisible=xaxis_v,
                xticksvisible=xaxis_v,
                xlabelvisible=xaxis_v,
                xticklabelsvisible=xaxis_v,
                yminorticksvisible=yaxis_v,
                yticksvisible=yaxis_v,
                ylabelvisible=yaxis_v,
                yticklabelsvisible=yaxis_v,
                xticklabelsize=28,
                yticklabelsize=28,
                xticks=[-30, -20, -10, 0, 10, 20, 30],
                yticks=[-30, -20, -10, 0, 10, 20, 30],
            )

            jldopen(path, "r") do jld2_file

                x, y, z = jld2_file[first(keys(jld2_file))]["simulation_001"]

                pf = heatmap!(ax, x, y, z; colorrange)

                if row == 1

                    Colorbar(
                        f[row, col],
                        pf,
                        label=colorbar_labels[col],
                        ticklabelsize=23,
                        ticks=3:1:9,
                        vertical=false,
                    )

                end

            end

            colgap!(f.layout, 50)
            colsize!(f.layout, col, Makie.Fixed(352.0f0))

        end

        Makie.save(joinpath(figures_path, Au6_MOL_label, "gas_density_maps_grid.png"), f)

    end

    rm(temp_folder; recursive=true)

    #############################################################################################
    # Profiles of the gas surface density and of the fractions, for the different gas components
    #############################################################################################

    temp_folder = joinpath(figures_path, Au6_MOL_label, "_gas_profiles")

    quantities = [
        :gas_mass,
        :ionized_mass,
        :atomic_mass,
        :molecular_mass,
        :stellar_mass,
    ]
    plot_params = GalaxyInspector.plotParams(:generic_area_density)
    filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
        :subhalo,
        plot_params.request,
    )

    grid = GalaxyInspector.CircularGrid(r2, 50)

    plotSnapshot(
        fill(Au6_MOL_path, length(quantities)),
        request,
        [lines!];
        output_path=temp_folder,
        base_filename="gas_density_profiles",
        slice=Au6_MOL_slice,
        filter_function,
        da_functions=[GalaxyInspector.daProfile],
        da_args=[(quantity, grid) for quantity in quantities],
        da_kwargs=[(; flat=true, total=true, cumulative=false, density=true)],
        transform_box=true,
        translation,
        rotation,
        x_unit=u"kpc",
        y_unit=plot_params.unit,
        save_figures=false,
        backup_results=true,
    )

    quantities = [:ionized_mass, :atomic_mass, :molecular_mass]
    plot_params = GalaxyInspector.plotParams(:generic_fraction)
    filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
        :subhalo,
        plot_params.request,
    )

    grid = GalaxyInspector.CircularGrid(r2, 50)

    plotSnapshot(
        fill(Au6_MOL_path, length(quantities)),
        request,
        [lines!];
        output_path=temp_folder,
        base_filename="gas_fractions_profiles",
        slice=Au6_MOL_slice,
        filter_function,
        da_functions=[GalaxyInspector.daProfile],
        da_args=[(quantity, grid) for quantity in quantities],
        da_kwargs=[(; flat=true, fractions=true)],
        transform_box=true,
        translation,
        rotation,
        x_unit=u"kpc",
        y_unit=plot_params.unit,
        save_figures=false,
        backup_results=true,
    )

    paths = joinpath.(temp_folder, ["gas_density_profiles.jld2", "gas_fractions_profiles.jld2"])

    with_theme(merge(theme_latexfonts(), GalaxyInspector.DEFAULT_THEME)) do

        f = Figure(size=(880, 1200),)

        ax_1 = CairoMakie.Axis(
            f[1, 1];
            xlabel=L"r \, [\mathrm{kpc}]",
            ylabel=L"\Sigma \, [\mathrm{M_\odot \, pc^{-2}}]",
            yscale=log10,
            xminorticksvisible=false,
            xticksvisible=false,
            xlabelvisible=false,
            xticklabelsvisible=false,
        )

        jldopen(paths[1], "r") do jld2_file

            first_address = first(keys(jld2_file))

            x_h, y_h = jld2_file[first_address]["simulation_001"]
            x_i, y_i = jld2_file[first_address]["simulation_002"]
            x_a, y_a = jld2_file[first_address]["simulation_003"]
            x_m, y_m = jld2_file[first_address]["simulation_004"]
            x_s, y_s = jld2_file[first_address]["simulation_005"]

            lines!(ax_1, x_h, y_h; color=:black, label="Gas density")
            lines!(ax_1, x_i, y_i; color=Makie.wong_colors()[1], label="Ionized density")
            lines!(ax_1, x_a, y_a; color=Makie.wong_colors()[4], label="Atomic density")
            lines!(ax_1, x_m, y_m; color=Makie.wong_colors()[3], label="Molecular density")
            lines!(ax_1, x_s, y_s; color=Makie.wong_colors()[2], label="Stellar density")

        end

        axislegend(ax_1, position=:rt, framevisible=false, fontsize=20, nbanks=1)

        ax_2 = CairoMakie.Axis(
            f[2, 1];
            xlabel=L"r \, [\mathrm{kpc}]",
            ylabel=L"\log_{10} \, f",
            aspect=nothing,
            xticks=0:10:100,
            limits=(-2, nothing, nothing, nothing),
        )

        jldopen(paths[2], "r") do jld2_file

            first_address = first(keys(jld2_file))

            x_i, y_i = jld2_file[first_address]["simulation_001"]
            x_a, y_a = jld2_file[first_address]["simulation_002"]
            x_m, y_m = jld2_file[first_address]["simulation_003"]

            lines!(ax_2, x_i, log10.(y_i); color=Makie.wong_colors()[1])
            lines!(ax_2, x_a, log10.(y_a); color=Makie.wong_colors()[4])
            lines!(ax_2, x_m, log10.(y_m); color=Makie.wong_colors()[3])

        end

        linkxaxes!(ax_1, ax_2)
        rowsize!(f.layout, 1, Relative(2 / 3))
        colsize!(f.layout, 1, Makie.Fixed(pixelarea(ax_1.scene)[].widths[2]))

        Makie.save(
            joinpath(figures_path, Au6_MOL_label, "gas_density_and_fractions_profiles.pdf"),
            f,
        )

    end

    rm(temp_folder; recursive=true)

    ###############################################################
    # Stellar metallicity and birth times histograms, at t = 1 Gyr
    ###############################################################

    time_list = [0.0, 0.65, 0.75, 0.85, 1.0] .* u"Gyr"

    last_time = ustrip(u"Gyr", time_list[end])

    filter_function, _, _, request = GalaxyInspector.selectFilter(
        :subhalo,
        Dict(:stars => ["GZ2 ", "GAGE"]),
    )

    data_dict = makeDataDict(Au6_MOL_path, 44, request)

    GalaxyInspector.filterData!(data_dict; filter_function)

    metallicities = GalaxyInspector.daStellarMetallictyHistogram(data_dict)[1]
    birth_ticks = data_dict[:stars]["GAGE"]
    birth_times = GalaxyInspector.computeTime(birth_ticks, data_dict[:snap_data].header)

    with_theme(merge(theme_latexfonts(), GalaxyInspector.DEFAULT_THEME)) do

        f = Figure(size=(880, 880), figure_padding=(0, 0, 5, 25),)

        ax = CairoMakie.Axis(
            f[1, 1];
            xlabel=L"\log_{10} \, Z_\star \, / \, Z_\odot",
            ylabel="counts",
            yscale=log10,
            limits=(nothing, nothing, exp10(0), exp10(3)),
            xticks=(
                [-4.9, -4, -3, -2, -1, 0],
                [
                    rich(
                        rich("Z", font=:italic, color=:blue),
                        subscript("⋆", offset=(0.15, 0)),
                        " = 0",
                        color=:blue,
                    ),
                    "-4",
                    "-3",
                    "-2",
                    "-1",
                    "0",
                ],
            ),
        )

        grid = GalaxyInspector.LinearGrid(exp10(-4), last_time, 25; log=true)

        plot_grid = GalaxyInspector.LinearGrid(exp10(-5), last_time, 30; log=true)

        ##############################
        # Metalic stars and INFO FILE
        ##############################

        println(info_file, "Stellar history up to $(last_time) Gyr:\n")

        zero_stars = Vector{Int64}(undef, length(time_list) - 1)

        for i in eachindex(time_list)[1:(end-1)]

            idxs = findall(t -> time_list[end-i] <= t < time_list[end-(i-1)], birth_times)

            isempty(idxs) && continue

            t1 = round(ustrip(u"Gyr", time_list[end-i]); sigdigits=3)
            t2 = round(ustrip(u"Gyr", time_list[end-(i-1)]); sigdigits=3)

            zero_stars[i] = count(iszero, metallicities[idxs])

            println(info_file, "\tTime range: $(t1) Gyr - $(t2) Gyr:\n")
            println(info_file, "\t\tNumber of stars:              $(length(idxs))\n")
            println(info_file, "\t\tNumber of stars with Z = 0:   $(zero_stars[i])\n")
            println(
                info_file,
                "\t\tExtrema of metallicity [Z⊙]: $(extrema(metallicities[idxs]))\n\n",
            )

            counts = [
                fill(NaN, length(plot_grid.grid) - length(grid.grid))...,
                GalaxyInspector.histogram1D(metallicities[idxs], grid)...,
            ]

            barplot!(
                ax,
                log10.(plot_grid.grid),
                counts;
                label=L"%$(t1) \, \mathrm{Gyr} \leq t \, < \, %$(t2) \, \mathrm{Gyr}",
                color=Makie.wong_colors()[i],
                strokecolor=Makie.wong_colors()[i],
                direction=:y,
                strokewidth=1,
                bar_labels=nothing,
                dodge_gap=0.0,
                gap=0.0,
                fillto=0.5,
            )

        end

        ###################
        # Metal free stars
        ###################

        nans = fill(NaN, length(plot_grid.grid) - 1)

        for i in eachindex(zero_stars)[end:-1:1]

            zs = zero_stars[i]

            barplot!(
                ax,
                log10.(plot_grid.grid),
                [iszero(zs) ? NaN : zs, nans...];
                color=Makie.wong_colors()[i],
                strokecolor=Makie.wong_colors()[i],
                direction=:y,
                strokewidth=1,
                bar_labels=nothing,
                dodge_gap=0.0,
                gap=0.0,
                fillto=0.5,
            )

        end

        println(info_file, "#"^100, "\n")

        axislegend(ax, position=:lt, nbanks=1)

        Makie.save(joinpath(figures_path, Au6_MOL_label, "metallicity_history_histogram.pdf"), f)

        f = Figure()

        ax = CairoMakie.Axis(
            f[1, 1];
            xlabel=L"\mathrm{Birth \,\, time \, [Gyr]}",
            ylabel="counts",
            yscale=log10,
            xticks=0.4:0.2:2.0,
        )

        grid = GalaxyInspector.LinearGrid(time_list[1], time_list[end], 40; log=false)

        for i in eachindex(time_list)[1:(end-1)]

            bt = filter(t -> time_list[end-i] <= t < time_list[end-(i-1)], birth_times)

            t1 = round(ustrip(u"Gyr", time_list[end-i]); sigdigits=2)
            t2 = round(ustrip(u"Gyr", time_list[end-(i-1)]); sigdigits=2)

            counts = Float64.(GalaxyInspector.histogram1D(bt, grid))

            replace!(counts, 0.0 => NaN)

            barplot!(
                ax,
                ustrip.(u"Gyr", grid.grid),
                counts;
                label=L"%$(t1) \, \mathrm{Gyr} \leq t \, < \, %$(t2) \, \mathrm{Gyr}",
                color=Makie.wong_colors()[i],
                strokecolor=Makie.wong_colors()[i],
                direction=:y,
                strokewidth=1,
                bar_labels=nothing,
                dodge_gap=0.0,
                gap=0.0,
                fillto=0.5,
            )

        end

        axislegend(ax, position=:lt, nbanks=1)

        Makie.save(joinpath(figures_path, Au6_MOL_label, "birth_time_history_histogram.pdf"), f)

    end

    ##########################################################################################
    # Gas components density maps (face-on projection evolution, from t = 2 Gyr to t = 5 Gyr)
    ##########################################################################################

    temp_folder = joinpath(figures_path, Au6_MOL_label, "_density_maps_evolution")
    simulation_table = GalaxyInspector.makeSimulationTable(Au6_MOL_path)

    grid = GalaxyInspector.CubicGrid(50.0u"kpc", 400)
    snap_list = [54, 60, 67, 73]
    quantities = [:stellar_mass, :molecular_mass, :atomic_mass, :ionized_mass]
    types = [:particles, :cells, :cells, :cells]
    times = ustrip.(u"Gyr", simulation_table[snap_list, :physical_times])

    colorbar_ranges = [(5.0, 11.0), (5.0, 11.0), (5.0, 11.0), (5.0, 11.0)]
    colorbar_ticks = [5:1:11, 5:1:11, 5:1:11, 5:1:11]
    colorbar_labels = [
        L"\log_{10} \, \Sigma_\star \,\, [\mathrm{M_\odot \, kpc^{-2}}]"
        L"\log_{10} \, \Sigma_\mathrm{H_2} \,\, [\mathrm{M_\odot \, kpc^{-2}}]"
        L"\log_{10} \, \Sigma_\mathrm{HI} \,\, [\mathrm{M_\odot \, kpc^{-2}}]"
        L"\log_{10} \, \Sigma_\mathrm{HII} \,\, [\mathrm{M_\odot \, kpc^{-2}}]"
    ]

    x_label = GalaxyInspector.getLabel("x", 0, u"kpc")
    y_label = GalaxyInspector.getLabel("y", 0, u"kpc")
    n_rows = length(quantities)
    n_cols = length(snap_list)
    x_size = 1700
    y_size = (x_size / n_cols) * n_rows - 110

    paths = joinpath.(
        temp_folder,
        vcat(
            ["stellar_mass_$(snap).jld2" for snap in snap_list],
            ["molecular_mass_$(snap).jld2" for snap in snap_list],
            ["atomic_mass_$(snap).jld2" for snap in snap_list],
            ["ionized_mass_$(snap).jld2" for snap in snap_list],
        ),
    )

    for slice in snap_list

        for (quantity, type, colorrange) in zip(quantities, types, colorbar_ranges)

            filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
                :all_subhalo,
                GalaxyInspector.plotParams(quantity).request,
            )

            GalaxyInspector.plotSnapshot(
                [Au6_MOL_path],
                request,
                [heatmap!];
                pf_kwargs=[(; colorrange)],
                output_path=temp_folder,
                base_filename="$(quantity)_$(slice)",
                slice,
                filter_function,
                da_functions=[GalaxyInspector.daDensity2DProjection],
                da_args=[(grid, quantity, type)],
                transform_box=true,
                translation,
                rotation,
                x_unit=u"kpc",
                y_unit=u"kpc",
                save_figures=false,
                backup_results=true,
            )

        end

    end

    with_theme(merge(theme_latexfonts(), GalaxyInspector.DEFAULT_THEME)) do

        f = Figure(size=(x_size, y_size), figure_padding=(5, 0, 0, 0))

        for (idx, path) in enumerate(paths)

            row = ceil(Int, idx / n_cols)
            col = mod1(idx, n_cols)

            xaxis_v = row == n_rows
            yaxis_v = col == 1

            filename = split(basename(path), '.')[1]
            snap_n = parse(Int, split(filename, '_')[3]) - 1

            ax = CairoMakie.Axis(
                f[row, col];
                xlabel=x_label,
                ylabel=y_label,
                xminorticksvisible=xaxis_v,
                xticksvisible=xaxis_v,
                xlabelvisible=xaxis_v,
                xticklabelsvisible=xaxis_v,
                yminorticksvisible=yaxis_v,
                yticksvisible=yaxis_v,
                ylabelvisible=yaxis_v,
                yticklabelsvisible=yaxis_v,
                xticklabelsize=28,
                yticklabelsize=28,
                xticks=[-20, -10, 0, 10, 20],
                yticks=[-20, -10, 0, 10, 20],
            )

            jldopen(path, "r") do jld2_file

                x, y, z = jld2_file["$(filename)_snap_$(lpad(snap_n, 3, '0'))"]["simulation_001"]

                pf = heatmap!(ax, x, y, z; colorrange=colorbar_ranges[row])

                if col == n_cols

                    Colorbar(
                        f[row, col+1],
                        pf,
                        label=colorbar_labels[row],
                        ticklabelsize=25,
                        ticks=colorbar_ticks[row],
                    )

                end

            end

            if row == 1

                time = times[col]

                if time < 1.0
                    time_stamp = round(time, digits=1)
                else
                    time_stamp = round(time, sigdigits=1)
                end

                text!(
                    ax,
                    0.0,
                    1.0;
                    text=L"t = %$(rpad(time_stamp, 3, '0')) \, \mathrm{Gyr}",
                    align=(:left, :top),
                    offset=(4, -2),
                    color=:white,
                    space=:relative,
                    fontsize=35,
                )

            end

            rowgap!(f.layout, 35)
            colgap!(f.layout, 35)

            rowsize!(f.layout, row, Makie.Fixed(pixelarea(ax.scene)[].widths[2]))

        end

        Makie.save(
            joinpath(
                figures_path,
                Au6_MOL_label,
                "density_maps_evolution_grid.png",
            ),
            f,
        )

    end

    rm(temp_folder; recursive=true)

    #############################################################
    # Evolution of the total mass and fractions of the different
    # gas components (within a sphere of radius r1)
    #############################################################

    r1_label = string(round(Int, ustrip(u"kpc", r1))) * "kpc"

    quantities = [:ionized_fraction, :atomic_fraction, :molecular_fraction]
    sim_labels = ["Ionized fraction", "Atomic fraction", "Molecular fraction"]

    x_plot_params = GalaxyInspector.plotParams(:physical_time)
    y_plot_params = GalaxyInspector.plotParams(:generic_fraction)

    temp_folder = joinpath(figures_path, Au6_MOL_label, "_gas_evolution")

    plotTimeSeries(
        fill(Au6_MOL_path, length(quantities)),
        [lines!];
        output_path=temp_folder,
        filename="gas_fractions_vs_physical_time_inside_$(r1_label)",
        output_format=".pdf",
        da_functions=[GalaxyInspector.daEvolution],
        da_args=[(:physical_time, quantity) for quantity in quantities],
        da_kwargs=[
            (;
                filter_mode=:subhalo,
                extra_filter=dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r1), (1, 1)),
                scaling=identity,
            ),
        ],
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=y_plot_params.exp_factor,
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        save_figures=false,
        backup_results=true,
    )

    quantities = [:stellar_mass, :hydrogen_mass, :ionized_mass, :atomic_mass, :molecular_mass]
    sim_labels = [
        "Stellar mass",
        "Hydrogen mass",
        "Ionized mass",
        "Atomic mass",
        "Molecular mass",
    ]

    y_plot_params = GalaxyInspector.plotParams(:generic_mass)

    plotTimeSeries(
        fill(Au6_MOL_path, length(quantities)),
        [lines!];
        output_path=temp_folder,
        filename="gas_masses_vs_physical_time_inside_$(r1_label)",
        output_format=".pdf",
        da_functions=[GalaxyInspector.daEvolution],
        da_args=[(:physical_time, quantity) for quantity in quantities],
        da_kwargs=[
            (;
                filter_mode=:subhalo,
                extra_filter=dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r1), (1, 1)),
                scaling=identity,
            ),
        ],
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=y_plot_params.exp_factor,
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        save_figures=false,
        backup_results=true,
    )

    paths = joinpath.(
        temp_folder,
        [
            "gas_masses_vs_physical_time_inside_40kpc.jld2",
            "gas_fractions_vs_physical_time_inside_40kpc.jld2",
        ],
    )

    # Starts at ~197 Myr to ignore initial very low fractions
    initial_snap = 20

    with_theme(merge(theme_latexfonts(), GalaxyInspector.DEFAULT_THEME)) do

        f = Figure(size=(880, 1200),)

        ax_1 = CairoMakie.Axis(
            f[1, 1];
            xlabel=L"t \, [\mathrm{Gyr}]",
            ylabel=L"M \, [\mathrm{10^{10} \, M_\odot}]",
            yscale=log10,
            xminorticksvisible=false,
            xticksvisible=false,
            xlabelvisible=false,
            xticklabelsvisible=false,
            limits=(nothing, nothing, 1.0e-5, nothing),
        )

        jldopen(paths[1], "r") do jld2_file

            first_address = first(keys(jld2_file))

            x_s, y_s = jld2_file[first_address]["simulation_001"]
            x_h, y_h = jld2_file[first_address]["simulation_002"]
            x_i, y_i = jld2_file[first_address]["simulation_003"]
            x_a, y_a = jld2_file[first_address]["simulation_004"]
            x_m, y_m = jld2_file[first_address]["simulation_005"]

            lines!(
                ax_1,
                x_s[initial_snap:end], y_s[initial_snap:end];
                color=Makie.wong_colors()[2],
                label="Stellar mass",
            )
            lines!(
                ax_1,
                x_h[initial_snap:end], y_h[initial_snap:end];
                color=:black,
                label="Hydrogen mass",
            )
            lines!(
                ax_1,
                x_i[initial_snap:end], y_i[initial_snap:end];
                color=Makie.wong_colors()[1],
                label="Ionized mass",
            )
            lines!(
                ax_1,
                x_a[initial_snap:end], y_a[initial_snap:end];
                color=Makie.wong_colors()[4],
                label="Atomic mass",
            )
            lines!(
                ax_1,
                x_m[initial_snap:end], y_m[initial_snap:end];
                color=Makie.wong_colors()[3],
                label="Molecular mass",
            )

        end

        axislegend(ax_1, position=:rb, framevisible=false, fontsize=20, nbanks=1)

        ax_2 = CairoMakie.Axis(
            f[2, 1];
            xlabel=L"t \, [\mathrm{Gyr}]",
            ylabel=L"\log_{10} \, f",
            aspect=nothing,
        )

        jldopen(paths[2], "r") do jld2_file

            first_address = first(keys(jld2_file))

            x_i, y_i = jld2_file[first_address]["simulation_001"]
            x_a, y_a = jld2_file[first_address]["simulation_002"]
            x_m, y_m = jld2_file[first_address]["simulation_003"]

            lines!(ax_2, x_i[initial_snap:end], log10.(y_i[initial_snap:end]); color=Makie.wong_colors()[1])
            lines!(ax_2, x_a[initial_snap:end], log10.(y_a[initial_snap:end]); color=Makie.wong_colors()[4])
            lines!(ax_2, x_m[initial_snap:end], log10.(y_m[initial_snap:end]); color=Makie.wong_colors()[3])

        end

        linkxaxes!(ax_1, ax_2)
        rowsize!(f.layout, 1, Relative(2 / 3))
        colsize!(f.layout, 1, Makie.Fixed(pixelarea(ax_1.scene)[].widths[2]))

        Makie.save(joinpath(figures_path, Au6_MOL_label, "gas_evolution_inside_40kpc.pdf"), f)

    end

    ##################################################################################
    # INFO FILE - Mean values for the gas fractions (from t_limit up to t = 13.8 Gyr)
    ##################################################################################

    t_limit = 3.0 # Gyr

    jldopen(paths[2], "r") do jld2_file

        first_address = first(keys(jld2_file))

        x_i, y_i = jld2_file[first_address]["simulation_001"]
        x_a, y_a = jld2_file[first_address]["simulation_002"]
        x_m, y_m = jld2_file[first_address]["simulation_003"]

        # Ionized fraction
        deleteat!(y_i, map(x -> x < t_limit, x_i))
        mean_fi = round(mean(y_i) * 100; sigdigits=3)
        std_fi = round(std(y_i) * 100; sigdigits=1)

        # Atomic fraction
        deleteat!(y_a, map(x -> x < t_limit, x_a))
        mean_fa = round(mean(y_a) * 100; sigdigits=3)
        std_fa = round(std(y_a) * 100; sigdigits=1)

        # Molecular fraction
        deleteat!(y_m, map(x -> x < t_limit, x_m))
        mean_fm = round(mean(y_m) * 100; sigdigits=3)
        std_fm = round(std(y_m) * 100; sigdigits=1)

        println(info_file, "Mean mass fractions inside 40 kpc (t = 3 Gyr - 13.8 Gyr):\n")
        println(info_file, "\tIonized fraction:   $(mean_fi)% ± $(std_fi)%")
        println(info_file, "\tAtomic fraction:    $(mean_fa)% ± $(std_fa)%")
        println(info_file, "\tMolecular fraction: $(mean_fm)% ± $(std_fm)%\n")
        println(info_file, "#"^100, "\n")

    end

    rm(temp_folder; force=true, recursive=true)

    ################################################################################################
    # For the cosmological simulations
    ################################################################################################

    ######################################################################################
    # Evolution of the total mass of the different gas components, SFR, and stellar mass
    # (within a sphere of radius r1)
    ######################################################################################

    x_plot_params = GalaxyInspector.plotParams(:physical_time)

    temp_folder = joinpath(figures_path, "comparison/_mass_evolution")

    simulations = [Au6_MOL_path, Au6_BLT_path, Au6_STD_path]
    sim_names = [Au6_MOL_label, Au6_BLT_label, Au6_STD_label]
    quantities = [
        [:stellar_mass, :ionized_mass, :neutral_mass, :molecular_mass],
        [:stellar_mass, :ionized_mass, :neutral_mass, :molecular_mass],
        [:stellar_mass, :ionized_mass, :neutral_mass],
    ]

    # Starts at ~197 Myr to ignore initial very low fractions
    initial_snap = 20

    n_panels = length(quantities[1]) + 1

    current_theme = merge(theme_latexfonts(), GalaxyInspector.DEFAULT_THEME)
    colors = current_theme[:palette][:color][]
    ylls = [exp10(-5.0), exp10(-2.5), exp10(-2.5), exp10(-5.0)]
    yhls = [exp10(1.5), nothing, nothing, nothing]

    for (i, simulation) in pairs(simulations)

        for quantity in quantities[i]

            y_plot_params = GalaxyInspector.plotParams(quantity)

            plotTimeSeries(
                [simulation],
                [lines!];
                output_path=temp_folder,
                filename="$(quantity)_$(basename(simulation))",
                output_format=".pdf",
                da_functions=[GalaxyInspector.daEvolution],
                da_args=[(:physical_time, quantity)],
                da_kwargs=[
                    (;
                        filter_mode=:subhalo,
                        extra_filter=dd -> GalaxyInspector.filterWithinSphere(
                            dd,
                            (0.0u"kpc", r1),
                            (1, 1),
                        ),
                        scaling=identity,
                    ),
                ],
                x_unit=x_plot_params.unit,
                y_unit=y_plot_params.unit,
                x_exp_factor=x_plot_params.exp_factor,
                y_exp_factor=y_plot_params.exp_factor,
                save_figures=false,
                backup_results=true,
            )

        end

        filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
            :subhalo,
            GalaxyInspector.plotParams(:sfr).request,
        )

        plotSnapshot(
            [simulation],
            request,
            [lines!];
            output_path=temp_folder,
            base_filename="sfr_$(basename(simulation))",
            slice=128,
            filter_function,
            da_functions=[GalaxyInspector.daStellarHistory],
            da_args=[()],
            da_kwargs=[
                (;
                    quantity=:sfr,
                    n_bins=80,
                    filter_function=dd -> GalaxyInspector.filterWithinSphere(
                        dd,
                        (0.0u"kpc", r1),
                        :zero,
                    ),
                ),
            ],
            transform_box=true,
            translation,
            rotation,
            x_unit=x_plot_params.unit,
            y_unit=GalaxyInspector.plotParams(:sfr).unit,
            save_figures=false,
            backup_results=true,
        )

    end

    jld2_paths = joinpath.(
        temp_folder,
        vcat(
            ["sfr_$(label).jld2" for label in sim_names],
            ["stellar_mass_$(label).jld2" for label in sim_names],
            ["ionized_mass_$(label).jld2" for label in sim_names],
            ["neutral_mass_$(label).jld2" for label in sim_names],
            ["molecular_mass_$(label).jld2" for label in sim_names[1:2]],
        ),
    )

    xlabel = LaTeXString(
        replace(
            x_plot_params.axis_label,
            "auto_label" => GalaxyInspector.getLabel(
                x_plot_params.var_name,
                x_plot_params.exp_factor,
                x_plot_params.unit,
            ),
        ),
    )

    with_theme(current_theme) do

        f = Figure(size=(880, 2300), figure_padding=(1, 15, 5, 20))

        ylabel = LaTeXString(
            replace(
                GalaxyInspector.plotParams(:sfr).axis_label,
                "auto_label" => GalaxyInspector.getLabel(
                    GalaxyInspector.plotParams(:sfr).var_name,
                    GalaxyInspector.plotParams(:sfr).exp_factor,
                    GalaxyInspector.plotParams(:sfr).unit,
                ),
            ),
        )

        ax = CairoMakie.Axis(
            f[1, 1];
            xlabel,
            ylabel,
            aspect=AxisAspect(1.7),
            yscale=log10,
            xlabelvisible=false,
            xticklabelsvisible=false,
            xticks=0:2:14,
            yticks=[exp10(-2.0), exp10(-1.0), exp10(0.0), exp10(1.0)],
            limits=(nothing, nothing, exp10(-2.5), exp10(1.5)),
        )

        for (path, label, color) in zip(jld2_paths[1:3], sim_names, colors)

            jldopen(path, "r") do jld2_file

                address = first(keys(jld2_file))

                x, y = jld2_file[address]["simulation_001"]

                lines!(ax, x, y; label, linestyle=:solid, color)

            end

        end

        axislegend(ax, position=:rb, framevisible=false, fontsize=20, nbanks=1)

        for (i, (quantity, yll, yhl)) in enumerate(zip(quantities[1], ylls, yhls))

            y_plot_params = GalaxyInspector.plotParams(quantity)

            ylabel = LaTeXString(
                replace(
                    y_plot_params.axis_label,
                    "auto_label" => GalaxyInspector.getLabel(
                        y_plot_params.var_name,
                        y_plot_params.exp_factor,
                        y_plot_params.unit,
                    ),
                ),
            )

            ax = CairoMakie.Axis(
                f[i+1, 1];
                xlabel,
                ylabel,
                aspect=AxisAspect(1.7),
                yscale=log10,
                xlabelvisible=i+1 == n_panels,
                xticklabelsvisible=i+1 == n_panels,
                xticks=0:2:14,
                limits=(nothing, nothing, yll, yhl),
            )

            paths = jld2_paths[(3*i+1):min(3*i+3, length(jld2_paths))]

            for (path, label, color) in zip(paths, sim_names, colors)

                jldopen(path, "r") do jld2_file

                    address = first(keys(jld2_file))

                    x, y = jld2_file[address]["simulation_001"]

                    lines!(
                        ax,
                        x[initial_snap:end],
                        y[initial_snap:end];
                        label,
                        linestyle=:solid,
                        color,
                    )

                end

            end

        end

        Makie.save(joinpath(figures_path, "comparison/mass_evolution.pdf"), f)

    end

    rm(temp_folder; recursive=true)

    #################################################################
    # Gas components density map (face-on projections at redshift 0)
    #################################################################

    temp_folder = joinpath(figures_path, "comparison/_density_maps")

    grid = GalaxyInspector.CubicGrid(65.0u"kpc", 400)
    types = [:particles, :cells, :cells, :cells]
    quantities = [
        [:stellar_mass, :ionized_mass, :neutral_mass, :molecular_mass],
        [:stellar_mass, :ionized_mass, :neutral_mass, :molecular_mass],
        [:stellar_mass, :ionized_mass, :neutral_mass]
    ]
    sim_names = [Au6_MOL_label, Au6_BLT_label, Au6_STD_label]

    colorbar_labels = [
        L"\log_{10} \, \Sigma_\star \,\, [\mathrm{M_\odot \, kpc^{-2}}]"
        L"\log_{10} \, \Sigma_\mathrm{HII} \,\, [\mathrm{M_\odot \, kpc^{-2}}]"
        L"\log_{10} \, \Sigma_\mathrm{HI + H_2} \,\, [\mathrm{M_\odot \, kpc^{-2}}]"
        L"\log_{10} \, \Sigma_\mathrm{H_2} \,\, [\mathrm{M_\odot \, kpc^{-2}}]"
    ]

    x_label = GalaxyInspector.getLabel("x", 0, u"kpc")
    y_label = GalaxyInspector.getLabel("y", 0, u"kpc")
    z_label = GalaxyInspector.getLabel("z", 0, u"kpc")
    n_rows = length(sim_names)
    n_cols = length(quantities[1])
    x_size = 1700
    y_size = (x_size / n_cols) * n_rows + 150
    colorranges = [(5.5, 10.5), (4.5, 8.5), (4.5, 8.5), (2.5, 8.5)]
    cticks = [6:1:10, 5:1:8, 5:1:8, 3:1:8]

    paths = joinpath.(
        temp_folder,
        vcat([["$(qty)_$(sim).jld2" for qty in quantities[i]] for (i, sim) in pairs(sim_names)]...),
    )

    for (i, (simulation, label)) in enumerate(zip(simulation_paths, sim_names))

        for (quantity, type) in zip(quantities[i], types)

            filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
                :all_subhalo,
                GalaxyInspector.plotParams(quantity).request,
            )

            GalaxyInspector.plotSnapshot(
                [simulation],
                request,
                [heatmap!];
                output_path=temp_folder,
                base_filename="$(quantity)_$(label)",
                slice=128,
                filter_function,
                da_functions=[GalaxyInspector.daDensity2DProjection],
                da_args=[(grid, quantity, type)],
                transform_box=true,
                translation,
                rotation,
                x_unit=u"kpc",
                y_unit=u"kpc",
                save_figures=false,
                backup_results=true,
            )

        end

    end

    with_theme(merge(theme_latexfonts(), GalaxyInspector.DEFAULT_THEME)) do

        f = Figure(size=(x_size, y_size), figure_padding=(5, 15, 0, 0))

        for (idx, path) in enumerate(paths)

            row = ceil(Int, idx / n_cols)
            col = mod1(idx, n_cols)

            xaxis_v = row == n_rows
            yaxis_v = col == 1

            ax = CairoMakie.Axis(
                f[row+1, col];
                xlabel=x_label,
                ylabel=(row == 1 ? y_label : z_label),
                xminorticksvisible=xaxis_v,
                xticksvisible=xaxis_v,
                xlabelvisible=xaxis_v,
                xticklabelsvisible=xaxis_v,
                yminorticksvisible=yaxis_v,
                yticksvisible=yaxis_v,
                ylabelvisible=yaxis_v,
                yticklabelsvisible=yaxis_v,
                xticklabelsize=28,
                yticklabelsize=28,
                xticks=[-30, -20, -10, 0, 10, 20, 30],
                yticks=[-30, -20, -10, 0, 10, 20, 30],
            )

            jldopen(path, "r") do jld2_file

                x, y, z = jld2_file[first(keys(jld2_file))]["simulation_001"]

                pf = heatmap!(ax, x, y, z; colorrange=colorranges[col])

                if col == 1

                    GalaxyInspector.ppAnnotation!(
                        f,
                        sim_names[row];
                        position=(0.04, 0.98),
                        color=:white,
                        fontsize=30,
                    )

                end

                if row == 1

                    Colorbar(
                        f[row, col],
                        pf,
                        label=colorbar_labels[col],
                        ticklabelsize=23,
                        ticks=cticks[col],
                        vertical=false,
                    )

                end

            end

            colgap!(f.layout, 20)

        end

        Makie.save(joinpath(figures_path, "comparison/gas_density_maps.png"), f)

    end

    rm(temp_folder; recursive=true)

    ################################################################################################
    # Resolved Kennicutt–Schmidt law - Scatter
    # Redshift 0 - High resolution
    ################################################################################################

    #             11, 11.5, 12, 12.5, 13, 13.5, 14 Gyr
    snap_list = [111, 114, 117, 120, 123, 126, 128]

    ####################
    # Molecular density
    ####################

    temp_folder = mkpath(joinpath(figures_path, "comparison/_ks_law_molecular"))

    kennicuttSchmidtLaw(
        [Au6_MOL_path],
        snap_list;
        quantity=:molecular_mass,
        grid_size=30.0u"kpc",
        gas_weights=nothing,
        measurements=true,
        output_file=joinpath(temp_folder, "Au6_MOL.png"),
        filter_mode=:subhalo,
        sim_labels=[Au6_MOL_label],
        theme=Theme(
            size=(850, 780),
            Axis=(
                limits=(4.0, 10.0, -4.3, 1.3),
                xticks=4:1:10,
                yticks=-4:1:1,
            ),
            Legend=(
                nbanks=1,
                labelsize=30,
                halign=:left,
                valign=:top,
                padding=(60, 0, 0, 0),
            ),
        ),
    )

    kennicuttSchmidtLaw(
        [Au6_BLT_path],
        snap_list;
        quantity=:molecular_mass,
        grid_size=30.0u"kpc",
        gas_weights=nothing,
        measurements=true,
        output_file=joinpath(temp_folder, "Au6_BLT.png"),
        filter_mode=:subhalo,
        sim_labels=[Au6_BLT_label],
        theme=Theme(
            size=(850, 780),
            Axis=(
                limits=(4.0, 10.0, -4.3, 1.3),
                xticks=4:1:10,
                yticks=-4:1:1,
            ),
            Legend=(
                nbanks=1,
                labelsize=30,
                halign=:left,
                valign=:top,
                padding=(60, 0, 0, 0),
            ),
        ),
    )

    GalaxyInspector.hvcatImages(
        2,
        joinpath.(temp_folder, ["Au6_MOL.png", "Au6_BLT.png"]);
        output_path=joinpath(figures_path, "comparison/molecular_bigiel2008.png"),
    )

    rm(temp_folder; recursive=true)

    ##################
    # Neutral density
    ##################

    temp_folder = mkpath(joinpath(figures_path, "comparison/_ks_law_neutral"))

    kennicuttSchmidtLaw(
        [Au6_MOL_path],
        snap_list;
        quantity=:neutral_mass,
        grid_size=30.0u"kpc",
        gas_weights=nothing,
        measurements=true,
        output_file=joinpath(temp_folder, "Au6_MOL.png"),
        filter_mode=:subhalo,
        sim_labels=[Au6_MOL_label],
        theme=Theme(
            size=(850, 780),
            Axis=(
                limits=(5.5, 10.0, -4.3, 1.3),
                xticks=6:1:10,
                yticks=-4:1:1,
            ),
            Legend=(
                nbanks=1,
                labelsize=30,
                halign=:left,
                valign=:top,
                padding=(60, 0, 0, 0),
            ),
        ),
    )

    kennicuttSchmidtLaw(
        [Au6_BLT_path],
        snap_list;
        quantity=:neutral_mass,
        grid_size=30.0u"kpc",
        gas_weights=nothing,
        measurements=true,
        output_file=joinpath(temp_folder, "Au6_BL.png"),
        filter_mode=:subhalo,
        sim_labels=[Au6_BLT_label],
        theme=Theme(
            size=(850, 780),
            Axis=(
                limits=(5.5, 10.0, -4.3, 1.3),
                xticks=6:1:10,
                yticks=-4:1:1,
            ),
            Legend=(
                nbanks=1,
                labelsize=30,
                halign=:left,
                valign=:top,
                padding=(60, 0, 0, 0),
            ),
        ),
    )

    kennicuttSchmidtLaw(
        [Au6_STD_path],
        snap_list;
        quantity=:neutral_mass,
        grid_size=30.0u"kpc",
        gas_weights=nothing,
        measurements=true,
        output_file=joinpath(temp_folder, "Au6_STD.png"),
        filter_mode=:subhalo,
        sim_labels=[Au6_STD_label],
        theme=Theme(
            size=(850, 780),
            Axis=(
                limits=(4.5, 10.5, -4.1, 1.1),
                xticks=5:1:10,
                yticks=-4:1:1,
            ),
            Legend=(
                nbanks=1,
                labelsize=30,
                halign=:left,
                valign=:top,
                padding=(60, 0, 0, 0),
            ),
        ),
    )

    GalaxyInspector.hvcatImages(
        3,
        joinpath.(temp_folder, ["Au6_MOL.png", "Au6_BL.png", "Au6_STD.png"]);
        output_path=joinpath(figures_path, "comparison/neutral_bigiel2008.png"),
    )

    rm(temp_folder; recursive=true)

    ################################################################################################
    # Comparison between different resolutions
    ################################################################################################

    ff = dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r1), (1, 1))

    temp_folder = joinpath(figures_path, "comparison/_resolution")

    timeSeries(
        [Au6_MOL_path, Au6_MOL_LR_path],
        :physical_time,
        :sfr;
        output_path=temp_folder,
        filter_mode=:subhalo,
        extra_filter=ff,
        sim_labels=[Au6_MOL_label, Au6_MOL_LR_label],
        backup_results=true,
    )

    timeSeries(
        [Au6_MOL_path, Au6_MOL_LR_path],
        :physical_time,
        :stellar_mass;
        output_path=temp_folder,
        filter_mode=:subhalo,
        extra_filter=ff,
        sim_labels=nothing,
        backup_results=true,
    )

    timeSeries(
        [Au6_MOL_path, Au6_MOL_LR_path],
        :physical_time,
        :stellar_metallicity;
        output_path=temp_folder,
        filter_mode=:subhalo,
        extra_filter=ff,
        sim_labels=nothing,
        backup_results=true,
    )

    timeSeries(
        [Au6_MOL_path, Au6_MOL_LR_path],
        :physical_time,
        :ionized_mass;
        output_path=temp_folder,
        filter_mode=:subhalo,
        extra_filter=ff,
        sim_labels=nothing,
        backup_results=true,
    )

    timeSeries(
        [Au6_MOL_path, Au6_MOL_LR_path],
        :physical_time,
        :atomic_mass;
        output_path=temp_folder,
        filter_mode=:subhalo,
        extra_filter=ff,
        sim_labels=nothing,
        backup_results=true,
    )

    timeSeries(
        [Au6_MOL_path, Au6_MOL_LR_path],
        :physical_time,
        :molecular_mass;
        output_path=temp_folder,
        filter_mode=:subhalo,
        extra_filter=ff,
        sim_labels=nothing,
        backup_results=true,
    )

    jld2_names = [
        "sfr_vs_physical_time",
        "stellar_mass_vs_physical_time",
        "stellar_metallicity_vs_physical_time",
        "ionized_mass_vs_physical_time",
        "atomic_mass_vs_physical_time",
        "molecular_mass_vs_physical_time",
    ]

    y_labels = [
        L"\mathrm{SFR \, [M_\odot \, yr^{-1}]}",
        L"M_\star \, \mathrm{[10^{10} \, M_\odot]}",
        L"Z_\star",
        L"M_\mathrm{HII} \, \mathrm{[10^{10} \, M_\odot]}",
        L"M_\mathrm{HI} \, \mathrm{[10^{10} \, M_\odot]}",
        L"M_\mathrm{H_2} \, \mathrm{[10^{10} \, M_\odot]}",
    ]

    n_cols = 3
    t_limit = 3.0

    with_theme(merge(theme_latexfonts(), GalaxyInspector.DEFAULT_THEME, Theme())) do

        f = Figure(size=(1700, 1130),)

        for idx in eachindex(jld2_names)

            row = ceil(Int, idx / n_cols)
            col = mod1(idx, n_cols)

            ax = CairoMakie.Axis(
                f[row, col];
                xlabel=L"t \, [\mathrm{Gyr}]",
                ylabel=y_labels[idx],
                yscale=log10,
            )

            jldopen(joinpath(temp_folder, "$(jld2_names[idx]).jld2"), "r") do jld2_file

                hr_data = jld2_file["$(jld2_names[idx])/simulation_001"]
                lr_data = jld2_file["$(jld2_names[idx])/simulation_002"]

                hr_t_idxs = map(x -> x > t_limit, hr_data[1])
                lr_t_idxs = map(x -> x > t_limit, lr_data[1])

                if row == 2 && col == 1
                    label_hr = Au6_MOL_label
                    label_lr = Au6_MOL_LR_label
                else
                    label_hr = nothing
                    label_lr = nothing
                end

                lines!(
                    ax,
                    lr_data[1][lr_t_idxs],
                    lr_data[2][lr_t_idxs];
                    color=Makie.wong_colors()[1],
                    label=label_lr,
                )

                lines!(
                    ax,
                    hr_data[1][hr_t_idxs],
                    hr_data[2][hr_t_idxs];
                    color=Makie.wong_colors()[2],
                    label=label_hr,
                )

                axislegend(ax, position=:rb, framevisible=false, fontsize=20, nbanks=1)

            end

        end

        Makie.save(joinpath(figures_path, "comparison/resolution_comparison.png"), f)

    end

    ################################
    # INFO FILE - Resolution ratios
    ################################

    paths = joinpath.(temp_folder, jld2_names .* ".jld2")

    component_labels = ["SFR", "Stellar mass", "Stellar metallicity", "HII", "HI", "H2"]

    println(
        info_file,
        "The maximum porcentual differences between Au6_MOL and Au6_MOL_LR are after 10Gyr:\n",
    )

    for (path, label) in zip(paths, component_labels)

        jldopen(path, "r") do jld2_file

            x_Au6_MOL, y_Au6_MOL = jld2_file[first(keys(jld2_file))]["simulation_001"]
            x_Au6_MOL_LR, y_Au6_MOL_LR = jld2_file[first(keys(jld2_file))]["simulation_002"]

            idxs = map(t -> t < 10.0, x_Au6_MOL)
            idxs_lr = map(t -> t < 10.0, x_Au6_MOL_LR)

            deleteat!(y_Au6_MOL, idxs)
            deleteat!(y_Au6_MOL_LR, idxs_lr)

            deltas = @. abs(y_Au6_MOL - y_Au6_MOL_LR) / y_Au6_MOL

            max_diff = round(maximum(deltas) * 100; sigdigits=3)

            println(info_file, "$(label): $(max_diff)%\n")

        end

    end

    println(info_file, "#"^100, "\n")

    rm(temp_folder; recursive=true)

    ################################################################################################
    # Report files
    ################################################################################################

    simulationReport(simulation_paths; output_path=report_path)

    snapshotReport(simulation_paths, [128]; output_path=report_path, filter_mode=:subhalo)

    ##############################################
    # INFO FILE - Insitu/Exsitu stellar fractions
    ##############################################

    filter_function, _, _, request = GalaxyInspector.selectFilter(
        :subhalo,
        Dict(:stars => ["ID  "]),
    )

    data_dict = GalaxyInspector.makeDataDict(Au6_MOL_path, Au6_MOL_slice, request)

    GalaxyInspector.filterData!(data_dict; filter_function)

    insitu = GalaxyInspector.filterData(
        data_dict;
        filter_function=dd->GalaxyInspector.filterByBirthPlace(dd, :exsitu),
    )

    exsitu = GalaxyInspector.filterData(
        data_dict;
        filter_function=dd->GalaxyInspector.filterByBirthPlace(dd, :insitu),
    )

    insitu_mass = sum(insitu[:stars]["MASS"]; init=0.0u"Msun")
    exsitu_mass = sum(exsitu[:stars]["MASS"]; init=0.0u"Msun")

    total_mass = insitu_mass + exsitu_mass

    println(
        info_file,
        "Stellar insitu fraction: $(round((insitu_mass / total_mass) * 100; sigdigits=3))%",
    )
    println(
        info_file,
        "Stellar exsitu fraction: $(round((exsitu_mass / total_mass) * 100; sigdigits=3))%\n",
    )
    println(info_file, "#"^100, "\n")

    #####################################
    # INFO FILE - Mean SFR for t > 3 Gyr
    #####################################

    filter_function, _, _, request = GalaxyInspector.selectFilter(
        :subhalo,
        Dict(:stars => ["GAGE", "MASS"]),
    )

    data_dict = makeDataDict(Au6_MOL_path, Au6_MOL_slice, request)

    ts, sfrs = GalaxyInspector.daStellarHistory(data_dict; filter_function)

    idxs = map(x -> x > 3.0u"Gyr", ts)

    sfr = round(ustrip.(u"Msun*yr^-1", mean(sfrs[idxs])); digits=2)

    println(info_file, "Mean SFR for t > 3 Gyr: $(sfr) M⊙ yr^-1\n")

    println(info_file, "#"^100, "\n")

    #######################################################################
    # INFO FILE - Change in R95 when going from Rmax = 40 to Rmax = 60 kpc
    #######################################################################

    filter_function, translation, _, request = GalaxyInspector.selectFilter(
        :subhalo,
        GalaxyInspector.mergeRequests(
            GalaxyInspector.plotParams(:gas_mass).request,
            GalaxyInspector.plotParams(:stellar_mass).request,
            GalaxyInspector.plotParams(:molecular_mass).request,
            GalaxyInspector.plotParams(:ionized_mass).request,
            GalaxyInspector.plotParams(:atomic_mass).request,
        ),
    )

    radial_limits = [40.0, 60.0] .* u"kpc"

    data_dict = makeDataDict(Au6_MOL_path, Au6_MOL_slice, request)

    GalaxyInspector.translateData!(data_dict, translation)

    GalaxyInspector.filterData!(data_dict; filter_function)

    stellar_masses = data_dict[:stars]["MASS"]
    gas_masses = data_dict[:gas]["MASS"]
    ionized_masses = GalaxyInspector.computeMass(data_dict, :ionized)
    atomic_masses = GalaxyInspector.computeMass(data_dict, :atomic)
    molecular_masses = GalaxyInspector.computeMass(data_dict, :molecular)

    stellar_R95 = Vector{Unitful.Length}(undef, length(radial_limits))
    gas_R95 = Vector{Unitful.Length}(undef, length(radial_limits))
    ionized_R95 = Vector{Unitful.Length}(undef, length(radial_limits))
    atomic_R95 = Vector{Unitful.Length}(undef, length(radial_limits))
    molecular_R95 = Vector{Unitful.Length}(undef, length(radial_limits))

    for (i, radial_limit) in enumerate(radial_limits)

        disc_idxs = GalaxyInspector.filterWithinSphere(data_dict, (0.0u"kpc", radial_limit), :zero)

        stellar_mass_inside = stellar_masses[disc_idxs[:stars]]
        gas_mass_inside = gas_masses[disc_idxs[:gas]]
        ionized_mass_inside = ionized_masses[disc_idxs[:gas]]
        atomic_mass_inside = atomic_masses[disc_idxs[:gas]]
        molecular_mass_inside = molecular_masses[disc_idxs[:gas]]

        stellar_R95[i] = GalaxyInspector.computeMassRadius(
            data_dict[:stars]["POS "][:, disc_idxs[:stars]],
            stellar_mass_inside;
            percent=95.0,
        )

        gas_R95[i] = GalaxyInspector.computeMassRadius(
            data_dict[:gas]["POS "][:, disc_idxs[:gas]],
            gas_mass_inside;
            percent=95.0,
        )

        ionized_R95[i] = GalaxyInspector.computeMassRadius(
            data_dict[:gas]["POS "][:, disc_idxs[:gas]],
            ionized_mass_inside;
            percent=95.0,
        )

        atomic_R95[i] = GalaxyInspector.computeMassRadius(
            data_dict[:gas]["POS "][:, disc_idxs[:gas]],
            atomic_mass_inside;
            percent=95.0,
        )

        molecular_R95[i] = GalaxyInspector.computeMassRadius(
            data_dict[:gas]["POS "][:, disc_idxs[:gas]],
            molecular_mass_inside;
            percent=95.0,
        )

    end

    stella_change = uconvert(
        Unitful.NoUnits,
        (stellar_R95[2] - stellar_R95[1]) / stellar_R95[1],
    ) * 100

    gas_change = uconvert(
        Unitful.NoUnits,
        (gas_R95[2] - gas_R95[1]) / gas_R95[1],
    ) * 100

    ionized_change = uconvert(
        Unitful.NoUnits,
        (ionized_R95[2] - ionized_R95[1]) / ionized_R95[1],
    ) * 100

    atomic_change = uconvert(
        Unitful.NoUnits,
        (atomic_R95[2] - atomic_R95[1]) / atomic_R95[1],
    ) * 100

    molecular_change = uconvert(
        Unitful.NoUnits,
        (molecular_R95[2] - molecular_R95[1]) / molecular_R95[1],
    ) * 100

    println(
        info_file,
        "Stellar R95 change   (40 kpc <-> 60 kpc): $(round(stella_change; digits=1))%",
    )
    println(
        info_file,
        "Gas R95 change       (40 kpc <-> 60 kpc): $(round(gas_change; digits=1))%",
    )
    println(
        info_file,
        "Ionized R95 change   (40 kpc <-> 60 kpc): $(round(ionized_change; digits=1))%",
    )
    println(
        info_file,
        "Atomic R95 change    (40 kpc <-> 60 kpc): $(round(atomic_change; digits=1))%",
    )
    println(
        info_file,
        "Molecular R95 change (40 kpc <-> 60 kpc): $(round(molecular_change; digits=1))%\n",
    )
    println(info_file, "#"^100, "\n")

    #############################################
    # INFO FILE - Atomic/ionized gas mass ratios
    #############################################

    filter_function, translation, _, request = GalaxyInspector.selectFilter(
        :subhalo,
        GalaxyInspector.mergeRequests(
            GalaxyInspector.plotParams(:ionized_mass).request,
            GalaxyInspector.plotParams(:atomic_mass).request,
        ),
    )

    data_dict = makeDataDict(Au6_MOL_path, Au6_MOL_slice, request)

    GalaxyInspector.translateData!(data_dict, translation)

    GalaxyInspector.filterData!(data_dict; filter_function)

    idx_within = GalaxyInspector.filterWithinSphere(
        data_dict,
        (0.0u"kpc", 60.0u"kpc"),
        :zero,
    )[:gas]

    idx_outside = GalaxyInspector.filterWithinSphere(
        data_dict,
        (60.0u"kpc", 200.0u"kpc"),
        :zero,
    )[:gas]

    ionized_masses = GalaxyInspector.computeMass(data_dict, :ionized)
    atomic_masses = GalaxyInspector.computeMass(data_dict, :atomic)
    gas_masses = data_dict[:gas]["MASS"]

    factor_within = uconvert(
        Unitful.NoUnits,
        sum(atomic_masses[idx_within]) / sum(ionized_masses[idx_within]),
    )

    println(
        info_file,
        "The atomic fraction dominates in the inner 60kpc over the ionized one, by a \
        factor of $(round(factor_within; sigdigits=3))\n",
    )

    factor_outside = uconvert(
        Unitful.NoUnits,
        sum(ionized_masses[idx_outside]) / sum(atomic_masses[idx_outside]),
    )

    println(
        info_file,
        "Outside the central region, the relative amounts of ionized and atomic material \
        remain approximately constant, with a ratio of ionized to atomic mass of \
        $(round(factor_outside; sigdigits=3))\n",
    )
    println(info_file, "#"^100, "\n")

    ##################################################
    # INFO FILE - Stellar mass ratios (at redshift 0)
    ##################################################

    filter_function, _, _, request = GalaxyInspector.selectFilter(
        :subhalo,
        GalaxyInspector.plotParams(:stellar_mass).request,
    )

    stellar_masses = Vector{Number}(undef, length(simulation_paths))

    for (i, path) in enumerate(simulation_paths)

        data_dict = makeDataDict(path, 128, request)

        GalaxyInspector.filterData!(data_dict; filter_function)

        stellar_masses[i] = sum(data_dict[:stars]["MASS"], init=0.0u"Msun")

    end

    factor_MOL = uconvert(
        Unitful.NoUnits,
        (stellar_masses[3] - stellar_masses[1]) / stellar_masses[3],
    ) * 100

    factor_BL = uconvert(
        Unitful.NoUnits,
        (stellar_masses[3] - stellar_masses[2]) / stellar_masses[3],
    ) * 100

    println(info_file, "Relative reduction in total stellar mass, at z = 0\n")
    println(info_file, "\tGoing from from Au6_STD to Au6_MOL: $(round(factor_MOL; sigdigits=3))%")
    println(info_file, "\tGoing from from Au6_STD to Au6_BL:  $(round(factor_BL; sigdigits=3))%")

    ################################################################################################
    # Mass table
    ################################################################################################

    simulations = [Au6_MOL_path, Au6_BLT_path, Au6_STD_path]
    row_labels = ["Stars", "Gas", "HII", L"\mathrm{HI + H_2}", L"\mathrm{H_2}"]
    header = [
        L"R_{95}",
        L"M_\mathrm{tot}",
        L"M_\mathrm{< 40 \, kpc}",
        L"R_{95}",
        L"M_\mathrm{tot}",
        L"M_\mathrm{< 40 \, kpc}",
        L"R_{95}",
        L"M_\mathrm{tot}",
        L"M_\mathrm{< 40 \, kpc}",
    ]
    units = [
        L"[\mathrm{kpc}]",
        L"[10^{10} \, \mathrm{M_\odot}]",
        L"[10^{10} \, \mathrm{M_\odot}]",
        L"[\mathrm{kpc}]",
        L"[10^{10} \, \mathrm{M_\odot}]",
        L"[10^{10} \, \mathrm{M_\odot}]",
        L"[\mathrm{kpc}]",
        L"[10^{10} \, \mathrm{M_\odot}]",
        L"[10^{10} \, \mathrm{M_\odot}]",
    ]

    values = Matrix{Float64}(undef, length(row_labels), 3 * length(simulations))

    filter_function, translation, _, request = GalaxyInspector.selectFilter(
        :subhalo,
        GalaxyInspector.mergeRequests(
            GalaxyInspector.plotParams(:stellar_mass).request,
            GalaxyInspector.plotParams(:gas_mass).request,
            GalaxyInspector.plotParams(:ionized_mass).request,
            GalaxyInspector.plotParams(:neutral_mass).request,
            GalaxyInspector.plotParams(:molecular_mass).request,
        )
    )

    massFormater(dd, qty) = round(
        ustrip(u"Msun", sum(GalaxyInspector.computeMass(dd, qty); init=0.0u"Msun") / exp10(10));
        sigdigits=3,
    )

    for (i, simulation) in pairs(simulations)

        data_dict = makeDataDict(simulation, 128, request)

        GalaxyInspector.filterData!(data_dict; filter_function)

        GalaxyInspector.translateData!(data_dict, translation)

        #################################
        # Mtot (within the main subhalo)
        #################################

        Mtot_stars = massFormater(data_dict, :stars)
        Mtot_gas = massFormater(data_dict, :gas)
        Mtot_hii = massFormater(data_dict, :ionized)
        Mtot_h2 = massFormater(data_dict, :molecular)
        Mtot_hn = massFormater(data_dict, :neutral)

        GalaxyInspector.filterData!(
            data_dict;
            filter_function=dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r1), :zero)
        )

        ###############
        # R95 (r < r1)
        ###############

        stellar_masses = GalaxyInspector.computeMass(data_dict, :stars)
        gas_masses = GalaxyInspector.computeMass(data_dict, :gas)
        ionized_masses = GalaxyInspector.computeMass(data_dict, :ionized)
        neutral_masses = GalaxyInspector.computeMass(data_dict, :neutral)
        molecular_masses = GalaxyInspector.computeMass(data_dict, :molecular)

        ########
        # Stars
        ########

        mass_radius_95 = GalaxyInspector.computeMassRadius(
            data_dict[:stars]["POS "],
            stellar_masses;
            percent=95.0,
        )

        R95_stars = round(ustrip(u"kpc", mass_radius_95); sigdigits=3)

        ############
        # Total gas
        ############

        mass_radius_95 = GalaxyInspector.computeMassRadius(
            data_dict[:gas]["POS "],
            gas_masses;
            percent=95.0,
        )

        R95_gas = round(ustrip(u"kpc", mass_radius_95); sigdigits=3)

        ##############
        # Ionized gas
        ##############

        mass_radius_95 = GalaxyInspector.computeMassRadius(
            data_dict[:gas]["POS "],
            ionized_masses;
            percent=95.0,
        )

        R95_hii = round(ustrip(u"kpc", mass_radius_95); sigdigits=3)

        ##############
        # Neutral gas
        ##############

        mass_radius_95 = GalaxyInspector.computeMassRadius(
            data_dict[:gas]["POS "],
            neutral_masses;
            percent=95.0,
        )

        R95_hn = round(ustrip(u"kpc", mass_radius_95); sigdigits=3)

        ################
        # Molecular gas
        ################

        mass_radius_95 = GalaxyInspector.computeMassRadius(
            data_dict[:gas]["POS "],
            molecular_masses;
            percent=95.0,
        )

        R95_h2 = round(ustrip(u"kpc", mass_radius_95); sigdigits=3)

        ################
        # Mtot (r < r1)
        ################

        M40_stars = massFormater(data_dict, :stars)
        M40_gas = massFormater(data_dict, :gas)
        M40_hii = massFormater(data_dict, :ionized)
        M40_hn = massFormater(data_dict, :neutral)
        M40_h2 = massFormater(data_dict, :molecular)

        ##############
        # Save values
        ##############

        values[1, 1+3*(i-1)] = R95_stars
        values[1, 2+3*(i-1)] = Mtot_stars
        values[1, 3+3*(i-1)] = M40_stars

        values[2, 1+3*(i-1)] = R95_gas
        values[2, 2+3*(i-1)] = Mtot_gas
        values[2, 3+3*(i-1)] = M40_gas

        values[3, 1+3*(i-1)] = R95_hii
        values[3, 2+3*(i-1)] = Mtot_hii
        values[3, 3+3*(i-1)] = M40_hii

        values[4, 1+3*(i-1)] = R95_hn
        values[4, 2+3*(i-1)] = Mtot_hn
        values[4, 3+3*(i-1)] = M40_hn

        if i == 3

            values[5, 1+3*(i-1)] = NaN
            values[5, 2+3*(i-1)] = NaN
            values[5, 3+3*(i-1)] = NaN

        else

            values[5, 1+3*(i-1)] = R95_h2
            values[5, 2+3*(i-1)] = Mtot_h2
            values[5, 3+3*(i-1)] = M40_h2

        end

    end

    pretty_table(
        table_file,
        replace(values, NaN => "--");
        vlines=:all,
        backend=Val(:latex),
        row_label_alignment=:c,
        alignment=:c,
        header=(header, units),
        row_labels,
    )

    ################################################################################################
    # Close files
    ################################################################################################

    close(info_file)
    close(table_file)

    if logging
        close(log_file)
    end

    return nothing

end

function (@main)(ARGS)

    LOGGING = true
    BASE_OUT_PATH = "./"
    BASE_SRC_PATH = "F:/simulations/2024-11-05"
    SIMULATIONS = ["Au6_MOL", "Au6_BLT", "Au6_STD", "Au6_MOL_LR"]
    LABELS = ["Au6_MOL", "Au6_BLT", "Au6_STD", "Au6_MOL_LR"]
    SIMULATION_PATHS = joinpath.(BASE_SRC_PATH, SIMULATIONS)

    R1 = 40.0u"kpc"
    R2 = 100.0u"kpc"
    R3 = 2.0u"kpc"

    lozano2024(SIMULATION_PATHS, LABELS, BASE_OUT_PATH, R1, R2, R3, LOGGING)

end
