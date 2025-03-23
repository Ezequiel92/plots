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

function proceedingAAA2024(
    simulation_paths::Vector{String},
    labels::Vector{String},
    base_out_path::String,
    r1::Unitful.Length,
    logging::Bool,
)::Nothing

    # Construct the necessary folders
    figures_path = mkpath(joinpath(base_out_path, "figures"))
    report_path = mkpath(joinpath(base_out_path, "reports"))

    # If requested, activate logging
    if logging
        log_file = open(joinpath(report_path, "logs.txt"), "w+")
        GalaxyInspector.setLogging!(logging; stream=log_file)
    end

    # Select the simulation paths
    Au6_MOL_path = simulation_paths[1]
    Au6_STD_path = simulation_paths[2]

    # Select the simulation labels
    Au6_MOL_label = labels[1]
    Au6_STD_label = labels[2]

    # Select the redshift 0 snapshot of each simulation
    Au6_MOL_z0_snap = GalaxyInspector.findClosestSnapshot(Au6_MOL_path, 14.0u"Gyr")
    Au6_STD_z0_snap = GalaxyInspector.findClosestSnapshot(Au6_STD_path, 14.0u"Gyr")

    z0_snaps = [Au6_MOL_z0_snap, Au6_STD_z0_snap]

    ######################################################################
    # SFR vs physical time (using an stellar age histogram at redshift 0)
    ######################################################################

    x_plot_params = GalaxyInspector.plotParams(:physical_time)
    y_plot_params = GalaxyInspector.plotParams(:sfr)

    filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
        :subhalo,
        y_plot_params.request,
    )

    plotSnapshot(
        [Au6_MOL_path, Au6_STD_path],
        request,
        [lines!];
        output_path=joinpath(figures_path),
        base_filename="sfr_vs_physical_time",
        output_format=".pdf",
        slice=min(z0_snaps...),
        filter_function,
        da_functions=[GalaxyInspector.daStellarHistory],
        da_kwargs=[
            (;
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
        y_unit=y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=y_plot_params.exp_factor,
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        yaxis_scale_func=log10,
        theme=Theme(
            fontsize=45,
            size=(880, 570),
            figure_padding=(1, 20, 5, 15),
            palette=(
                color=[Makie.wong_colors()[3], Makie.wong_colors()[4]],
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
                padding=(0, 25, 10, 0),
                labelsize=40,
                rowgap=-4,
            ),
        ),
        sim_labels=[Au6_MOL_label, Au6_STD_label],
    )

    #########################################################################################
    # Stellar density maps (face-on/edge-on projections with velocity fields, at redshift 0)
    #########################################################################################

    temp_folder = joinpath(figures_path, "_stellar_maps")

    box_size = 65.0u"kpc"
    limit = ustrip(u"kpc", box_size / 2.0)

    grid_hm = GalaxyInspector.CubicGrid(box_size, 350)
    grid_vf = GalaxyInspector.SquareGrid(box_size, 25)

    filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
        :subhalo,
        GalaxyInspector.plotParams(:stellar_mass).request,
    )

    for (simulation, slice) in zip(simulation_paths, z0_snaps)

        plotSnapshot(
            [simulation, simulation],
            request,
            [heatmap!, arrows!];
            output_path=joinpath(temp_folder, basename(simulation)),
            base_filename="$(basename(simulation))_xy",
            output_format=".pdf",
            slice,
            filter_function,
            da_functions=[GalaxyInspector.daDensity2DProjection, GalaxyInspector.daVelocityField],
            da_args=[(grid_hm, :stellar_mass, :particles), (grid_vf, :stars)],
            da_kwargs=[(; projection_plane=:xy),(; projection_plane=:xy),],
            transform_box=true,
            translation,
            rotation,
            x_unit=u"kpc",
            y_unit=u"kpc",
            save_figures=false,
            backup_results=true,
        )

        plotSnapshot(
            [simulation, simulation],
            request,
            [heatmap!, arrows!];
            output_path=joinpath(temp_folder, basename(simulation)),
            base_filename="$(basename(simulation))_xz",
            output_format=".pdf",
            slice,
            filter_function,
            da_functions=[GalaxyInspector.daDensity2DProjection, GalaxyInspector.daVelocityField],
            da_args=[(grid_hm, :stellar_mass, :particles), (grid_vf, :stars)],
            da_kwargs=[(; projection_plane=:xz),(; projection_plane=:xz),],
            transform_box=true,
            translation,
            rotation,
            x_unit=u"kpc",
            y_unit=u"kpc",
            save_figures=false,
            backup_results=true,
        )

    end

    paths = joinpath.(
        temp_folder,
        vcat(
            [
                [
                    "$(simulation)/$(simulation)_xy.jld2",
                    "$(simulation)/$(simulation)_xz.jld2",
                ] for simulation in basename.(simulation_paths)
            ]...,
        ),
    )

    current_theme = merge(
        Theme(
            fontsize=45,
            Heatmap=(
                colorrange=(5, 10),
                colormap=:nipy_spectral,
                nan_color=ColorSchemes.nipy_spectral[1],
            ),
        ),
        GalaxyInspector.DEFAULT_THEME,
        theme_latexfonts(),
    )

    with_theme(current_theme) do

        f = Figure(size=(1200, 740), figure_padding=(5, 5, 0, 15))

        ax_1 = CairoMakie.Axis(
            f[1, 1];
            ylabel=L"y \, [\mathrm{kpc}]",
            aspect=AxisAspect(1),
            limits=(-limit, limit, -limit, limit),
            xminorticksvisible=false,
            xticksvisible=false,
            xminorgridvisible=false,
            xlabelvisible=false,
            xticklabelsvisible=false,
            yticks=[-30, -20, -10, 0, 10, 20, 30],
        )

        ax_2 = CairoMakie.Axis(
            f[1, 2];
            ylabel=L"y \, [\mathrm{kpc}]",
            aspect=AxisAspect(1),
            limits=(-limit, limit, -limit, limit),
            xminorticksvisible=false,
            xticksvisible=false,
            xminorgridvisible=false,
            xlabelvisible=false,
            xticklabelsvisible=false,
            yminorticksvisible=false,
            yticksvisible=false,
            yminorgridvisible=false,
            ylabelvisible=false,
            yticklabelsvisible=false,
            yticks=[-30, -20, -10, 0, 10, 20, 30],
        )

        jldopen(paths[1], "r") do jld2_file

            first_address = first(keys(jld2_file))

            x_hm, y_hm, z_hm       = jld2_file[first_address]["simulation_001"]
            x_ar, y_ar, u_ar, v_ar = jld2_file[first_address]["simulation_002"]

            heatmap!(ax_1, x_hm, y_hm, z_hm)
            arrows!(ax_1, x_ar, y_ar, u_ar, v_ar)

            text!(
                ax_1,
                0.04,
                0.98;
                text=Au6_MOL_label,
                align=(:left, :top),
                color=:white,
                space=:relative,
                fontsize=50,
            )

        end

        jldopen(paths[3], "r") do jld2_file

            first_address = first(keys(jld2_file))

            x_hm, y_hm, z_hm       = jld2_file[first_address]["simulation_001"]
            x_ar, y_ar, u_ar, v_ar = jld2_file[first_address]["simulation_002"]

            pf = heatmap!(ax_2, x_hm, y_hm, z_hm)
            arrows!(ax_2, x_ar, y_ar, u_ar, v_ar)

            text!(
                ax_2,
                0.04,
                0.98;
                text=Au6_STD_label,
                align=(:left, :top),
                color=:white,
                space=:relative,
                fontsize=50,
            )

            Colorbar(
                f[1, 3];
                ticks=5:1:10,
                colormap=:nipy_spectral,
                colorrange=(5, 10),
                label=L"\mathrm{log}_{10} \, \Sigma_\star \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
            )

        end

        ax_3 = CairoMakie.Axis(
            f[2, 1];
            xlabel=L"x \, [\mathrm{kpc}]",
            ylabel=L"z \, [\mathrm{kpc}]",
            aspect=DataAspect(),
            limits=(-ustrip(u"kpc", box_size) / 2, ustrip(u"kpc", box_size) / 2, -10.0, 10.0,),
            xticks=[-30, -15, 0, 15, 30],
            yticks=[-10, 0, 10,],
        )

        ax_4 = CairoMakie.Axis(
            f[2, 2];
            xlabel=L"x \, [\mathrm{kpc}]",
            ylabel=L"z \, [\mathrm{kpc}]",
            aspect=DataAspect(),
            yminorticksvisible=false,
            yticksvisible=false,
            yminorgridvisible=false,
            ylabelvisible=false,
            yticklabelsvisible=false,
            limits=(-ustrip(u"kpc", box_size) / 2, ustrip(u"kpc", box_size) / 2, -10.0, 10.0,),
            xticks=[-30, -15, 0, 15, 30],
        )

        jldopen(paths[2], "r") do jld2_file

            first_address = first(keys(jld2_file))

            x_hm, y_hm, z_hm       = jld2_file[first_address]["simulation_001"]
            x_ar, y_ar, u_ar, v_ar = jld2_file[first_address]["simulation_002"]

            heatmap!(ax_3, x_hm, y_hm, z_hm)
            arrows!(ax_3, x_ar, y_ar, u_ar, v_ar)

        end

        jldopen(paths[4], "r") do jld2_file

            first_address = first(keys(jld2_file))

            x_hm, y_hm, z_hm       = jld2_file[first_address]["simulation_001"]
            x_ar, y_ar, u_ar, v_ar = jld2_file[first_address]["simulation_002"]

            heatmap!(ax_4, x_hm, y_hm, z_hm)
            arrows!(ax_4, x_ar, y_ar, u_ar, v_ar)

        end

        rowsize!(f.layout, 2, Relative(1/3.8))
        rowsize!(f.layout, 1, Makie.Fixed(pixelarea(ax_1.scene)[].widths[2]))

        colgap!(f.layout, 30)

        Makie.save(joinpath(figures_path, "stellar_maps.png"), f)

    end

    rm(temp_folder; force=true, recursive=true)

    ############################################################
    # Profiles of the gas surface density and of the fractions,
    # for the different gas components (at redshift 0)
    ############################################################

    temp_folder = joinpath(figures_path, "_gas_profiles")

    grid = GalaxyInspector.CircularGrid(r1, 15)

    quantities = [:gas_mass, :stellar_mass]
    plot_params = GalaxyInspector.plotParams(:generic_area_density)

    filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
        :subhalo,
        plot_params.request,
    )

    for (path, slice) in zip(simulation_paths, z0_snaps)

        plotSnapshot(
            fill(path, length(quantities)),
            request,
            [lines!];
            output_path=joinpath(temp_folder, basename(path)),
            base_filename="density_profiles",
            slice,
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

    end

    quantities = [:ionized_mass, :neutral_mass]
    plot_params = GalaxyInspector.plotParams(:generic_fraction)

    filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
        :subhalo,
        plot_params.request,
    )

    plotSnapshot(
        fill(Au6_MOL_path, length(quantities)),
        request,
        [lines!];
        output_path=joinpath(temp_folder, basename(Au6_MOL_path)),
        base_filename="fractions_profiles",
        slice=Au6_MOL_z0_snap,
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

    quantities = [:ionized_mass, :neutral_mass]

    plotSnapshot(
        fill(Au6_STD_path, length(quantities)),
        request,
        [lines!];
        output_path=joinpath(temp_folder, basename(Au6_STD_path)),
        base_filename="fractions_profiles",
        slice=Au6_STD_z0_snap,
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

    density_paths = joinpath.(temp_folder, basename.(simulation_paths), "density_profiles.jld2")
    fraction_paths = joinpath.(temp_folder, basename.(simulation_paths), "fractions_profiles.jld2")

    current_theme = merge(
        Theme(fontsize=45,),
        GalaxyInspector.DEFAULT_THEME,
        theme_latexfonts(),
    )

    with_theme(current_theme) do

        f = Figure(size=(880, 1600),)

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

        jldopen(density_paths[1], "r") do jld2_file

            first_address = first(keys(jld2_file))

            x_h, y_h = jld2_file[first_address]["simulation_001"]
            x_s, y_s = jld2_file[first_address]["simulation_002"]

            lines!(ax_1, x_h, y_h; linestyle=:solid, color=:black, label="Au6_MOL - Gas density")
            lines!(
                ax_1,
                x_s,
                y_s;
                linestyle=:solid,
                color=Makie.wong_colors()[2],
                label="Au6_MOL - Stellar density",
            )

        end

        jldopen(density_paths[2], "r") do jld2_file

            first_address = first(keys(jld2_file))

            x_h, y_h = jld2_file[first_address]["simulation_001"]
            x_s, y_s = jld2_file[first_address]["simulation_002"]

            lines!(ax_1, x_h, y_h; linestyle=:dash, color=:black, label="Au6_STD - Gas density")
            lines!(
                ax_1,
                x_s,
                y_s;
                linestyle=:dash,
                color=Makie.wong_colors()[2],
                label="Au6_STD - Stellar density",
            )

        end

        axislegend(ax_1, position=:rt, framevisible=false, fontsize=30, nbanks=1, rowgap=-5)

        ax_2 = CairoMakie.Axis(
            f[2, 1];
            xlabel=L"r \, [\mathrm{kpc}]",
            ylabel=L"f",
            aspect=nothing,
            xticks=0:5:40,
            limits=(nothing, nothing, -0.1, 1.3),
        )

        jldopen(fraction_paths[1], "r") do jld2_file

            first_address = first(keys(jld2_file))

            x_i, y_i = jld2_file[first_address]["simulation_001"]
            x_n, y_n = jld2_file[first_address]["simulation_002"]

            lines!(
                ax_2,
                x_i,
                y_i;
                linestyle=:solid,
                color=Makie.wong_colors()[1],
                label="Au6_MOL - Ionized fraction",
            )
            lines!(
                ax_2,
                x_n,
                y_n;
                linestyle=:solid,
                color=Makie.wong_colors()[4],
                label="Au6_MOL - Neutral fraction",
            )

        end

        jldopen(fraction_paths[2], "r") do jld2_file

            first_address = first(keys(jld2_file))

            x_i, y_i = jld2_file[first_address]["simulation_001"]
            x_n, y_n = jld2_file[first_address]["simulation_002"]

            lines!(
                ax_2,
                x_i,
                y_i;
                linestyle=:dash,
                color=Makie.wong_colors()[1],
                label="Au6_STD - Ionized fraction",
            )
            lines!(
                ax_2,
                x_n,
                y_n;
                linestyle=:dash,
                color=Makie.wong_colors()[4],
                label="Au6_STD - Neutral fraction",
            )

        end

        axislegend(ax_2, position=:rt, framevisible=false, fontsize=30, nbanks=1, rowgap=-5)

        linkxaxes!(ax_1, ax_2)
        rowsize!(f.layout, 1, Relative(1/2))
        colsize!(f.layout, 1, Makie.Fixed(pixelarea(ax_1.scene)[].widths[2]))

        Makie.save(joinpath(figures_path, "gas_density_and_fractions_profiles.pdf"), f)

    end

    rm(temp_folder; recursive=true)

    ################################################################################################
    # Close files
    ################################################################################################

    if logging
        close(log_file)
    end

    return nothing

end

function (@main)(ARGS)

    # If logging into a file will be enable
    LOGGING = true

    # Output folder
    BASE_OUT_PATH = "./"

    # Simulation folders
    BASE_SRC_PATH = "F:/simulations/lozano_2025/"
    SIMULATIONS = ["Au6_MOL_test18", "test_cosmo_volker_06"]
    SIMULATION_PATHS = joinpath.(BASE_SRC_PATH, SIMULATIONS)

    # Simulation labels
    LABELS = ["Au6_MOL", "Au6_STD"]

    # Characteristic radii
    R1 = 40.0u"kpc"

    proceedingAAA2024(SIMULATION_PATHS, LABELS, BASE_OUT_PATH, R1, LOGGING)

end
