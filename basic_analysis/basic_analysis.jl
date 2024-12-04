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
    ProgressMeter,
    QuadGK,
    Rotations,
    Statistics,
    StatsBase,
    Unitful,
    UnitfulAstro

push!(LOAD_PATH, "../../codes/GalaxyInspector/src/")
using GalaxyInspector

function basic_analysis(
    simulation_path::String,
    base_out_path::String,
    r1::Unitful.Length,
    r2::Unitful.Length,
    logging::Bool,
)::Nothing

    figures_path = mkpath(joinpath(base_out_path, "figures"))
    report_path = mkpath(joinpath(base_out_path, "reports"))
    log_file = open(joinpath(report_path, "logs.txt"), "w+")

    GalaxyInspector.setLogging!(logging; stream=log_file)

    n_snapshots = GalaxyInspector.countSnapshot(simulation_path)

    ###############
    # Report files
    ###############

    simulationReport([simulation_path]; output_path=report_path)

    snapshotReport(
        [simulation_path],
        n_snapshots;
        output_path=report_path,
        filter_mode=:subhalo,
        halo_idx=1,
        subhalo_rel_idx=1,
    )

    ##################################
    # Stellar surface density profile
    ##################################

    densityProfile(
        [simulation_path],
        n_snapshots,
        :stellar_mass;
        cumulative=false,
        yscale=log10,
        radius=r1,
        n_bins=40,
        output_path=figures_path,
        filter_mode=:subhalo,
        sim_labels=nothing,
        theme=Theme(palette=(color=[Makie.wong_colors()[2]],),),
    )

    #######################
    # SFR vs physical time
    #######################

    stellarHistory(
        [simulation_path],
        n_snapshots,
        :sfr;
        y_log=true,
        n_bins=70,
        output_path=figures_path,
        filter_mode=:subhalo,
        sim_labels=nothing,
        theme=Theme(palette=(color=[Makie.wong_colors()[2]],),),
    )

    ################################
    # Stellar mass vs physical time
    ################################

    stellarHistory(
        [simulation_path],
        n_snapshots,
        :stellar_mass;
        y_log=true,
        n_bins=70,
        output_path=figures_path,
        filter_mode=:subhalo,
        sim_labels=nothing,
        theme=Theme(palette=(color=[Makie.wong_colors()[2]],),),
    )

    ##########################################
    # Circularity histogram (radio separated)
    ##########################################

    plot_params = GalaxyInspector.plotParams(:stellar_circularity)
    filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
        :subhalo,
        plot_params.request,
    )
    grid = GalaxyInspector.LinearGrid(-2.0, 2.0, 200)

    da_ff = [
        dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r1), :zero),
        dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r2), :zero),
        dd -> GalaxyInspector.filterWithinSphere(dd, (r2, r1), :zero),
    ]

    r1_label = string(round(Int, ustrip(u"kpc", r1)))
    r2_label = string(round(Int, ustrip(u"kpc", r2)))

    plotSnapshot(
        [simulation_path, simulation_path, simulation_path],
        request,
        [lines!];
        output_path=figures_path,
        base_filename="circularity_histogram_radio_separated",
        output_format=".png",
        slice=n_snapshots,
        filter_function,
        da_functions=[GalaxyInspector.daLineHistogram],
        da_args=[(:stellar_circularity, grid, :stars)],
        da_kwargs=[
            (; filter_function=da_ff[1], norm=14232),
            (; filter_function=da_ff[2], norm=14232),
            (; filter_function=da_ff[3], norm=14232),
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
            size=(880, 880),
            palette=(
                color=[:gray65, :orangered2, :navy],
                linestyle=[:solid],
            ),
            Legend=(
                nbanks=1,
                halign=:left,
                valign=:top,
                padding=(15, 0, 0, 10),
                labelsize=25,
                rowgap=-4,
            ),
        ),
        sim_labels=[
            L"r \,\, \le \,\, %$(r1_label) \, \mathrm{kpc}",
            L"r \,\, \le \,\, %$(r2_label) \, \mathrm{kpc}",
            L"%$(r2_label) \, \mathrm{kpc} \,\, < \,\, r \,\, \le \,\, %$(r1_label) \, \mathrm{kpc}",
        ],
    )

    ##########################################
    # Efficiency per free-fall time histogram
    ##########################################

    lineHistogram(
        [simulation_path],
        n_snapshots,
        :stellar_eff,
        :stars,
        (1.0e-4, 1.0);
        n_bins=100,
        log=true,
        output_path=figures_path,
        filter_mode=:subhalo,
        extra_filter=dd -> GalaxyInspector.intersectFilters(
            GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r1), :zero),
            GalaxyInspector.filterOldStars(dd),
        ),
        ff_request=Dict(:gas => ["FRAC"]),
        sim_labels=nothing,
        theme=Theme(palette=(color=[Makie.wong_colors()[2]],),),
    )

    lineHistogram(
        [simulation_path],
        n_snapshots,
        :gas_eff,
        :gas,
        (1.0e-4, 1.0);
        n_bins=100,
        log=true,
        output_path=figures_path,
        filter_mode=:subhalo,
        extra_filter=dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r1), :zero),
        ff_request=Dict(:gas => ["FRAC"]),
        sim_labels=nothing,
        theme=Theme(palette=(color=[:black],),),
    )

    ################################################################################
    # Gas components density map (face-on/edge-on projections at the last snapshot)
    ################################################################################

    temp_folder = joinpath(figures_path, "density_maps")

    grid = GalaxyInspector.CubicGrid(65.0u"kpc", 400)
    projections = [:xy, :xz]
    quantities = [:stellar_mass, :gas_mass, :molecular_mass, :atomic_mass, :ionized_mass]
    types = [:particles, :cells, :cells, :cells, :cells]

    colorbar_labels = [
        L"\log_{10} \, \Sigma_\star \,\, [\mathrm{M_\odot \, kpc^{-2}}]"
        L"\log_{10} \, \Sigma_\mathrm{gas} \,\, [\mathrm{M_\odot \, kpc^{-2}}]"
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
    y_size = (x_size / n_cols) * n_rows + 180.0

    paths = joinpath.(
        temp_folder,
        [
            "stellar_mass_xy.jld2",
            "gas_mass_xy.jld2",
            "molecular_mass_xy.jld2",
            "atomic_mass_xy.jld2",
            "ionized_mass_xy.jld2",
            "stellar_mass_xz.jld2",
            "gas_mass_xz.jld2",
            "molecular_mass_xz.jld2",
            "atomic_mass_xz.jld2",
            "ionized_mass_xz.jld2",
        ]
    )

    for projection_plane in projections

        for (quantity, type) in zip(quantities, types)

            filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
                :all_subhalo,
                GalaxyInspector.plotParams(quantity).request,
            )

            GalaxyInspector.plotSnapshot(
                [simulation_path],
                request,
                [heatmap!];
                output_path=temp_folder,
                base_filename="$(quantity)_$(projection_plane)",
                slice=n_snapshots,
                filter_function,
                da_functions=[GalaxyInspector.daDensity2DProjection],
                da_args=[(grid, quantity, type)],
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

        f = Figure(size=(x_size, y_size), figure_padding=(5, 10, 0, 0))

        for (idx, path) in enumerate(paths)

            row = ceil(Int, idx / n_cols)
            col = mod1(idx, n_cols)

            xaxis_v = row == 2
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
                xticklabelsize=23,
                yticklabelsize=23,
                xticks=[-30, -20, -10, 0, 10, 20, 30],
                yticks=[-30, -20, -10, 0, 10, 20, 30],
            )

            jldopen(path, "r") do jld2_file

                x, y, z = jld2_file[first(keys(jld2_file))]["simulation_001"]

                pf = heatmap!(ax, x, y, z)

                if row == 1

                    Colorbar(
                        f[row, col],
                        pf,
                        label=colorbar_labels[col],
                        labelsize=30,
                        ticklabelsize=20,
                        vertical=false,
                    )

                end

            end

            colgap!(f.layout, 20)

        end

        Makie.save(joinpath(figures_path, "density_maps_grid.png"), f)

    end

    rm(temp_folder; recursive=true)

    close(log_file)

    return nothing

end

function (@main)(ARGS)

    LOGGING = true
    BASE_OUT_PATH = "./"
    SIMULATION_PATH = "F:/simulations/Au6_MOL_test/Au6_MOL_test18"
    R1 = 40.0u"kpc"
    R2 = 2.0u"kpc"

    basic_analysis(SIMULATION_PATH, BASE_OUT_PATH, R1, R2, LOGGING)

end
