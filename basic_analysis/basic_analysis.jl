cd(@__DIR__)

using Pkg

Pkg.activate("../")
Pkg.instantiate()

using CSV,
    CairoMakie,
    ColorSchemes,
    Colors,
    DataFrames,
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
    logging::Bool,
    r1::Unitful.Length,
    r2::Unitful.Length,
    r3::Unitful.Length;
    gas_evolution::Bool=false,
    label::String=basename(simulation_path),
)::Nothing

    ################################################################################################
    # Preparations and checks
    ################################################################################################

    # Create the output folders
    figures_path = mkpath(joinpath(base_out_path, "$(basename(simulation_path))_figures"))
    report_path  = mkpath(joinpath(base_out_path, "$(basename(simulation_path))_reports"))

    # Activate logging
    if logging
        log_file = open(joinpath(report_path, "logs.txt"), "w+")
        GalaxyInspector.setLogging!(logging; stream=log_file)
    end

    # Number of snapshots
    n_snapshots = GalaxyInspector.countSnapshot(simulation_path)

    # Number of snapshots as a string
    n_snaps_str = lpad(string(n_snapshots - 1), 3, "0")

    (
        !iszero(n_snapshots) ||
        throw(ArgumentError("basic_analysis: $(simulation_path) has no snapshots"))
    )

    # Check if the simulation is cosmological
    cosmological = GalaxyInspector.isSimCosmological(simulation_path)
    if cosmological
        filter_mode = :subhalo
    else
        filter_mode = :all_stellar
    end

    # Set default theme
    default_theme = merge(GalaxyInspector.DEFAULT_THEME, theme_latexfonts())

    # ################################################################################################
    # # Report files
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(log_file, "# Report files")
    #     println(log_file, "#"^100, "\n")
    # end

    # simulationReport([simulation_path]; output_path=report_path)

    # snapshotReport(
    #     [simulation_path],
    #     [n_snapshots];
    #     output_path=report_path,
    #     filter_mode,
    #     halo_idx=1,
    #     subhalo_rel_idx=1,
    # )

    # ################################################################################################
    # # Stellar surface density profile of the last snapshot, comparison with Agertz et al. (2021)
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(
    #         log_file,
    #         "# Stellar surface density profile of the last snapshot, comparison with Agertz et al. \
    #         (2021)",
    #     )
    #     println(log_file, "#"^100, "\n")
    # end

    # filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
    #     filter_mode,
    #     GalaxyInspector.plotParams(:stellar_mass).request,
    # )

    # grid = GalaxyInspector.CircularGrid(25.0u"kpc", 25)

    # y_label = GalaxyInspector.getLabel(L"\Sigma_\star", 0, u"Msun * kpc^-2")

    # plotSnapshot(
    #     [simulation_path],
    #     request,
    #     [lines!];
    #     output_path=figures_path,
    #     base_filename="stellar_mass_density",
    #     slice=n_snapshots,
    #     filter_function,
    #     da_functions=[GalaxyInspector.daProfile],
    #     da_args=[(:stellar_mass, grid)],
    #     da_kwargs=[(; flat=true, total=true, cumulative=false, density=true)],
    #     post_processing=GalaxyInspector.ppAgertz2021!,
    #     pp_kwargs=(;
    #         galaxies=[:all, "MW"],
    #         colors=[Makie.wong_colors()[4], Makie.wong_colors()[1]],
    #         linestyle=:dash,
    #         y_log=false,
    #     ),
    #     transform_box=true,
    #     translation,
    #     rotation,
    #     x_unit=u"kpc",
    #     xaxis_label="auto_label",
    #     y_unit=u"Msun * kpc^-2",
    #     yaxis_label=y_label,
    #     xaxis_var_name=L"r",
    #     yaxis_scale_func=log10,
    #     theme=Theme(
    #         palette=(linestyle=[:solid], color=[Makie.wong_colors()[2]]),
    #         Legend=(nbanks=1, valign=:top, padding=(0, 10, 0, 15)),
    #         Lines=(linewidth=3,),
    #         Scatter=(markersize=20,),
    #         Band=(alpha=0.7,),
    #     ),
    #     sim_labels=[label],
    # )

    # ################################################################################################
    # # SFR vs. physical time
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(log_file, "# SFR vs. physical time")
    #     println(log_file, "#"^100, "\n")
    # end

    # timeSeries(
    #     [simulation_path],
    #     :physical_time,
    #     :sfr;
    #     y_log=true,
    #     cumulative=false,
    #     fraction=false,
    #     slice=(:),
    #     output_path=figures_path,
    #     filter_mode,
    #     sim_labels=nothing,
    #     theme=Theme(palette=(color=[Makie.wong_colors()[2]],),),
    # )

    # ################################################################################################
    # # Stellar mass vs. physical time
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(log_file, "# Stellar mass vs. physical time")
    #     println(log_file, "#"^100, "\n")
    # end

    # timeSeries(
    #     [simulation_path],
    #     :physical_time,
    #     :stellar_mass;
    #     y_log=true,
    #     cumulative=false,
    #     fraction=false,
    #     slice=(:),
    #     output_path=figures_path,
    #     filter_mode,
    #     sim_labels=nothing,
    #     theme=Theme(palette=(color=[Makie.wong_colors()[2]],),),
    # )

    # ################################################################################################
    # # Circularity histogram of the last snapshot, radio separated
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(log_file, "# Circularity histogram of the last snapshot, radio separated")
    #     println(log_file, "#"^100, "\n")
    # end

    # plot_params = GalaxyInspector.plotParams(:stellar_circularity)
    # filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
    #     filter_mode,
    #     plot_params.request,
    # )
    # grid = GalaxyInspector.LinearGrid(-2.0, 2.0, 200)

    # da_ff = [
    #     dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r1), :zero),
    #     dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r2), :zero),
    #     dd -> GalaxyInspector.filterWithinSphere(dd, (r2, r1), :zero),
    # ]

    # r1_label = string(round(Int, ustrip(u"kpc", r1)))
    # r2_label = string(round(Int, ustrip(u"kpc", r2)))

    # plotSnapshot(
    #     [simulation_path, simulation_path, simulation_path],
    #     request,
    #     [lines!];
    #     output_path=figures_path,
    #     base_filename="circularity_histogram",
    #     slice=n_snapshots,
    #     filter_function,
    #     da_functions=[GalaxyInspector.daLineHistogram],
    #     da_args=[(:stellar_circularity, grid, :stellar)],
    #     da_kwargs=[
    #         (; filter_function=da_ff[1], norm=1),
    #         (; filter_function=da_ff[2], norm=1),
    #         (; filter_function=da_ff[3], norm=1),
    #     ],
    #     transform_box=true,
    #     translation,
    #     rotation,
    #     x_unit=plot_params.unit,
    #     x_exp_factor=plot_params.exp_factor,
    #     xaxis_label=plot_params.axis_label,
    #     xaxis_var_name=plot_params.var_name,
    #     yaxis_var_name=L"\mathrm{Normalized \,\, counts}",
    #     theme=Theme(
    #         size=(880, 880),
    #         palette=(
    #             color=[:gray65, :orangered2, :navy],
    #             linestyle=[:solid],
    #         ),
    #         Legend=(
    #             nbanks=1,
    #             halign=:left,
    #             valign=:top,
    #             padding=(15, 0, 0, 10),
    #             labelsize=25,
    #             rowgap=-4,
    #         ),
    #     ),
    #     sim_labels=[
    #         L"r \,\, \le \,\, %$(r1_label) \, \mathrm{kpc}",
    #         L"r \,\, \le \,\, %$(r2_label) \, \mathrm{kpc}",
    #         L"%$(r2_label) \, \mathrm{kpc} \,\, < \,\, r \,\, \le \,\, %$(r1_label) \, \mathrm{kpc}",
    #     ],
    # )

    # ################################################################################################
    # # Efficiency per free-fall time histograms of the last snapshot
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(log_file, "# Efficiency per free-fall time histograms of the last snapshot")
    #     println(log_file, "#"^100, "\n")
    # end

    # plot_params = GalaxyInspector.plotParams(:stellar_eff)

    # filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
    #     filter_mode,
    #     GalaxyInspector.mergeRequests(
    #         plot_params.request,
    #         GalaxyInspector.plotParams(:stellar_eff).request,
    #         GalaxyInspector.plotParams(:gas_eff).request,
    #         Dict(:stellar => ["GAGE"]),
    #     ),
    # )

    # grid = GalaxyInspector.LinearGrid(1.0e-4, 1.0, 100; log=true)

    # plotSnapshot(
    #     [simulation_path, simulation_path],
    #     request,
    #     [lines!];
    #     output_path=figures_path,
    #     base_filename="stellar_eff_histogram_all_stars",
    #     slice=n_snapshots,
    #     filter_function,
    #     da_functions=[GalaxyInspector.daLineHistogram, GalaxyInspector.daLineHistogram],
    #     da_args=[(:stellar_eff, grid, :stellar), (:gas_eff, grid, :gas)],
    #     da_kwargs=[
    #         (;
    #             filter_function=dd -> GalaxyInspector.intersectFilters(
    #                 GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r1), :zero),
    #             ),
    #             norm=0,
    #         )
    #     ],
    #     transform_box=true,
    #     translation,
    #     rotation,
    #     x_unit=plot_params.unit,
    #     x_exp_factor=plot_params.exp_factor,
    #     xaxis_label=plot_params.axis_label,
    #     xaxis_var_name=plot_params.var_name,
    #     yaxis_var_name=L"\mathrm{Normalized \,\, counts}",
    #     xaxis_scale_func=log10,
    #     theme=Theme(
    #         palette=(linestyle=[:solid], color=[Makie.wong_colors()[2], :black]),
    #         Legend=(nbanks=1, halign=:left, valign=:top, padding=(40, 0, 0, 0)),
    #     ),
    #     sim_labels=["All stars", "Gas"],
    # )

    # plotSnapshot(
    #     [simulation_path, simulation_path],
    #     request,
    #     [lines!];
    #     pf_kwargs=[(;)],
    #     output_path=figures_path,
    #     base_filename="stellar_eff_histogram_young_stars",
    #     slice=n_snapshots,
    #     filter_function,
    #     da_functions=[GalaxyInspector.daLineHistogram, GalaxyInspector.daLineHistogram],
    #     da_args=[(:stellar_eff, grid, :stellar), (:gas_eff, grid, :gas)],
    #     da_kwargs=[
    #         (;
    #             filter_function=dd -> GalaxyInspector.intersectFilters(
    #                 GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r1), :zero),
    #                 GalaxyInspector.filterOldStars(dd),
    #             ),
    #             norm=0,
    #         )
    #     ],
    #     transform_box=true,
    #     translation,
    #     rotation,
    #     x_unit=plot_params.unit,
    #     x_exp_factor=plot_params.exp_factor,
    #     xaxis_label=plot_params.axis_label,
    #     xaxis_var_name=plot_params.var_name,
    #     yaxis_var_name=L"\mathrm{Normalized \,\, counts}",
    #     xaxis_scale_func=log10,
    #     theme=Theme(
    #         palette=(linestyle=[:solid], color=[Makie.wong_colors()[2], :black]),
    #         Legend=(nbanks=1, halign=:left, valign=:top, padding=(40, 0, 0, 0)),
    #     ),
    #     sim_labels=["Young stars", "Gas"],
    # )

    ################################################################################################
    # Profiles for the last snapshot, comparison with Mollá et al. (2015)
    ################################################################################################

    for quantity in [
        # :stellar_area_density,
        :molecular_area_density,
        :sfr_area_density,
        :atomic_area_density,
        # :O_stellar_abundance,
        # :N_stellar_abundance,
        # :C_stellar_abundance,
    ]

        if logging
            println(log_file, "#"^100)
            println(
                log_file,
                "# $(quantity) profile for the last snapshot, comparison with Mollá et al. (2015)",
            )
            println(log_file, "#"^100, "\n")
        end

        compareMolla2015(
            [simulation_path],
            n_snapshots,
            quantity;
            output_path=joinpath(figures_path, "Molla2015"),
            filter_mode,
            sim_labels=[label],
            theme=Theme(
                size=(1500, 880),
                palette=(linestyle=[:solid],),
                Axis=(aspect=nothing,),
                Legend=(halign=:right, valign=:top),
            ),
        )

    end

    # ################################################################################################
    # # Resolved Kennicutt–Schmidt law of the last snapshot, scatter with square grid
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(
    #         log_file,
    #         "# Resolved Kennicutt–Schmidt law of the last snapshot, scatter with square grid",
    #     )
    #     println(log_file, "#"^100, "\n")
    # end

    # ###################
    # # Molecular KS law
    # ###################

    # if logging
    #     println(log_file, "#"^50)
    #     println(log_file, "# Molecular KS law")
    #     println(log_file, "#"^50, "\n")
    # end

    # kennicuttSchmidtLaw(
    #     [simulation_path],
    #     n_snapshots;
    #     quantity=:molecular_mass,
    #     reduce_grid=:square,
    #     grid_size=30.0u"kpc",
    #     bin_size=1.5u"kpc",
    #     gas_weights=nothing,
    #     post_processing=GalaxyInspector.ppSun2023!,
    #     pp_kwargs=(; color=Makie.wong_colors()[2]),
    #     fit=true,
    #     output_file=joinpath(figures_path, "_ks_law/sun2023_molecular.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 0),),
    #         Axis=(
    #             limits=(3.5, 8.5, -4.5, 0.5),
    #             xticks=4:1:8,
    #             yticks=-4:1:0,
    #         ),
    #     ),
    # )

    # kennicuttSchmidtLaw(
    #     [simulation_path],
    #     n_snapshots;
    #     quantity=:molecular_mass,
    #     reduce_grid=:square,
    #     grid_size=30.0u"kpc",
    #     bin_size=0.8u"kpc",
    #     gas_weights=nothing,
    #     post_processing=GalaxyInspector.ppLeroy2008!,
    #     pp_kwargs=(; color=Makie.wong_colors()[2]),
    #     fit=true,
    #     output_file=joinpath(figures_path, "_ks_law/leroy2008_molecular.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 0),),
    #         Axis=(
    #             limits=(3.5, 8.5, -4.5, 0.5),
    #             xticks=4:1:8,
    #             yticks=-4:1:0,
    #         ),
    #     ),
    # )

    # molecular_paths = [
    #     joinpath(figures_path, "_ks_law/leroy2008_molecular.png"),
    #     joinpath(figures_path, "_ks_law/sun2023_molecular.png"),
    # ]

    # GalaxyInspector.hvcatImages(
    #     1,
    #     molecular_paths;
    #     output_path=joinpath(figures_path, "_ks_law/molecular.png"),
    # )

    # ################
    # # Atomic KS law
    # ################

    # if logging
    #     println(log_file, "#"^50)
    #     println(log_file, "# Atomic KS law")
    #     println(log_file, "#"^50, "\n")
    # end

    # kennicuttSchmidtLaw(
    #     [simulation_path],
    #     n_snapshots;
    #     quantity=:atomic_mass,
    #     reduce_grid=:square,
    #     grid_size=30.0u"kpc",
    #     bin_size=0.6u"kpc",
    #     gas_weights=nothing,
    #     post_processing=GalaxyInspector.ppBigiel2010!,
    #     pp_kwargs=(; galaxy=:all, quantity=:atomic, color=Makie.wong_colors()[2]),
    #     fit=false,
    #     output_file=joinpath(figures_path, "_ks_law/bigiel2010_atomic.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 20),),
    #         Axis=(
    #             limits=(6.0, 8.0, -4.5, 0.5),
    #             xticks=6:0.5:8,
    #             yticks=-4:1:0,
    #         ),
    #     ),
    # )

    # kennicuttSchmidtLaw(
    #     [simulation_path],
    #     n_snapshots;
    #     quantity=:atomic_mass,
    #     reduce_grid=:square,
    #     grid_size=30.0u"kpc",
    #     bin_size=0.8u"kpc",
    #     gas_weights=nothing,
    #     post_processing=GalaxyInspector.ppLeroy2008!,
    #     pp_kwargs=(; quantity=:atomic, color=Makie.wong_colors()[2]),
    #     fit=false,
    #     output_file=joinpath(figures_path, "_ks_law/leroy2008_atomic.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 20),),
    #         Axis=(
    #             limits=(6.0, 8.0, -4.5, 0.5),
    #             xticks=6:0.5:8,
    #             yticks=-4:1:0,
    #         ),
    #     ),
    # )

    # atomic_paths = [
    #     joinpath(figures_path, "_ks_law/leroy2008_atomic.png"),
    #     joinpath(figures_path, "_ks_law/bigiel2010_atomic.png"),
    # ]

    # GalaxyInspector.hvcatImages(
    #     1,
    #     atomic_paths;
    #     output_path=joinpath(figures_path, "_ks_law/atomic.png"),
    # )

    # #################
    # # Neutral KS law
    # #################

    # if logging
    #     println(log_file, "#"^50)
    #     println(log_file, "# Neutral KS law")
    #     println(log_file, "#"^50, "\n")
    # end

    # kennicuttSchmidtLaw(
    #     [simulation_path],
    #     n_snapshots;
    #     quantity=:neutral_mass,
    #     reduce_grid=:square,
    #     grid_size=30.0u"kpc",
    #     bin_size=0.6u"kpc",
    #     gas_weights=nothing,
    #     post_processing=GalaxyInspector.ppBigiel2010!,
    #     pp_kwargs=(; galaxy=:all, quantity=:neutral, color=Makie.wong_colors()[2]),
    #     fit=false,
    #     output_file=joinpath(figures_path, "_ks_law/bigiel2010_neutral.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 20),),
    #         Axis=(
    #             limits=(6.5, 8.5, -4.5, 0.5),
    #             xticks=7:1:8,
    #             yticks=-4:1:0,
    #         ),
    #     ),
    # )

    # kennicuttSchmidtLaw(
    #     [simulation_path],
    #     n_snapshots;
    #     quantity=:neutral_mass,
    #     reduce_grid=:square,
    #     grid_size=30.0u"kpc",
    #     bin_size=0.8u"kpc",
    #     gas_weights=nothing,
    #     post_processing=GalaxyInspector.ppLeroy2008!,
    #     pp_kwargs=(; quantity=:neutral, color=Makie.wong_colors()[2]),
    #     fit=false,
    #     output_file=joinpath(figures_path, "_ks_law/leroy2008_neutral.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 20),),
    #         Axis=(
    #             limits=(6.5, 8.5, -4.5, 0.5),
    #             xticks=7:1:8,
    #             yticks=-4:1:0,
    #         ),
    #     ),
    # )

    # neutral_paths = [
    #     joinpath(figures_path, "_ks_law/leroy2008_neutral.png"),
    #     joinpath(figures_path, "_ks_law/bigiel2010_neutral.png"),
    # ]

    # GalaxyInspector.hvcatImages(
    #     1,
    #     neutral_paths;
    #     output_path=joinpath(figures_path, "_ks_law/neutral.png"),
    # )

    # ###################
    # # Total gas KS law
    # ###################

    # if logging
    #     println(log_file, "#"^50)
    #     println(log_file, "# Total gas KS law")
    #     println(log_file, "#"^50, "\n")
    # end

    # kennicuttSchmidtLaw(
    #     [simulation_path],
    #     n_snapshots;
    #     quantity=:gas_mass,
    #     reduce_grid=:square,
    #     grid_size=30.0u"kpc",
    #     bin_size=0.6u"kpc",
    #     gas_weights=nothing,
    #     post_processing=GalaxyInspector.ppBigiel2010!,
    #     pp_kwargs=(; galaxy=:all, quantity=:neutral, color=Makie.wong_colors()[2]),
    #     fit=false,
    #     output_file=joinpath(figures_path, "_ks_law/bigiel2010_total.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 20),),
    #         Axis=(
    #             limits=(6.8, 8.5, -4.5, 0.5),
    #             xticks=7:1:8,
    #             yticks=-4:1:0,
    #         ),
    #     ),
    # )

    # kennicuttSchmidtLaw(
    #     [simulation_path],
    #     n_snapshots;
    #     quantity=:gas_mass,
    #     reduce_grid=:square,
    #     grid_size=30.0u"kpc",
    #     bin_size=0.8u"kpc",
    #     gas_weights=nothing,
    #     post_processing=GalaxyInspector.ppLeroy2008!,
    #     pp_kwargs=(; quantity=:neutral, color=Makie.wong_colors()[2]),
    #     fit=false,
    #     output_file=joinpath(figures_path, "_ks_law/leroy2008_total.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 20),),
    #         Axis=(
    #             limits=(6.8, 8.5, -4.5, 0.5),
    #             xticks=7:1:8,
    #             yticks=-4:1:0,
    #         ),
    #     ),
    # )

    # total_paths = [
    #     joinpath(figures_path, "_ks_law/leroy2008_total.png"),
    #     joinpath(figures_path, "_ks_law/bigiel2010_total.png"),
    # ]

    # GalaxyInspector.hvcatImages(
    #     1,
    #     total_paths;
    #     output_path=joinpath(figures_path, "_ks_law/total.png"),
    # )

    # ###############################
    # # Final image with all KS laws
    # ###############################

    # gas_paths = joinpath.(
    #     figures_path,
    #     "_ks_law",
    #     ["molecular.png", "atomic.png", "neutral.png", "total.png"],
    # )

    # GalaxyInspector.hvcatImages(
    #     4,
    #     gas_paths;
    #     output_path=joinpath(figures_path, "ks_law.png"),
    # )

    # rm(joinpath(figures_path, "_ks_law"); recursive=true, force=true)




















    # ################################################################################################
    # # Stellar density maps of the last snapshot, face-on/edge-on projections
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(
    #         log_file,
    #         "# Stellar density maps of the last snapshot, face-on/edge-on projections",
    #     )
    #     println(log_file, "#"^100, "\n")
    # end

    # temp_folder = joinpath(figures_path, "_stellar_density_maps")

    # grid = GalaxyInspector.CubicGrid(r3, 400)

    # x_label = GalaxyInspector.getLabel("x", 0, u"kpc")
    # y_label = GalaxyInspector.getLabel("y", 0, u"kpc")
    # z_label = GalaxyInspector.getLabel("z", 0, u"kpc")

    # y_limits = [ceil(ustrip(u"kpc", r3) / 2.0), 12]
    # x_limits = ceil(ustrip(u"kpc", r3) / 2.0)

    # paths = joinpath.(temp_folder, ["stellar_mass_xy.jld2", "stellar_mass_xz.jld2"])

    # for projection_plane in [:xy, :xz]

    #     filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
    #         filter_mode,
    #         GalaxyInspector.plotParams(:stellar_mass).request,
    #     )

    #     GalaxyInspector.plotSnapshot(
    #         [simulation_path],
    #         request,
    #         [heatmap!];
    #         output_path=temp_folder,
    #         base_filename="stellar_mass_$(projection_plane)",
    #         slice=n_snapshots,
    #         filter_function,
    #         da_functions=[GalaxyInspector.daDensity2DProjection],
    #         da_args=[(grid, :stellar_mass, :particles)],
    #         da_kwargs=[(; projection_plane)],
    #         transform_box=true,
    #         translation,
    #         rotation,
    #         x_unit=u"kpc",
    #         y_unit=u"kpc",
    #         save_figures=false,
    #         backup_results=true,
    #     )

    # end

    # with_theme(default_theme) do

    #     f = Figure(size=(880, 1300), figure_padding=(5, 20, 0, 0))

    #     for (row, path) in pairs(paths)

    #         xaxis_v = row == 2

    #         ax = CairoMakie.Axis(
    #             f[row+1, 1];
    #             xlabel=x_label,
    #             ylabel=(row == 1 ? y_label : z_label),
    #             xminorticksvisible=xaxis_v,
    #             xticksvisible=xaxis_v,
    #             xlabelvisible=xaxis_v,
    #             xticklabelsvisible=xaxis_v,
    #             xticklabelsize=35,
    #             yticklabelsize=35,
    #             xticks=[-30, -20, -10, 0, 10, 20, 30],
    #             yticks=[-30, -20, -10, 0, 10, 20, 30],
    #             limits=(-x_limits, x_limits, -y_limits[row], y_limits[row]),
    #             aspect=DataAspect(),
    #         )

    #         jldopen(path, "r") do jld2_file

    #             x, y, z = jld2_file[first(keys(jld2_file))]["simulation_001"]

    #             pf = heatmap!(ax, x, y, z; colorrange=(6.5, 10.5))

    #             if row == 1

    #                 Colorbar(
    #                     f[row, 1],
    #                     pf,
    #                     label=L"\log_{10} \, \Sigma_* \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    #                     ticklabelsize=28,
    #                     ticks=7:0.5:10,
    #                     vertical=false,
    #                 )

    #             end

    #         end

    #     end

    #     rowsize!(f.layout, 3, Relative(0.3f0))

    #     Makie.save(joinpath(figures_path, "stellar_density_maps.png"), f)

    # end

    # rm(temp_folder; recursive=true)

    # ##########################################################################
    # # Gas components density map (last snapshot, face-on/edge-on projections)
    # ##########################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(
    #         log_file,
    #         "# Gas components density map (last snapshot, face-on/edge-on projections)",
    #     )
    #     println(log_file, "#"^100, "\n")
    # end

    # temp_folder = joinpath(figures_path, "_density_maps")

    # grid = GalaxyInspector.CubicGrid(r3, 400)
    # projections = [:xy, :xz]
    # quantities = [:gas_mass, :molecular_mass, :atomic_mass, :ionized_mass, :dust_mass]

    # colorbar_labels = [
    #     L"\log_{10} \, \Sigma_\mathrm{gas} \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    #     L"\log_{10} \, \Sigma_\mathrm{H_2} \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    #     L"\log_{10} \, \Sigma_\mathrm{HI} \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    #     L"\log_{10} \, \Sigma_\mathrm{HII} \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    #     L"\log_{10} \, \Sigma_\mathrm{dust} \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    # ]

    # x_label = GalaxyInspector.getLabel("x", 0, u"kpc")
    # y_label = GalaxyInspector.getLabel("y", 0, u"kpc")
    # z_label = GalaxyInspector.getLabel("z", 0, u"kpc")
    # n_rows = length(projections)
    # n_cols = length(quantities)
    # x_size = 2100
    # y_size = (x_size / n_cols) * n_rows + 180.0

    # paths = Vector{String}(undef, n_rows * n_cols)

    # for (j, quantity) in pairs(quantities)

    #     for (i, projection_plane) in pairs(projections)

    #         filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
    #             cosmological ? :all_subhalo : :all_stellar,
    #             GalaxyInspector.plotParams(quantity).request,
    #         )

    #         GalaxyInspector.plotSnapshot(
    #             [simulation_path],
    #             request,
    #             [heatmap!];
    #             output_path=temp_folder,
    #             base_filename="$(quantity)_$(projection_plane)",
    #             slice=n_snapshots,
    #             filter_function,
    #             da_functions=[GalaxyInspector.daDensity2DProjection],
    #             da_args=[(grid, quantity, :cells)],
    #             da_kwargs=[(; projection_plane)],
    #             transform_box=true,
    #             translation,
    #             rotation,
    #             x_unit=u"kpc",
    #             y_unit=u"kpc",
    #             save_figures=false,
    #             backup_results=true,
    #         )

    #         paths[j + n_cols * (i - 1)] = joinpath(
    #             temp_folder,
    #             "$(quantity)_$(projection_plane).jld2",
    #         )

    #     end

    # end

    # with_theme(default_theme) do

    #     f = Figure(size=(x_size, y_size), figure_padding=(0, 0, 0, 0))

    #     for (idx, path) in enumerate(paths)

    #         row = ceil(Int, idx / n_cols)
    #         col = mod1(idx, n_cols)

    #         xaxis_v = row == 2
    #         yaxis_v = col == 1

    #         ax = CairoMakie.Axis(
    #             f[row+1, col];
    #             xlabel=x_label,
    #             ylabel=(row == 1 ? y_label : z_label),
    #             xminorticksvisible=xaxis_v,
    #             xticksvisible=xaxis_v,
    #             xlabelvisible=xaxis_v,
    #             xticklabelsvisible=xaxis_v,
    #             yminorticksvisible=yaxis_v,
    #             yticksvisible=yaxis_v,
    #             ylabelvisible=yaxis_v,
    #             yticklabelsvisible=yaxis_v,
    #             xticklabelsize=28,
    #             yticklabelsize=28,
    #             xticks=[-30, -20, -10, 0, 10, 20, 30],
    #             yticks=[-30, -20, -10, 0, 10, 20, 30],
    #         )

    #         jldopen(path, "r") do jld2_file

    #             x, y, z = jld2_file[first(keys(jld2_file))]["simulation_001"]

    #             pf = heatmap!(ax, x, y, z)

    #             if row == 1

    #                 Colorbar(
    #                     f[row, col],
    #                     pf,
    #                     label=colorbar_labels[col],
    #                     ticklabelsize=23,
    #                     ticks=2:1:8,
    #                     vertical=false,
    #                 )

    #             end

    #         end

    #         colgap!(f.layout, 50)
    #         colsize!(f.layout, col, Makie.Fixed(352.0f0))

    #     end

    #     Makie.save(joinpath(figures_path, "gas_density_maps.png"), f)

    # end

    # rm(temp_folder; recursive=true)

    # ##############################################################
    # # Face-on density maps of different quantities (last snapshot)
    # ##############################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(log_file, "# Face-on density maps of different quantities (last snapshot)")
    #     println(log_file, "#"^100, "\n")
    # end

    # temp_folder = joinpath(figures_path, "_face_on_maps")

    # densityMap(
    #     [simulation_path],
    #     n_snapshots;
    #     quantities=[:gas_mass],
    #     types=[:cells],
    #     output_path=temp_folder,
    #     filter_mode,
    #     projection_planes=[:xy],
    #     box_size=r3,
    #     pixel_length=r3 / 400.0,
    #     theme=Theme(
    #         size=(910, 760),
    #         figure_padding=(20, 20, 30, 20),
    #         Colorbar=(
    #             label=L"\mathrm{log}_{10} \, \Sigma_\mathrm{gas} \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    #         ),
    #         Axis=(xticks=[-30, -20, -10, 0, 10, 20, 30], yticks=[-30, -20, -10, 0, 10, 20, 30]),
    #     ),
    #     title="Gas density",
    #     colorbar=true,
    #     colorrange=nothing,
    # )

    # densityMap(
    #     [simulation_path],
    #     n_snapshots;
    #     quantities=[:stellar_mass],
    #     types=[:particles],
    #     output_path=temp_folder,
    #     filter_mode,
    #     projection_planes=[:xy],
    #     box_size=r3,
    #     pixel_length=r3 / 400.0,
    #     theme=Theme(
    #         size=(910, 760),
    #         figure_padding=(20, 20, 30, 20),
    #         Colorbar=(
    #             label=L"\mathrm{log}_{10} \, \Sigma_\star \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    #         ),
    #         Axis=(xticks=[-30, -20, -10, 0, 10, 20, 30], yticks=[-30, -20, -10, 0, 10, 20, 30]),
    #     ),
    #     title="Stellar density",
    #     colorbar=true,
    #     colorrange=nothing,
    # )

    # temperatureMap(
    #     [simulation_path],
    #     n_snapshots;
    #     type=:cells,
    #     output_path=temp_folder,
    #     filter_mode,
    #     projection_planes=[:xy],
    #     box_size=r3,
    #     pixel_length=r3 / 400.0,
    #     theme=Theme(
    #         size=(910, 760),
    #         figure_padding=(20, 20, 30, 20),
    #         Colorbar=(label=L"\mathrm{log}_{10} \, T \,\, [\mathrm{K}]",),
    #         Axis=(xticks=[-30, -20, -10, 0, 10, 20, 30], yticks=[-30, -20, -10, 0, 10, 20, 30]),
    #     ),
    #     title="Gas temperature",
    #     colorbar=true,
    #     colorrange=nothing,
    # )

    # metallicityMap(
    #     [simulation_path],
    #     n_snapshots;
    #     components=[:gas],
    #     types=[:cells],
    #     output_path=temp_folder,
    #     filter_mode,
    #     projection_planes=[:xy],
    #     box_size=r3,
    #     pixel_length=r3 / 400.0,
    #     theme=Theme(
    #         size=(910, 760),
    #         figure_padding=(20, 20, 30, 20),
    #         Colorbar=(label=L"\mathrm{log}_{10} \, Z_\mathrm{gas} \, [\mathrm{Z_\odot}]",),
    #         Axis=(xticks=[-30, -20, -10, 0, 10, 20, 30], yticks=[-30, -20, -10, 0, 10, 20, 30]),
    #     ),
    #     title="Gas metallicity",
    #     colorbar=true,
    #     colorrange=nothing,
    # )

    # metallicityMap(
    #     [simulation_path],
    #     n_snapshots;
    #     components=[:stellar],
    #     types=[:particles],
    #     output_path=temp_folder,
    #     filter_mode,
    #     projection_planes=[:xy],
    #     box_size=r3,
    #     pixel_length=r3 / 400.0,
    #     theme=Theme(
    #         size=(910, 760),
    #         figure_padding=(20, 20, 30, 20),
    #         Colorbar=(label=L"\mathrm{log}_{10} \, Z_\star \, [\mathrm{Z_\odot}]",),
    #         Axis=(xticks=[-30, -20, -10, 0, 10, 20, 30], yticks=[-30, -20, -10, 0, 10, 20, 30]),
    #     ),
    #     title="Stellar metallicity",
    #     colorbar=true,
    #     colorrange=(-0.1, 0.1),
    # )

    # gasSFRMap(
    #     [simulation_path],
    #     n_snapshots;
    #     type=:cells,
    #     output_path=temp_folder,
    #     filter_mode,
    #     projection_planes=[:xy],
    #     box_size=r3,
    #     pixel_length=r3 / 400.0,
    #     theme=Theme(
    #         size=(910, 760),
    #         figure_padding=(20, 20, 30, 20),
    #         Colorbar=(
    #             label=L"\mathrm{log}_{10} \, SFR_\mathrm{gas} \, [\mathrm{M_\odot \, yr^{-1}}]",
    #         ),
    #         Axis=(xticks=[-30, -20, -10, 0, 10, 20, 30], yticks=[-30, -20, -10, 0, 10, 20, 30]),
    #     ),
    #     title="Gas SFR",
    #     colorbar=true,
    #     colorrange=nothing,
    # )

    # GalaxyInspector.hvcatImages(
    #     3,
    #     joinpath.(
    #         temp_folder,
    #         [
    #             "$(basename(simulation_path))_gas_mass_xy_density_map_snap_$(n_snaps_str).png",
    #             "$(basename(simulation_path))_stellar_mass_xy_density_map_snap_$(n_snaps_str).png",
    #             "$(basename(simulation_path))_xy_gas_sfr_map_snap_$(n_snaps_str).png",
    #             "$(basename(simulation_path))_xy_temperature_map_snap_$(n_snaps_str).png",
    #             "$(basename(simulation_path))_gas_xy_metallicity_map_snap_$(n_snaps_str).png",
    #             "$(basename(simulation_path))_stellar_xy_metallicity_map_snap_$(n_snaps_str).png",
    #         ]
    #     );
    #     output_path=joinpath(figures_path, "face_on_maps.png"),
    # )

    # rm(temp_folder, recursive=true, force=true)


    # if gas_evolution

    #     #################################################################################
    #     # Evolution of the masses and of the fractions, for the different gas components
    #     # (last snapshot, within a sphere of radius r1)
    #     #################################################################################

    #     if logging
    #         println(log_file, "#"^100)
    #         println(
    #             log_file,
    #             "# Evolution of the masses and of the fractions, for the different gas components",
    #             "\n# (last snapshot, within a sphere of radius r1)",
    #         )
    #         println(log_file, "#"^100, "\n")
    #     end

    #     r1_label = string(round(Int, ustrip(u"kpc", r1))) * "kpc"

    #     quantities = [:ionized_fraction, :atomic_fraction, :molecular_fraction]
    #     sim_labels = ["Ionized fraction", "Atomic fraction", "Molecular fraction"]

    #     x_plot_params = GalaxyInspector.plotParams(:physical_time)
    #     y_plot_params = GalaxyInspector.plotParams(:generic_fraction)

    #     temp_folder = joinpath(figures_path, "_gas_evolution")

    #     # Starts at ~200 Myr to ignore initial very low fractions
    #     initial_snap = GalaxyInspector.findClosestSnapshot(simulation_path, 0.2u"Gyr")

    #     plotTimeSeries(
    #         fill(simulation_path, length(quantities)),
    #         [lines!];
    #         output_path=temp_folder,
    #         filename="gas_fractions_evolution",
    #         slice=initial_snap:n_snapshots,
    #         da_functions=[GalaxyInspector.daEvolution],
    #         da_args=[(:physical_time, quantity) for quantity in quantities],
    #         da_kwargs=[
    #             (;
    #                 filter_mode,
    #                 extra_filter=dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r1), :zero),
    #                 scaling=identity,
    #             ),
    #         ],
    #         x_unit=x_plot_params.unit,
    #         y_unit=y_plot_params.unit,
    #         x_exp_factor=x_plot_params.exp_factor,
    #         y_exp_factor=y_plot_params.exp_factor,
    #         xaxis_label=x_plot_params.axis_label,
    #         yaxis_label=y_plot_params.axis_label,
    #         xaxis_var_name=x_plot_params.var_name,
    #         yaxis_var_name=y_plot_params.var_name,
    #         save_figures=false,
    #         backup_results=true,
    #     )

    #     quantities = [:stellar_mass, :gas_mass, :ionized_mass, :atomic_mass, :molecular_mass]
    #     sim_labels = [
    #         "Stellar mass",
    #         "Gas mass",
    #         "Ionized mass",
    #         "Atomic mass",
    #         "Molecular mass",
    #     ]

    #     y_plot_params = GalaxyInspector.plotParams(:generic_mass)

    #     plotTimeSeries(
    #         fill(simulation_path, length(quantities)),
    #         [lines!];
    #         output_path=temp_folder,
    #         filename="gas_masses_evolution",
    #         slice=initial_snap:n_snapshots,
    #         da_functions=[GalaxyInspector.daEvolution],
    #         da_args=[(:physical_time, quantity) for quantity in quantities],
    #         da_kwargs=[
    #             (;
    #                 filter_mode,
    #                 extra_filter=dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", r1), :zero),
    #                 scaling=identity,
    #             ),
    #         ],
    #         x_unit=x_plot_params.unit,
    #         y_unit=y_plot_params.unit,
    #         x_exp_factor=x_plot_params.exp_factor,
    #         y_exp_factor=y_plot_params.exp_factor,
    #         xaxis_label=x_plot_params.axis_label,
    #         yaxis_label=y_plot_params.axis_label,
    #         xaxis_var_name=x_plot_params.var_name,
    #         yaxis_var_name=y_plot_params.var_name,
    #         save_figures=false,
    #         backup_results=true,
    #     )

    #     paths = joinpath.(
    #         temp_folder,
    #         [
    #             "gas_masses_evolution.jld2",
    #             "gas_fractions_evolution.jld2",
    #         ],
    #     )

    #     current_theme = merge(Theme(palette=(linestyle=[:solid],),), default_theme)

    #     with_theme(current_theme) do

    #         f = Figure(size=(880, 1200),)

    #         ax_1 = CairoMakie.Axis(
    #             f[1, 1];
    #             xlabel=L"t \, [\mathrm{Gyr}]",
    #             ylabel=L"\log_{10} \, M \, [\mathrm{10^{10} \, M_\odot}]",
    #             xminorticksvisible=false,
    #             xticksvisible=false,
    #             xlabelvisible=false,
    #             xticklabelsvisible=false,
    #             limits=(nothing, nothing, -2.2, nothing),
    #         )

    #         jldopen(paths[1], "r") do jld2_file

    #             first_address = first(keys(jld2_file))

    #             x_s, y_s = jld2_file[first_address]["simulation_001"]
    #             x_h, y_h = jld2_file[first_address]["simulation_002"]
    #             x_i, y_i = jld2_file[first_address]["simulation_003"]
    #             x_a, y_a = jld2_file[first_address]["simulation_004"]
    #             x_m, y_m = jld2_file[first_address]["simulation_005"]

    #             lines!(ax_1, x_s, log10.(y_s); color=Makie.wong_colors()[2], label="Stellar mass")
    #             lines!(ax_1, x_h, log10.(y_h); color=:black, label="Gas mass")
    #             lines!(ax_1, x_i, log10.(y_i); color=Makie.wong_colors()[1], label="Ionized mass" )
    #             lines!(ax_1, x_a, log10.(y_a); color=Makie.wong_colors()[4], label="Atomic mass")
    #             lines!(ax_1, x_m, log10.(y_m); color=Makie.wong_colors()[3], label="Molecular mass")

    #         end

    #         axislegend(ax_1, position=:rb, framevisible=false, nbanks=1)

    #         ax_2 = CairoMakie.Axis(
    #             f[2, 1];
    #             xlabel=L"t \, [\mathrm{Gyr}]",
    #             ylabel=L"\log_{10} \, f",
    #             aspect=nothing,
    #             limits=(nothing, nothing, -3, nothing),
    #         )

    #         jldopen(paths[2], "r") do jld2_file

    #             first_address = first(keys(jld2_file))

    #             x_i, y_i = jld2_file[first_address]["simulation_001"]
    #             x_a, y_a = jld2_file[first_address]["simulation_002"]
    #             x_m, y_m = jld2_file[first_address]["simulation_003"]

    #             lines!(ax_2, x_i, log10.(y_i); color=Makie.wong_colors()[1])
    #             lines!(ax_2, x_a, log10.(y_a); color=Makie.wong_colors()[4])
    #             lines!(ax_2, x_m, log10.(y_m); color=Makie.wong_colors()[3])

    #         end

    #         linkxaxes!(ax_1, ax_2)
    #         rowsize!(f.layout, 1, Relative(2 / 3))
    #         colsize!(f.layout, 1, Makie.Fixed(pixelarea(ax_1.scene)[].widths[2]))

    #         Makie.save(joinpath(figures_path, "gas_evolution_inside_$(r1_label).png"), f)

    #     end

    #     rm(temp_folder, recursive=true, force=true)

    # end

    ##############
    # Close files
    ##############

    logging && close(log_file)

    return nothing

end

function comparison(
    simulation_paths::Vector{String},
    base_out_path::String,
    r1::Unitful.Length,
    logging::Bool;
    labels::Vector{String}=basename.(simulation_paths),
)::Nothing

    # Create the output folders
    report_path = mkpath(joinpath(base_out_path, "comparison"))
    figures_path = mkpath(joinpath(report_path, "figures"))

    # Activate logging
    if logging
        log_file = open(joinpath(report_path, "logs.txt"), "w+")
        GalaxyInspector.setLogging!(logging; stream=log_file)
    end

    temp_folder = joinpath(figures_path, "_mass_evolution")

    # Select the closest snapshot to readshift 0
    snaps = GalaxyInspector.findClosestSnapshot.(simulation_paths, 14.0u"Gyr")

    # Starts at ~200 Myr to ignore initial very low fractions
    initial_snaps = GalaxyInspector.findClosestSnapshot.(simulation_paths, 0.2u"Gyr")

    n_sims= length(simulation_paths)
    @assert n_sims == length(labels) "Number of simulations and labels must match."

    #########################################################################################
    # Compute the evolution of the total mass of the different gas components and of the SFR
    # (within a sphere of radius r1)
    #########################################################################################

    if logging
        println(log_file, "#"^100)
        println(
            log_file,
            "# Compute the evolution of the total mass of the different gas components and of the SFR",
            "\n(within a sphere of radius r1)",
        )
        println(log_file, "#"^100, "\n")
    end

    quantities = [:ionized_mass, :atomic_mass, :molecular_mass]
    n_panels   = length(quantities) + 1

    x_plot_params = GalaxyInspector.plotParams(:physical_time)

    for (simulation, z0_snap) in zip(simulation_paths, snaps)

        simulation_table = GalaxyInspector.makeSimulationTable(simulation)

        # Check if the simulation is cosmological
        cosmological = GalaxyInspector.isCosmological(
            first(skipmissing(simulation_table[!, :snapshot_paths]))
        )

        if cosmological
            filter_mode = :subhalo
        else
            filter_mode = :all_stellar
        end

        for quantity in quantities

            if logging
                println(log_file, "#"^100)
                println(
                    log_file,
                    "# Computing the evolution of $(quantity) for $(basename(simulation))",
                )
                println(log_file, "#"^100, "\n")
            end

            y_plot_params = GalaxyInspector.plotParams(quantity)

            plotTimeSeries(
                [simulation],
                [lines!];
                output_path=temp_folder,
                filename="$(quantity)_$(basename(simulation))",
                da_functions=[GalaxyInspector.daEvolution],
                da_args=[(:physical_time, quantity)],
                da_kwargs=[
                    (;
                        filter_mode,
                        extra_filter=dd -> GalaxyInspector.filterWithinSphere(
                            dd,
                            (0.0u"kpc", r1),
                            :zero,
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

        y_plot_params = GalaxyInspector.plotParams(:sfr)

        filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
            filter_mode,
            y_plot_params.request,
        )

        if logging
            println(log_file, "#"^100)
            println(
                log_file,
                "# Computing the evolution of the SFR for $(basename(simulation))",
            )
            println(log_file, "#"^100, "\n")
        end

        plotSnapshot(
            [simulation],
            request,
            [lines!];
            output_path=temp_folder,
            base_filename="sfr_$(basename(simulation))",
            slice=z0_snap,
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
            y_unit=y_plot_params.unit,
            x_exp_factor=x_plot_params.exp_factor,
            y_exp_factor=y_plot_params.exp_factor,
            save_figures=false,
            backup_results=true,
        )

    end

    ###################################################################################
    # Plot the evolution of the masses and of the SFR for the different gas components
    # (within a sphere of radius r1)
    ###################################################################################

    jld2_paths = joinpath.(
        temp_folder,
        vcat(
            ["sfr_$(label).jld2" for label in basename.(simulation_paths)],
            ["ionized_mass_$(label).jld2" for label in basename.(simulation_paths)],
            ["atomic_mass_$(label).jld2" for label in basename.(simulation_paths)],
            ["molecular_mass_$(label).jld2" for label in basename.(simulation_paths)],
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

    current_theme = merge(
        Theme(palette=(linestyle=[:solid],),),
        GalaxyInspector.DEFAULT_THEME,
        theme_latexfonts(),
    )

    colors = current_theme[:palette][:color][]
    y_lows = [-3.2, exp10(-2.5), exp10(-2.5), exp10(-10.5)]
    x_limits = (-0.1, 14.0)
    xticks = 0.0:2.0:14.0

    with_theme(current_theme) do

        f = Figure(size=(880, 1840), figure_padding=(1, 15, 5, 20))

        #########################
        # Plot the SFR evolution
        #########################

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
            ylabel=L"$\log_{10}$ %$(ylabel)",
            aspect=AxisAspect(1.7),
            xlabelvisible=false,
            xticklabelsvisible=false,
            xticks,
            yticks=-3:1:1,
            limits=(x_limits[1], x_limits[2], y_lows[1], nothing),
        )

        for (path, label, color) in zip(jld2_paths[1:n_sims], labels, colors)

            jldopen(path, "r") do jld2_file

                address = first(keys(jld2_file))

                x, y = jld2_file[address]["simulation_001"]

                lines!(ax, x, log10.(y); label, color)

            end

        end

        axislegend(ax, position=:rb, framevisible=false, nbanks=1)

        ##############################
        # Plot the gas mass evolution
        ##############################

        for (i, (quantity, y_low)) in enumerate(zip(quantities, y_lows[2:end]))

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
                xticks,
                limits=(x_limits[1], x_limits[2], y_low, nothing),
            )

            paths = jld2_paths[(n_sims*i+1):(n_sims*i+n_sims)]

            for (path, label, color, initial_snap) in zip(paths, labels, colors, initial_snaps)

                jldopen(path, "r") do jld2_file

                    address = first(keys(jld2_file))

                    x, y = jld2_file[address]["simulation_001"]

                    lines!(ax, x[initial_snap:end], y[initial_snap:end]; label, color)

                end

            end

        end

        Makie.save(joinpath(figures_path, "mass_evolution.png"), f)

    end

    rm(temp_folder; recursive=true)

    ##############
    # Close files
    ##############

    logging && close(log_file)

    return nothing

end

function (@main)(ARGS)

    # If logging into a file will be enable
    LOGGING = true

    # Output folder
    BASE_OUT_PATH = "./"

    # Characteristic radii
    R1 = 40.0u"kpc"
    R2 = 2.0u"kpc"
    R3 = 65u"kpc"

    # If the evolution of the masses and of the fractions will be done (slow to run)
    GAS_EVOLUTION = true

    SIMULATIONS = [
        "F:/simulations/isolated/SFM_064",
    ]

    LABELS = [
        "SFM_064",
    ]

    # SIMULATIONS = [
    #     "F:/simulations/current/test_dust_05",
    #     "F:/simulations/current/test_dust_08",
    #     "F:/simulations/current/test_dust_09",
    #     "F:/simulations/current/test_dust_10",
    # ]

    # LABELS = [
    #     "No LWB",
    #     "Incatasciato et al. 2023",
    #     "Ahn et al. 2009",
    #     "Visbal et al. 2014",
    # ]

    # if length(SIMULATIONS) > 1
    #     comparison(SIMULATIONS, BASE_OUT_PATH, R1, LOGGING; labels=LABELS)
    # end

    for (simulation, label) in zip(SIMULATIONS, LABELS)
        basic_analysis(
            simulation,
            BASE_OUT_PATH,
            LOGGING,
            R1,
            R2,
            R3;
            gas_evolution=GAS_EVOLUTION,
            label,
        )
    end

end
