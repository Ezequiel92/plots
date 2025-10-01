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
    logging::Bool;
    label::AbstractString=basename(simulation_path),
)::Nothing

    ################################################################################################
    # Preparations and checks
    ################################################################################################

    # Create the output folder
    output_path = mkpath(joinpath(base_out_path, "$(basename(simulation_path))"))

    # Activate logging
    if logging
        log_file = open(joinpath(output_path, "logs.txt"), "w+")
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
        filter_mode  = :subhalo
        extra_filter = dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", R1), :zero)
    else
        filter_mode  = :all
        extra_filter = GalaxyInspector.filterNothing
    end

    # Set default theme
    default_theme = merge(GalaxyInspector.DEFAULT_THEME, theme_latexfonts())

    # Characteristic radii
    R1 = 40.0u"kpc"
    R2 = 2.0u"kpc"
    R3 = 65.0u"kpc"

    # ################################################################################################
    # # Report files
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(log_file, "# Report files")
    #     println(log_file, "#"^100, "\n")
    # end

    # simulationReport([simulation_path]; output_path)

    # snapshotReport(
    #     [simulation_path],
    #     [n_snapshots];
    #     output_path,
    #     filter_mode,
    #     halo_idx=1,
    #     subhalo_rel_idx=1,
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
    #     y_log=false,
    #     cumulative=false,
    #     fraction=false,
    #     output_path,
    #     filter_mode,
    #     smooth=n_snapshots ÷ 5,
    #     sim_labels=nothing,
    #     theme=Theme(
    #         palette=(color=[Makie.wong_colors()[2]],),
    #         size=(1320, 880),
    #         Axis=(aspect=nothing,),
    #     ),
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
    #     y_log=false,
    #     cumulative=false,
    #     fraction=false,
    #     output_path,
    #     filter_mode,
    #     sim_labels=nothing,
    #     theme=Theme(
    #         palette=(color=[Makie.wong_colors()[2]],),
    #         size=(1320, 880),
    #         Axis=(aspect=nothing,),
    #     ),
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
    #     cosmological ? :subhalo : :all_stellar,
    #     GalaxyInspector.plotParams(:stellar_mass).request,
    # )

    # grid = GalaxyInspector.CircularGrid(25.0u"kpc", 25)

    # y_label = GalaxyInspector.getLabel(L"\Sigma_\star", 0, u"Msun * kpc^-2")

    # plotSnapshot(
    #     [simulation_path],
    #     request,
    #     [lines!];
    #     output_path,
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
    #     dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", R1), :zero),
    #     dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", R2), :zero),
    #     dd -> GalaxyInspector.filterWithinSphere(dd, (R2, R1), :zero),
    # ]

    # R1_label = string(round(Int, ustrip(u"kpc", R1)))
    # R2_label = string(round(Int, ustrip(u"kpc", R2)))

    # plotSnapshot(
    #     [simulation_path, simulation_path, simulation_path],
    #     request,
    #     [lines!];
    #     output_path,
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
    #         L"r \,\, \le \,\, %$(R1_label) \, \mathrm{kpc}",
    #         L"r \,\, \le \,\, %$(R2_label) \, \mathrm{kpc}",
    #         L"%$(R2_label) \, \mathrm{kpc} \,\, < \,\, r \,\, \le \,\, %$(R1_label) \, \mathrm{kpc}",
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
    #     output_path,
    #     base_filename="stellar_eff_histogram_all_stars",
    #     slice=n_snapshots,
    #     filter_function,
    #     da_functions=[GalaxyInspector.daLineHistogram, GalaxyInspector.daLineHistogram],
    #     da_args=[(:stellar_eff, grid, :stellar), (:gas_eff, grid, :gas)],
    #     da_kwargs=[(; norm=0)],
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
    #     output_path,
    #     base_filename="stellar_eff_histogram_young_stars",
    #     slice=n_snapshots,
    #     filter_function,
    #     da_functions=[GalaxyInspector.daLineHistogram, GalaxyInspector.daLineHistogram],
    #     da_args=[(:stellar_eff, grid, :stellar), (:gas_eff, grid, :gas)],
    #     da_kwargs=[(; filter_function=GalaxyInspector.filterByStellarAge, norm=0)],
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

    # ################################################################################################
    # # Profiles for the last snapshot, comparison with Mollá et al. (2015)
    # ################################################################################################

    # for quantity in [
    #     :stellar_area_density,
    #     :molecular_area_density,
    #     :sfr_area_density,
    #     :atomic_area_density,
    #     :O_stellar_abundance,
    #     :N_stellar_abundance,
    #     :C_stellar_abundance,
    # ]

    #     if logging
    #         println(log_file, "#"^100)
    #         println(
    #             log_file,
    #             "# $(quantity) profile for the last snapshot, comparison with Mollá et al. (2015)",
    #         )
    #         println(log_file, "#"^100, "\n")
    #     end

    #     compareMolla2015(
    #         [simulation_path],
    #         n_snapshots,
    #         quantity;
    #         output_path=joinpath(output_path, "Molla2015"),
    #         filter_mode=cosmological ? :subhalo : :all_stellar,
    #         sim_labels=[label],
    #         theme=Theme(
    #             size=(1500, 880),
    #             figure_padding=(10, 15, 5, 15),
    #             palette=(linestyle=[:solid],),
    #             Axis=(aspect=nothing,),
    #             Legend=(halign=:right, valign=:top, nbanks=1),
    #         ),
    #     )

    # end

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
    #     output_file=joinpath(output_path, "_ks_law/sun2023_molecular.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 0),),
    #         Axis=(
    #             limits=(4.5, 9.5, -4.5, 0.5),
    #             xticks=5:1:9,
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
    #     output_file=joinpath(output_path, "_ks_law/leroy2008_molecular.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 0),),
    #         Axis=(
    #             limits=(4.5, 9.5, -4.5, 0.5),
    #             xticks=5:1:9,
    #             yticks=-4:1:0,
    #         ),
    #     ),
    # )

    # molecular_paths = [
    #     joinpath(output_path, "_ks_law/leroy2008_molecular.png"),
    #     joinpath(output_path, "_ks_law/sun2023_molecular.png"),
    # ]

    # GalaxyInspector.hvcatImages(
    #     1,
    #     molecular_paths;
    #     output_path=joinpath(output_path, "_ks_law/molecular.png"),
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
    #     output_file=joinpath(output_path, "_ks_law/bigiel2010_atomic.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 20),),
    #         Axis=(
    #             limits=(5.8, 8.2, -4.5, 0.5),
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
    #     output_file=joinpath(output_path, "_ks_law/leroy2008_atomic.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 20),),
    #         Axis=(
    #             limits=(5.8, 8.2, -4.5, 0.5),
    #             xticks=6:0.5:8,
    #             yticks=-4:1:0,
    #         ),
    #     ),
    # )

    # atomic_paths = [
    #     joinpath(output_path, "_ks_law/leroy2008_atomic.png"),
    #     joinpath(output_path, "_ks_law/bigiel2010_atomic.png"),
    # ]

    # GalaxyInspector.hvcatImages(
    #     1,
    #     atomic_paths;
    #     output_path=joinpath(output_path, "_ks_law/atomic.png"),
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
    #     output_file=joinpath(output_path, "_ks_law/bigiel2010_neutral.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 20),),
    #         Axis=(
    #             limits=(6.4, 8.6, -4.5, 0.5),
    #             xticks=6.5:0.5:8.5,
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
    #     output_file=joinpath(output_path, "_ks_law/leroy2008_neutral.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 20),),
    #         Axis=(
    #             limits=(6.4, 8.6, -4.5, 0.5),
    #             xticks=6.5:0.5:8.5,
    #             yticks=-4:1:0,
    #         ),
    #     ),
    # )

    # neutral_paths = [
    #     joinpath(output_path, "_ks_law/leroy2008_neutral.png"),
    #     joinpath(output_path, "_ks_law/bigiel2010_neutral.png"),
    # ]

    # GalaxyInspector.hvcatImages(
    #     1,
    #     neutral_paths;
    #     output_path=joinpath(output_path, "_ks_law/neutral.png"),
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
    #     output_file=joinpath(output_path, "_ks_law/bigiel2010_total.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 20),),
    #         Axis=(
    #             limits=(6.4, 8.6, -4.5, 0.5),
    #             xticks=6.5:0.5:8.5,
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
    #     output_file=joinpath(output_path, "_ks_law/leroy2008_total.png"),
    #     filter_mode,
    #     sim_labels=[label],
    #     theme=Theme(
    #         Legend=(padding=(10, 0, 0, 20),),
    #         Axis=(
    #             limits=(6.4, 8.6, -4.5, 0.5),
    #             xticks=6.5:0.5:8.5,
    #             yticks=-4:1:0,
    #         ),
    #     ),
    # )

    # total_paths = [
    #     joinpath(output_path, "_ks_law/leroy2008_total.png"),
    #     joinpath(output_path, "_ks_law/bigiel2010_total.png"),
    # ]

    # GalaxyInspector.hvcatImages(
    #     1,
    #     total_paths;
    #     output_path=joinpath(output_path, "_ks_law/total.png"),
    # )

    # ###############################
    # # Final image with all KS laws
    # ###############################

    # gas_paths = joinpath.(
    #     output_path,
    #     "_ks_law",
    #     ["molecular.png", "atomic.png", "neutral.png", "total.png"],
    # )

    # GalaxyInspector.hvcatImages(
    #     4,
    #     gas_paths;
    #     output_path=joinpath(output_path, "ks_law.png"),
    # )

    # rm(joinpath(output_path, "_ks_law"); recursive=true, force=true)

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

    # temp_folder = joinpath(output_path, "_stellar_density_maps")

    # grid = GalaxyInspector.CubicGrid(R3, 400)
    # half_box_size = ustrip(u"kpc", R3) / 2.0

    # x_label = GalaxyInspector.getLabel("x", 0, u"kpc")
    # y_label = GalaxyInspector.getLabel("y", 0, u"kpc")
    # z_label = GalaxyInspector.getLabel("z", 0, u"kpc")

    # x_limits = half_box_size
    # y_limits = [half_box_size, 12]

    # filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
    #     cosmological ? :subhalo : :all_stellar,
    #     GalaxyInspector.plotParams(:stellar_mass).request,
    # )

    # for projection_plane in [:xy, :xz]

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

    #     f = Figure(size=(880, 1300), figure_padding=(5, 25, 0, 0))

    #     # Color range
    #     min_color = Inf
    #     max_color = -Inf

    #     # Compute a good color range
    #     for projection_plane in [:xy, :xz]

    #         path = joinpath(temp_folder, "stellar_mass_$(projection_plane).jld2")

    #         jldopen(path, "r") do jld2_file

    #             _, _, z = jld2_file[first(keys(jld2_file))]["simulation_001"]

    #             if !all(isnan, z)

    #                 min_Σ, max_Σ = extrema(filter(!isnan, z))

    #                 floor_Σ = floor(min_Σ)
    #                 ceil_Σ  = ceil(max_Σ)

    #                 if floor_Σ < min_color
    #                     min_color = floor_Σ
    #                 end
    #                 if ceil_Σ > max_color
    #                     max_color = ceil_Σ
    #                 end

    #             end

    #         end

    #     end

    #     # Maximun tick for the axes
    #     tick = floor(half_box_size; sigdigits=1)

    #     for (row, projection_plane) in pairs([:xy, :xz])

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
    #             xticks=-tick:10:tick,
    #             yticks=-tick:10:tick,
    #             limits=(-x_limits, x_limits, -y_limits[row], y_limits[row]),
    #             aspect=DataAspect(),
    #         )

    #         path = joinpath(temp_folder, "stellar_mass_$(projection_plane).jld2")

    #         jldopen(path, "r") do jld2_file

    #             x, y, z = jld2_file[first(keys(jld2_file))]["simulation_001"]

    #             pf = heatmap!(ax, x, y, z; colorrange=(min_color + 0.5, max_color - 0.5))

    #             if row == 1

    #                 Colorbar(
    #                     f[row, 1],
    #                     pf,
    #                     label=L"\log_{10} \, \Sigma_* \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    #                     ticklabelsize=28,
    #                     ticks=min_color:0.5:max_color,
    #                     vertical=false,
    #                 )

    #             end

    #         end

    #     end

    #     rowsize!(f.layout, 3, Relative(0.3f0))

    #     save(joinpath(output_path, "stellar_density_maps.png"), f)

    # end

    # rm(temp_folder; recursive=true)

    # ################################################################################################
    # # Gas components density map of the last snapshot, face-on/edge-on projections
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(
    #         log_file,
    #         "# Gas components density map of the last snapshot, face-on/edge-on projections",
    #     )
    #     println(log_file, "#"^100, "\n")
    # end

    # temp_folder = joinpath(output_path, "_density_maps")

    # grid = GalaxyInspector.CubicGrid(R3, 400)
    # half_box_size = ustrip(u"kpc", R3) / 2.0

    # projections = [:xy, :xz]
    # quantities = [:gas_mass, :molecular_mass, :atomic_mass, :ionized_mass, :dust_mass]

    # colorbar_labels = [
    #     L"\log_{10} \, \Sigma_\mathrm{gas} \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    #     L"\log_{10} \, \Sigma_\mathrm{H2} \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
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
    # y_size = 1020

    # jld2_paths = Vector{String}(undef, n_rows * n_cols)

    # for (i, quantity) in pairs(quantities)

    #     for (j, projection_plane) in pairs(projections)

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

    #         jld2_paths[i + n_cols * (j - 1)] = joinpath(
    #             temp_folder,
    #             "$(quantity)_$(projection_plane).jld2",
    #         )

    #     end

    # end

    # with_theme(default_theme) do

    #     f = Figure(size=(x_size, y_size), figure_padding=(5, 10, 5, 0))

    #     # Compute a good color range
    #     min_color = Inf
    #     max_color = -Inf
    #     MIN_Σ     = 3.0

    #     for path in jld2_paths

    #         jldopen(path, "r") do jld2_file

    #             _, _, z = jld2_file[first(keys(jld2_file))]["simulation_001"]

    #             if !all(isnan, z)

    #                 min_Σ, max_Σ = extrema(filter(!isnan, z))

    #                 floor_Σ = floor(min_Σ)
    #                 ceil_Σ  = ceil(max_Σ)

    #                 if floor_Σ < min_color
    #                     min_color = floor_Σ
    #                 end
    #                 if ceil_Σ > max_color
    #                     max_color = ceil_Σ
    #                 end

    #             end

    #         end

    #     end

    #     if min_color < MIN_Σ
    #         min_color = MIN_Σ
    #     end

    #     # Maximun tick for the axes
    #     tick = floor(half_box_size; sigdigits=1)

    #     for (idx, path) in enumerate(jld2_paths)

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
    #             xticks=-tick:10:tick,
    #             yticks=-tick:10:tick,
    #             limits=(-half_box_size, half_box_size, -half_box_size, half_box_size),
    #             aspect=DataAspect(),
    #         )

    #         jldopen(path, "r") do jld2_file

    #             x, y, z = jld2_file[first(keys(jld2_file))]["simulation_001"]

    #             pf = heatmap!(ax, x, y, z; colorrange=(min_color + 0.5, max_color - 0.5))

    #             if row == 1

    #                 Colorbar(
    #                     f[row, col],
    #                     pf,
    #                     label=colorbar_labels[col],
    #                     ticklabelsize=23,
    #                     ticks=min_color:1:max_color,
    #                     vertical=false,
    #                 )

    #             end

    #         end

    #         colgap!(f.layout, 40)

    #     end

    #     save(joinpath(output_path, "gas_density_maps.png"), f)

    # end

    # rm(temp_folder; recursive=true)

    # ################################################################################################
    # # Face-on density maps for different quantities of the last snapshot
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(log_file, "# Face-on density maps for different quantities of the last snapshot")
    #     println(log_file, "#"^100, "\n")
    # end

    # temp_folder = joinpath(output_path, "_face_on_maps")

    # box_size       = R3
    # pixel_length   = box_size / 400.0
    # half_box_size  = ustrip(u"kpc", box_size) / 2.0
    # tick           = floor(half_box_size; sigdigits=1)
    # size           = (910, 760)
    # figure_padding = (20, 20, 30, 20)

    # densityMap(
    #     [simulation_path],
    #     n_snapshots;
    #     quantities=[:gas_mass],
    #     types=[:cells],
    #     output_path=temp_folder,
    #     filter_mode=cosmological ? :subhalo : :all_stellar,
    #     projection_planes=[:xy],
    #     box_size,
    #     pixel_length,
    #     theme=Theme(;
    #         size,
    #         figure_padding,
    #         Colorbar=(
    #             label=L"\log_{10} \, \Sigma_\text{gas} \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    #         ),
    #         Axis=(
    #             limits=(-half_box_size, half_box_size, -half_box_size, half_box_size),
    #             xticks=-tick:10:tick,
    #             yticks=-tick:10:tick,
    #         ),
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
    #     filter_mode=cosmological ? :subhalo : :all_stellar,
    #     projection_planes=[:xy],
    #     box_size,
    #     pixel_length,
    #     theme=Theme(;
    #         size,
    #         figure_padding,
    #         Colorbar=(
    #             label=L"\log_{10} \, \Sigma_\star \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    #         ),
    #         Axis=(
    #             limits=(-half_box_size, half_box_size, -half_box_size, half_box_size),
    #             xticks=-tick:10:tick,
    #             yticks=-tick:10:tick,
    #         ),
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
    #     filter_mode=cosmological ? :subhalo : :all_stellar,
    #     projection_planes=[:xy],
    #     box_size,
    #     pixel_length,
    #     theme=Theme(;
    #         size,
    #         figure_padding,
    #         Colorbar=(label=L"\log_{10} \, T \,\, [\mathrm{K}]",),
    #         Axis=(
    #             limits=(-half_box_size, half_box_size, -half_box_size, half_box_size),
    #             xticks=-tick:10:tick,
    #             yticks=-tick:10:tick,
    #         ),
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
    #     filter_mode=cosmological ? :subhalo : :all_stellar,
    #     projection_planes=[:xy],
    #     box_size,
    #     pixel_length,
    #     theme=Theme(;
    #         size,
    #         figure_padding,
    #         Colorbar=(label=L"\log_{10} \, Z_\text{gas} \, [\mathrm{Z_\odot}]",),
    #         Axis=(
    #             limits=(-half_box_size, half_box_size, -half_box_size, half_box_size),
    #             xticks=-tick:10:tick,
    #             yticks=-tick:10:tick,
    #         ),
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
    #     filter_mode=cosmological ? :subhalo : :all_stellar,
    #     projection_planes=[:xy],
    #     box_size,
    #     pixel_length,
    #     theme=Theme(;
    #         size,
    #         figure_padding,
    #         Colorbar=(label=L"\log_{10} \, Z_\star \, [\mathrm{Z_\odot}]",),
    #         Axis=(
    #             limits=(-half_box_size, half_box_size, -half_box_size, half_box_size),
    #             xticks=-tick:10:tick,
    #             yticks=-tick:10:tick,
    #         ),
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
    #     filter_mode=cosmological ? :subhalo : :all_stellar,
    #     projection_planes=[:xy],
    #     box_size,
    #     pixel_length,
    #     theme=Theme(;
    #         size,
    #         figure_padding,
    #         Colorbar=(
    #             label=L"\log_{10} \, \text{SFR}_\text{gas} \, [\mathrm{M_\odot \, yr^{-1}}]",
    #         ),
    #         Axis=(
    #             limits=(-half_box_size, half_box_size, -half_box_size, half_box_size),
    #             xticks=-tick:10:tick,
    #             yticks=-tick:10:tick,
    #         ),
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
    #     output_path=joinpath(output_path, "face_on_maps.png"),
    # )

    # rm(temp_folder, recursive=true, force=true)

    # ################################################################################################
    # # Stellar density video, face-on/edge-on projections
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(log_file, "# Stellar density video, face-on/edge-on projections")
    #     println(log_file, "#"^100, "\n")
    # end

    # quantity = :stellar_mass
    # label    = L"\log_{10} \, \Sigma_* \,\, [\mathrm{M_\odot \, kpc^{-2}}]"

    # temp_folder = joinpath(output_path, "_$(quantity)_maps_video")

    # grid = GalaxyInspector.CubicGrid(R3, 400)
    # half_box_size = ustrip(u"kpc", R3) / 2.0

    # # Maximun tick for the axes
    # tick = floor(half_box_size; sigdigits=1)

    # x_label = GalaxyInspector.getLabel("x", 0, u"kpc")
    # y_label = GalaxyInspector.getLabel("y", 0, u"kpc")
    # z_label = GalaxyInspector.getLabel("z", 0, u"kpc")

    # x_limits = half_box_size
    # y_limits = [half_box_size, 12]

    # filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
    #     cosmological ? :subhalo : :all_stellar,
    #     GalaxyInspector.plotParams(quantity).request,
    # )

    # simulation_table = GalaxyInspector.makeSimulationTable(simulation_path)

    # # Read the last snapshot
    # last_dd = makeDataDict(
    #     simulation_path,
    #     n_snapshots,
    #     request,
    #     simulation_table,
    # )

    # # Filter the last snapshot
    # GalaxyInspector.filterData!(last_dd; filter_function)

    # # Compute the translation for the last snapshot
    # last_origin = GalaxyInspector.computeCenter(last_dd, translation)
    # last_vcm    = GalaxyInspector.computeVcm(last_dd, translation)

    # # Translate the last snapshot
    # GalaxyInspector.translateData!(last_dd, (last_origin, last_vcm))

    # # Compute the rotation for the last snapshot
    # rotation_matrix = GalaxyInspector.computeRotation(last_dd, rotation)

    # # Rotate the last snapshot
    # GalaxyInspector.rotateData!(last_dd, rotation_matrix)

    # # Compute a good color range
    # _, _, z = GalaxyInspector.daDensity2DProjection(last_dd, grid, quantity, :particles)

    # if !all(isnan, z)

    #     min_Σ, max_Σ = extrema(filter(!isnan, z))

    #     min_color = floor(min_Σ)
    #     max_color = ceil(max_Σ)

    #     colorrange = (min_color + 0.5, max_color - 0.5)
    #     colorticks = min_color:0.5:max_color

    # else

    #     colorrange = (0.0, 1.0)
    #     colorticks = 0:0.5:1

    # end

    # for projection_plane in [:xy, :xz]

    #     GalaxyInspector.plotSnapshot(
    #         [simulation_path],
    #         request,
    #         [heatmap!];
    #         output_path=temp_folder,
    #         base_filename="stellar_mass_$(projection_plane)",
    #         filter_function,
    #         da_functions=[GalaxyInspector.daDensity2DProjection],
    #         da_args=[(grid, quantity, :particles)],
    #         da_kwargs=[(; projection_plane)],
    #         transform_box=true,
    #         translation=(last_origin, last_vcm),
    #         rotation=rotation_matrix,
    #         x_unit=u"kpc",
    #         y_unit=u"kpc",
    #         save_figures=false,
    #         backup_results=true,
    #     )

    # end

    # prog_bar = Progress(
    #     nrow(simulation_table),
    #     dt=0.5,
    #     desc="Analyzing and plotting the data... ",
    #     color=:blue,
    #     barglyphs=BarGlyphs("|#  |"),
    # )

    # with_theme(default_theme) do

    #     f = Figure(size=(880, 1300), figure_padding=(5, 25, 0, 0))

    #     # Initialize the animation stream
    #     vs = VideoStream(f; framerate=20)

    #     iterator = zip(simulation_table[!, :numbers], simulation_table[!, :physical_times])

    #     for (snap_n, time) in iterator

    #         for (row, proj_plane) in pairs([:xy, :xz])

    #             xaxis_v = row == 2

    #             ax = CairoMakie.Axis(
    #                 f[row+1, 1];
    #                 xlabel=x_label,
    #                 ylabel=(row == 1 ? y_label : z_label),
    #                 xminorticksvisible=xaxis_v,
    #                 xticksvisible=xaxis_v,
    #                 xlabelvisible=xaxis_v,
    #                 xticklabelsvisible=xaxis_v,
    #                 xticklabelsize=35,
    #                 yticklabelsize=35,
    #                 xticks=-tick:10:tick,
    #                 yticks=-tick:10:tick,
    #                 limits=(-x_limits, x_limits, -y_limits[row], y_limits[row]),
    #                 aspect=DataAspect(),
    #             )

    #             address = "$(quantity)_$(proj_plane)_$(GalaxyInspector.SNAP_BASENAME)_$(snap_n)"
    #             path    = joinpath(temp_folder, "$(quantity)_$(proj_plane).jld2")

    #             jldopen(path, "r") do jld2_file

    #                 x, y, z = jld2_file[address * "/simulation_001"]

    #                 pf = heatmap!(ax, x, y, z; colorrange)

    #                 if row == 1

    #                     Colorbar(
    #                         f[row, 1],
    #                         pf;
    #                         label,
    #                         ticklabelsize=28,
    #                         ticks=colorticks,
    #                         vertical=false,
    #                     )

    #                     c_t = ustrip(u"Gyr", time)

    #                     if c_t < 1.0
    #                         time_stamp = round(c_t; digits=2)
    #                     else
    #                         time_stamp = round(c_t; sigdigits=3)
    #                     end

    #                     text!(
    #                         ax,
    #                         0.73,
    #                         0.98;
    #                         text=L"t = %$(rpad(time_stamp, 4, '0')) \, \text{Gyr}",
    #                         align=(:left, :top),
    #                         color=:white,
    #                         space=:relative,
    #                         fontsize=35,
    #                     )

    #                 end

    #             end

    #         end

    #         rowsize!(f.layout, 3, Relative(0.3f0))

    #         # Add the figure as a frame to the animation stream
    #         recordframe!(vs)

    #         GalaxyInspector.cleanPlot!(f)

    #         next!(prog_bar)

    #     end

    #     save(joinpath(output_path, "$(quantity)_map.mkv"), vs)

    # end

    # rm(temp_folder; recursive=true)

    # ################################################################################################
    # # Stellar and gas density video, face-on/edge-on projections
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(log_file, "# Stellar and gas density video, face-on/edge-on projections")
    #     println(log_file, "#"^100, "\n")
    # end

    # temp_folder = joinpath(output_path, "_density_maps_video")

    # quantities = [:gas_mass, :molecular_mass, :stellar_mass]
    # labels = [
    #     L"\log_{10} \, \Sigma_\text{gas} \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    #     L"\log_{10} \, \Sigma_\text{H2} \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    #     L"\log_{10} \, \Sigma_* \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
    # ]
    # types = [:cells, :particles, :particles]

    # grid = GalaxyInspector.CubicGrid(R3, 400)
    # half_box_size = ustrip(u"kpc", R3) / 2.0

    # # Maximun tick for the axes
    # tick = floor(half_box_size; sigdigits=1)

    # x_label = GalaxyInspector.getLabel("x", 0, u"kpc")
    # y_label = GalaxyInspector.getLabel("y", 0, u"kpc")
    # z_label = GalaxyInspector.getLabel("z", 0, u"kpc")

    # x_limits = half_box_size
    # y_limits = [half_box_size, 12]

    # requests = [GalaxyInspector.plotParams(quantity).request for quantity in quantities]

    # filter_function, translation, rotation, request = GalaxyInspector.selectFilter(
    #     cosmological ? :subhalo : :all_stellar,
    #     GalaxyInspector.mergeRequests(requests...),
    # )

    # simulation_table = GalaxyInspector.makeSimulationTable(simulation_path)

    # #############################
    # # Compute a good color range
    # #############################

    # # Read the last snapshot
    # last_dd = makeDataDict(
    #     simulation_path,
    #     n_snapshots,
    #     request,
    #     simulation_table,
    # )

    # # Filter the last snapshot
    # GalaxyInspector.filterData!(last_dd; filter_function)

    # # Compute the translation for the last snapshot
    # last_origin = GalaxyInspector.computeCenter(last_dd, translation)
    # last_vcm    = GalaxyInspector.computeVcm(last_dd, translation)

    # # Translate the last snapshot
    # GalaxyInspector.translateData!(last_dd, (last_origin, last_vcm))

    # # Compute the rotation for the last snapshot
    # rotation_matrix = GalaxyInspector.computeRotation(last_dd, rotation)

    # # Rotate the last snapshot
    # GalaxyInspector.rotateData!(last_dd, rotation_matrix)

    # ######################
    # # Stellar color range
    # ######################

    # _, _, z = GalaxyInspector.daDensity2DProjection(last_dd, grid, :stellar_mass, :particles)

    # if !all(isnan, z)

    #     min_Σ, max_Σ = extrema(filter(!isnan, z))

    #     min_color = floor(min_Σ)
    #     max_color = ceil(max_Σ)

    #     stellar_colorrange = (min_color + 0.5, max_color - 0.5)
    #     stellar_colorticks = min_color:0.5:max_color

    # else

    #     stellar_colorrange = (0.0, 1.0)
    #     stellar_colorticks = 0:0.5:1

    # end

    # ##################
    # # Gas color range
    # ##################

    # MIN_Σ = 4.0

    # _, _, z = GalaxyInspector.daDensity2DProjection(last_dd, grid, :molecular_mass, :cells)

    # if !all(isnan, z)

    #     min_Σ, max_Σ = extrema(filter(!isnan, z))

    #     min_color = floor(min_Σ)
    #     max_color = ceil(max_Σ)

    #     if min_color < MIN_Σ
    #         min_color = MIN_Σ
    #     end

    #     gas_colorrange = (min_color - 0.5, max_color + 0.5)
    #     gas_colorticks = (min_color - 0.5):1.0:(max_color + 0.5)

    # else

    #     gas_colorrange = (0.0, 1.0)
    #     gas_colorticks = 0:0.5:1

    # end

    # colorranges = [gas_colorrange, gas_colorrange, stellar_colorrange]
    # colortickss = [gas_colorticks, gas_colorticks, stellar_colorticks]

    # for (quantity, type) in zip(quantities, types)

    #     for projection_plane in [:xy, :xz]

    #         GalaxyInspector.plotSnapshot(
    #             [simulation_path],
    #             request,
    #             [heatmap!];
    #             output_path=temp_folder,
    #             base_filename="$(quantity)_$(projection_plane)",
    #             filter_function,
    #             da_functions=[GalaxyInspector.daDensity2DProjection],
    #             da_args=[(grid, quantity, type)],
    #             da_kwargs=[(; projection_plane)],
    #             transform_box=true,
    #             translation=(last_origin, last_vcm),
    #             rotation=rotation_matrix,
    #             x_unit=u"kpc",
    #             y_unit=u"kpc",
    #             save_figures=false,
    #             backup_results=true,
    #         )

    #     end

    # end

    # prog_bar = Progress(
    #     nrow(simulation_table),
    #     dt=0.5,
    #     desc="Analyzing and plotting the data... ",
    #     color=:blue,
    #     barglyphs=BarGlyphs("|#  |"),
    # )

    # with_theme(default_theme) do

    #     f = Figure(size=(2400, 1300), figure_padding=(5, 25, 0, 0))

    #     # Initialize the animation stream
    #     vs = VideoStream(f; framerate=20)

    #     snap_iterator = zip(simulation_table[!, :numbers], simulation_table[!, :physical_times])
    #     col_iterator = enumerate(zip(quantities, labels, colorranges, colortickss))

    #     for (snap_n, time) in snap_iterator

    #         for (col, (quantity, label, colorrange, colorticks)) in col_iterator

    #             for (row, proj_plane) in pairs([:xy, :xz])

    #                 xaxis_v = row == 2
    #                 yaxis_v = col == 1

    #                 ax = CairoMakie.Axis(
    #                     f[row+1, col];
    #                     xlabel=x_label,
    #                     ylabel=(row == 1 ? y_label : z_label),
    #                     xminorticksvisible=xaxis_v,
    #                     xticksvisible=xaxis_v,
    #                     xlabelvisible=xaxis_v,
    #                     xticklabelsvisible=xaxis_v,
    #                     yminorticksvisible=yaxis_v,
    #                     yticksvisible=yaxis_v,
    #                     ylabelvisible=yaxis_v,
    #                     yticklabelsvisible=yaxis_v,
    #                     xticklabelsize=35,
    #                     yticklabelsize=35,
    #                     xticks=-tick:10:tick,
    #                     yticks=-tick:10:tick,
    #                     limits=(-x_limits, x_limits, -y_limits[row], y_limits[row]),
    #                     aspect=DataAspect(),
    #                 )

    #                 address = "$(quantity)_$(proj_plane)_$(GalaxyInspector.SNAP_BASENAME)_$(snap_n)"
    #                 path = joinpath(temp_folder, "$(quantity)_$(proj_plane).jld2")

    #                 jldopen(path, "r") do jld2_file

    #                     x, y, z = jld2_file[address * "/simulation_001"]

    #                     pf = heatmap!(ax, x, y, z; colorrange)

    #                     if row == 1

    #                         Colorbar(
    #                             f[row, col],
    #                             pf;
    #                             label,
    #                             ticklabelsize=28,
    #                             ticks=colorticks,
    #                             vertical=false,
    #                         )

    #                         c_t = ustrip(u"Gyr", time)

    #                         if c_t < 1.0
    #                             time_stamp = round(c_t; digits=2)
    #                         else
    #                             time_stamp = round(c_t; sigdigits=3)
    #                         end

    #                         text!(
    #                             ax,
    #                             0.68,
    #                             0.98;
    #                             text=L"t = %$(rpad(time_stamp, 4, '0')) \, \text{Gyr}",
    #                             align=(:left, :top),
    #                             color=:white,
    #                             space=:relative,
    #                             fontsize=40,
    #                         )

    #                     end

    #                 end

    #             end

    #             rowsize!(f.layout, 3, Relative(0.3f0))

    #         end

    #         colgap!(f.layout, 70)

    #         # Add the figure as a frame to the animation stream
    #         recordframe!(vs)

    #         GalaxyInspector.cleanPlot!(f)

    #         next!(prog_bar)

    #     end

    #     save(joinpath(output_path, "density_maps.mkv"), vs)

    # end

    # # rm(temp_folder; recursive=true)

    # ################################################################################################
    # # Evolution of the masses and of the fractions, for the different gas components
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(
    #         log_file,
    #         "# Evolution of the masses and of the fractions, for the different gas components",
    #     )
    #     println(log_file, "#"^100, "\n")
    # end

    # quantities = [
    #     :ionized_fraction,
    #     :atomic_fraction,
    #     :molecular_fraction,
    #     :dust_fraction,
    #     :metal_gas_fraction,
    # ]

    # x_plot_params = GalaxyInspector.plotParams(:physical_time)
    # y_plot_params = GalaxyInspector.plotParams(:generic_fraction)

    # temp_folder = joinpath(output_path, "_gas_evolution")

    # plotTimeSeries(
    #     fill(simulation_path, length(quantities)),
    #     [lines!];
    #     output_path=temp_folder,
    #     filename="gas_fractions_evolution",
    #     da_functions=[GalaxyInspector.daEvolution],
    #     da_args=[(:physical_time, quantity) for quantity in quantities],
    #     da_kwargs=[(;
    #         filter_mode,
    #         extra_filter=dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", R1), :zero),
    #         scaling=identity,
    #     )],
    #     x_unit=x_plot_params.unit,
    #     y_unit=y_plot_params.unit,
    #     x_exp_factor=x_plot_params.exp_factor,
    #     y_exp_factor=y_plot_params.exp_factor,
    #     xaxis_label=x_plot_params.axis_label,
    #     yaxis_label=y_plot_params.axis_label,
    #     xaxis_var_name=x_plot_params.var_name,
    #     yaxis_var_name=y_plot_params.var_name,
    #     save_figures=false,
    #     backup_results=true,
    # )

    # quantities = [
    #     :stellar_mass,
    #     :ionized_mass,
    #     :atomic_mass,
    #     :molecular_mass,
    #     :dust_mass,
    #     :ode_metal_mass,
    # ]

    # y_plot_params = GalaxyInspector.plotParams(:generic_mass)

    # plotTimeSeries(
    #     fill(simulation_path, length(quantities)),
    #     [lines!];
    #     output_path=temp_folder,
    #     filename="gas_masses_evolution",
    #     da_functions=[GalaxyInspector.daEvolution],
    #     da_args=[(:physical_time, quantity) for quantity in quantities],
    #     da_kwargs=[(;
    #         filter_mode,
    #         extra_filter=dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", R1), :zero),
    #         scaling=identity,
    #     )],
    #     x_unit=x_plot_params.unit,
    #     y_unit=y_plot_params.unit,
    #     x_exp_factor=x_plot_params.exp_factor,
    #     y_exp_factor=y_plot_params.exp_factor,
    #     xaxis_label=x_plot_params.axis_label,
    #     yaxis_label=y_plot_params.axis_label,
    #     xaxis_var_name=x_plot_params.var_name,
    #     yaxis_var_name=y_plot_params.var_name,
    #     save_figures=false,
    #     backup_results=true,
    # )

    # paths = joinpath.(temp_folder, ["gas_masses_evolution.jld2", "gas_fractions_evolution.jld2"])

    # current_theme = merge(Theme(palette=(linestyle=[:solid],),), default_theme)

    # with_theme(current_theme) do

    #     f = Figure(size=(880, 1650), figure_padding=(5, 10, 5, 10))

    #     ax_1 = CairoMakie.Axis(
    #         f[1, 1];
    #         xlabel=L"t \, [\mathrm{Gyr}]",
    #         ylabel=L"\log_{10} \, M \, [\mathrm{10^{10} \, M_\odot}]",
    #         # xminorticksvisible=false,
    #         # xticksvisible=false,
    #         xlabelvisible=false,
    #         xticklabelsvisible=false,
    #         limits=(-0.1, nothing, -4.2, 1.2),
    #         yticks=-4:1:1,
    #     )

    #     jldopen(paths[1], "r") do jld2_file

    #         first_address = first(keys(jld2_file))

    #         x_s, y_s = jld2_file[first_address]["simulation_001"]
    #         x_i, y_i = jld2_file[first_address]["simulation_002"]
    #         x_a, y_a = jld2_file[first_address]["simulation_003"]
    #         x_m, y_m = jld2_file[first_address]["simulation_004"]
    #         x_d, y_d = jld2_file[first_address]["simulation_005"]
    #         x_Z, y_Z = jld2_file[first_address]["simulation_006"]

    #         lines!(ax_1, x_s, log10.(y_s); color=Makie.wong_colors()[2], label="Stellar mass")
    #         lines!(ax_1, x_i, log10.(y_i); color=Makie.wong_colors()[1], label="Ionized mass" )
    #         lines!(ax_1, x_a, log10.(y_a); color=Makie.wong_colors()[4], label="Atomic mass")
    #         lines!(ax_1, x_m, log10.(y_m); color=Makie.wong_colors()[3], label="Molecular mass")
    #         lines!(ax_1, x_d, log10.(y_d); color=Makie.wong_colors()[6], label="Dust mass")
    #         lines!(ax_1, x_Z, log10.(y_Z); color=Makie.wong_colors()[5], label="Metals mass")

    #         # z ~ 6.0 (t = 0.95 Gyr)
    #         vlines!(ax_1, 0.95, color=:gray25)

    #         # 10^7 M⊙
    #         hlines!(ax_1, -3.0, color=:gray25)

    #     end

    #     axislegend(ax_1, position=:rb, framevisible=false, nbanks=2)

    #     ax_2 = CairoMakie.Axis(
    #         f[2, 1];
    #         xlabel=L"t \, [\mathrm{Gyr}]",
    #         ylabel=L"\log_{10} \, f",
    #         aspect=nothing,
    #         limits=(-0.1, nothing, -5.2, 0.2),
    #         yticks=-5:1:0,
    #     )

    #     jldopen(paths[2], "r") do jld2_file

    #         first_address = first(keys(jld2_file))

    #         x_i, y_i = jld2_file[first_address]["simulation_001"]
    #         x_a, y_a = jld2_file[first_address]["simulation_002"]
    #         x_m, y_m = jld2_file[first_address]["simulation_003"]
    #         x_d, y_d = jld2_file[first_address]["simulation_004"]
    #         x_Z, y_Z = jld2_file[first_address]["simulation_005"]

    #         lines!(ax_2, x_i, log10.(y_i); color=Makie.wong_colors()[1])
    #         lines!(ax_2, x_a, log10.(y_a); color=Makie.wong_colors()[4])
    #         lines!(ax_2, x_m, log10.(y_m); color=Makie.wong_colors()[3])
    #         lines!(ax_2, x_d, log10.(y_d); color=Makie.wong_colors()[6])
    #         lines!(ax_2, x_Z, log10.(y_Z); color=Makie.wong_colors()[5])

    #     end

    #     linkxaxes!(ax_1, ax_2)

    #     save(joinpath(output_path, "gas_evolution.png"), f)

    # end

    # rm(temp_folder, recursive=true, force=true)

    ################################################################################################
    # Dust to stellar mass ratio
    ################################################################################################

    y_plot_params = GalaxyInspector.plotParams(:dust_stellar_mass_ratio)
    x_plot_params = GalaxyInspector.plotParams(:physical_time)

    plotTimeSeries(
        [simulation_path],
        [scatter!];
        output_path,
        filename="dust_to_stellar_mass_ratio",
        da_functions=[GalaxyInspector.daEvolution],
        da_args=[(:physical_time, :dust_stellar_mass_ratio)],
        da_kwargs=[(;
            filter_mode,
            extra_filter=dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", R1), :zero),
            scaling=identity,
        )],
        post_processing=GalaxyInspector.ppHorizontalFlags!,
        pp_args=([exp10(-2.0), exp10(-3.0), exp10(-4.0)],),
        pp_kwargs=(; colors=[:gray65], line_styles=[:dash]),
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=y_plot_params.exp_factor,
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        yaxis_scale_func=log10,
    )

    ##############
    # Close files
    ##############

    logging && close(log_file)

    return nothing

end

function comparison(
    simulation_paths::Vector{String},
    output_path::String,
    logging::Bool;
    labels::Vector{<:AbstractString}=basename.(simulation_paths),
)::Nothing

    # Number of simulations
    n_sims = length(simulation_paths)
    @assert n_sims == length(labels) "Number of simulations and labels must match."

    # Create the output folder
    mkpath(output_path)

    # Activate logging
    if logging
        log_file = open(joinpath(output_path, "logs.txt"), "w+")
        GalaxyInspector.setLogging!(logging; stream=log_file)
    end

    temp_folder = joinpath(output_path, "_mass_evolution")

    # Number of snapshots
    n_snapshots = minimum(GalaxyInspector.countSnapshot.(simulation_paths))

    # Characteristic radii
    R1 = 40.0u"kpc"

    # Check if the simulations are cosmological
    cosmological = all(GalaxyInspector.isSimCosmological, simulation_paths)
    if cosmological
        filter_mode  = :subhalo
        extra_filter = dd -> GalaxyInspector.filterWithinSphere(dd, (0.0u"kpc", R1), :zero)
    else
        filter_mode  = :all
        extra_filter = GalaxyInspector.filterNothing
    end

    ################################################################################################
    # Evolution of the gas components and the SFR
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# Evolution of the gas components and the SFR")
        println(log_file, "#"^100, "\n")
    end

    quantities   = [:ionized_mass, :atomic_mass, :molecular_mass, :dust_mass, :sfr, :stellar_mass]
    smooths      = [0, 0, 0, 0, n_snapshots ÷ 5, 0]
    n_quantities = length(quantities)

    jld2_paths = Vector{String}(undef, n_quantities * n_sims)

    for (i, simulation) in pairs(simulation_paths)

        sim_name = basename(simulation)

        for (j, (quantity, smooth)) in enumerate(zip(quantities, smooths))

            if logging
                println(log_file, "#"^100)
                println(
                    log_file,
                    "# Computing the evolution of $(quantity) for $(sim_name)",
                )
                println(log_file, "#"^100, "\n")
            end

            y_plot_params = GalaxyInspector.plotParams(quantity)

            plotTimeSeries(
                [simulation],
                [lines!];
                output_path=temp_folder,
                filename="$(quantity)_$(sim_name)",
                da_functions=[GalaxyInspector.daEvolution],
                da_args=[(:physical_time, quantity)],
                da_kwargs=[(; filter_mode, extra_filter, smooth, scaling=identity)],
                x_unit=u"Gyr",
                y_unit=y_plot_params.unit,
                y_exp_factor=y_plot_params.exp_factor,
                save_figures=false,
                backup_results=true,
            )

            jld2_paths[(i - 1) * n_quantities + j] = joinpath(
                temp_folder,
                "$(quantity)_$(sim_name).jld2",
            )

        end

    end

    current_theme = merge(
        Theme(palette=(linestyle=[:solid],),),
        GalaxyInspector.DEFAULT_THEME,
        theme_latexfonts(),
    )

    colors = current_theme[:palette][:color][]

    with_theme(current_theme) do

        f = Figure(size=(2750, 1300), figure_padding=(5, 10, 5, 10))

        for (i, quantity) in pairs(quantities)

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
                f[ceil(Int, i / 3), mod1(i, 3)];
                xlabel=L"t \, [\mathrm{Gyr}]",
                ylabel,
                aspect=AxisAspect(1.4),
            )

            for (simulation_path, label, color) in zip(simulation_paths, labels, colors)

                jld2_path = joinpath(temp_folder, "$(quantity)_$(basename(simulation_path)).jld2")

                jldopen(jld2_path, "r") do jld2_file

                    address = first(keys(jld2_file))

                    x, y = jld2_file[address]["simulation_001"]

                    lines!(ax, x, y; label, color, linewidth=3)

                end

            end

            axislegend(ax, position=:rb, framevisible=false, nbanks=1)

        end

        save(joinpath(output_path, "mass_evolution.png"), f)

    end

    rm(temp_folder; recursive=true)

    ################################################################################################
    # Profiles for the last snapshot, comparison with Mollá et al. (2015)
    ################################################################################################

    for quantity in [
        :stellar_area_density,
        :molecular_area_density,
        :sfr_area_density,
        :atomic_area_density,
        :O_stellar_abundance,
        :N_stellar_abundance,
        :C_stellar_abundance,
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
            simulation_paths,
            n_snapshots,
            quantity;
            output_path=joinpath(output_path, "Molla2015"),
            filter_mode,
            sim_labels=labels,
            theme=Theme(
                size=(1500, 880),
                figure_padding=(10, 15, 5, 15),
                palette=(linestyle=[:solid],),
                Axis=(aspect=nothing,),
                Legend=(halign=:right, valign=:top, nbanks=1),
            ),
        )

    end

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

    SIMULATIONS =
    vcat(
        joinpath.(
            "F:/simulations/isolated/",
            [
                # "SFM_064",
                # "STD_064",
                # "BLT_064",
                # "SFM_032",
                "SFM_064_T00",
                # "SFM_064_T01",
                # "SFM_064_T02",
                # "SFM_064_T03",
                # "SFM_064_T04",
                # "SFM_064_T05",
                # "SFM_064_T06",
                # "SFM_064_T07",
                # "SFM_064_T08",
                # "SFM_128",
            ]
        ),
        ["F:/simulations/current/test_dust_11"],
    )

    for simulation in SIMULATIONS
        basic_analysis(simulation, BASE_OUT_PATH, LOGGING)
    end

    # ##################################################
    # # SF model tests
    # ##################################################

    # comparison(
    #     joinpath.("F:/simulations/isolated/", ["SFM_064", "STD_064", "BLT_064"]),
    #     joinpath(BASE_OUT_PATH, "comparison/sf_model"),
    #     LOGGING;
    #     labels=["SFM", "STD", "BLT"],
    # )

    # ##################################################
    # # LW background tests
    # ##################################################

    # comparison(
    #     joinpath.("F:/simulations/isolated/", ["SFM_064_T00", "SFM_064_T01", "SFM_064_T02"]),
    #     joinpath(BASE_OUT_PATH, "comparison/lw_background"),
    #     LOGGING;
    #     labels=["Incatasciato et al. (2023)", "Ahn et al. (2009)", "Visbal et al. (2014)"],
    # )

    # ##################################################
    # # Model equation tests
    # ##################################################

    # comparison(
    #     joinpath.(
    #         "F:/simulations/isolated/",
    #         ["SFM_064_T03", "SFM_064_T04", "SFM_064_T05", "SFM_064_T06", "SFM_064"]
    #     ),
    #     joinpath(BASE_OUT_PATH, "comparison/model_equation"),
    #     LOGGING;
    #     labels=[
    #         "Paper 1",
    #         "Paper 1 + dust",
    #         "Paper 1 + dust + shielding",
    #         "Paper 1 + dust + shielding + LW",
    #         "Paper 2",
    #     ],
    # )

    # ##################################################
    # # Resolution tests
    # ##################################################

    # comparison(
    #     joinpath.("F:/simulations/isolated/", ["SFM_032", "SFM_064", "SFM_128"]),
    #     joinpath(BASE_OUT_PATH, "comparison/resolution"),
    #     LOGGING;
    #     labels=[L"32^3", L"64^3", L"128^3"],
    # )

    # ##################################################
    # # Isolated IC test
    # ##################################################

    # comparison(
    #     joinpath.("F:/simulations/isolated/", ["SFM_064", "SFM_064_T08"]),
    #     joinpath(BASE_OUT_PATH, "comparison/isolated_ic"),
    #     LOGGING;
    #     labels=["OMEGA_MAT_1.5", "OMEGA_MAT_1.0"],
    # )

end
