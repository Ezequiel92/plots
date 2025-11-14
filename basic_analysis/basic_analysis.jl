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
        trans_mode   = :stellar_subhalo
        filter_mode  = :subhalo
        extra_filter = dd -> GalaxyInspector.filterBySphere(dd, (0.0u"kpc", R1), :zero)
    else
        trans_mode   = :stellar_box
        filter_mode  = :sphere
        extra_filter = GalaxyInspector.filterNothing
    end

    # Check if the simulation has our SF model
    sfm = GalaxyInspector.isSimSFM(simulation_path)

    # Characteristic radii
    R1 = 40.0u"kpc"
    R2 = 2.0u"kpc"
    R3 = 65.0u"kpc"

    ################################################################################################
    # Report files
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# Report files")
        println(log_file, "#"^100, "\n")
    end

    simulationReport([simulation_path]; output_path)

    snapshotReport(
        [simulation_path],
        [n_snapshots];
        output_path,
        trans_mode,
        filter_mode,
    )

    ################################################################################################
    # SFR vs. physical time
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# SFR vs. physical time")
        println(log_file, "#"^100, "\n")
    end

    timeSeries(
        [simulation_path],
        :physical_time,
        :sfr;
        xlog=false,
        ylog=false,
        output_path,
        trans_mode,
        filter_mode,
        smooth=n_snapshots > 300 ? n_snapshots ÷ 4 : 0,
        theme=Theme(
            palette=(color=[GalaxyInspector.WONG_ORANGE],),
            size=(1320, 880),
            Axis=(aspect=nothing, xticks=0:14),
        ),
    )

    ################################################################################################
    # Stellar mass vs. physical time
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# Stellar mass vs. physical time")
        println(log_file, "#"^100, "\n")
    end

    timeSeries(
        [simulation_path],
        :physical_time,
        :stellar_mass;
        xlog=false,
        ylog=false,
        output_path,
        trans_mode,
        filter_mode,
        theme=Theme(
            palette=(color=[GalaxyInspector.WONG_ORANGE],),
            size=(1320, 880),
            Axis=(aspect=nothing, xticks=0:14),
        ),
    )

    ################################################################################################
    # Stellar surface density profile of the last snapshot, comparison with Agertz et al. (2021)
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(
            log_file,
            "# Stellar surface density profile of the last snapshot, comparison with Agertz et al. \
            (2021)",
        )
        println(log_file, "#"^100, "\n")
    end

    compareAgertz2021(
        [simulation_path],
        n_snapshots;
        output_path,
        trans_mode,
        filter_mode,
        sim_labels=[label],
    )

    ################################################################################################
    # Circularity histogram of the last snapshot, radio separated
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# Circularity histogram of the last snapshot, radio separated")
        println(log_file, "#"^100, "\n")
    end

    circularityHistogram(
        [simulation_path],
        n_snapshots;
        R_in=R2,
        R_out=R1,
        output_path,
        trans_mode,
        filter_mode,
    )

    ################################################################################################
    # Efficiency per free-fall time histograms of the last snapshot
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# Efficiency per free-fall time histograms of the last snapshot")
        println(log_file, "#"^100, "\n")
    end

    efficiencyHistogram(
        [simulation_path],
        n_snapshots;
        range=(1.0e-4, 1.0),
        output_path,
        filename="eff_histogram_young_stars",
        trans_mode,
        filter_mode,
        stellar_ff=GalaxyInspector.filterByStellarAge,
        ff_request=Dict(:stellar => ["GAGE"]),
        labels=["Young stars", "Gas"]
    )

    ################################################################################################
    # Profiles for the last snapshot, comparison with Mollá et al. (2015)
    ################################################################################################

    if sfm

        molla_quantities = [
            :stellar_area_density,
            :ode_molecular_area_density,
            :sfr_area_density,
            :ode_atomic_area_density,
            :O_stellar_abundance,
            :N_stellar_abundance,
            :C_stellar_abundance,
        ]

    else

        molla_quantities = [
            :stellar_area_density,
            :br_molecular_area_density,
            :sfr_area_density,
            :br_atomic_area_density,
            :O_stellar_abundance,
            :N_stellar_abundance,
            :C_stellar_abundance,
        ]

    end

    for quantity in molla_quantities

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
            output_path=joinpath(output_path, "Molla2015"),
            trans_mode,
            filter_mode,
            sim_labels=[label],
            theme=Theme(Legend=(halign=:right, valign=:top, nbanks=1)),
        )

    end

    ################################################################################################
    # Resolved Kennicutt–Schmidt law of the last snapshot, scatter with square grid
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(
            log_file,
            "# Resolved Kennicutt–Schmidt law of the last snapshot, scatter with square grid",
        )
        println(log_file, "#"^100, "\n")
    end

    ###################
    # Molecular KS law
    ###################

    if logging
        println(log_file, "#"^50)
        println(log_file, "# Molecular KS law")
        println(log_file, "#"^50, "\n")
    end

    kennicuttSchmidtLaw(
        [simulation_path],
        n_snapshots;
        quantity=:molecular,
        grid_size=30.0u"kpc",
        bin_size=1.5u"kpc",
        post_processing=GalaxyInspector.ppSun2023!,
        pp_kwargs=(; color=GalaxyInspector.WONG_ORANGE),
        fit=true,
        output_file=joinpath(output_path, "_ks_law/sun2023_molecular.png"),
        trans_mode,
        filter_mode,
        sim_labels=[label],
        theme=Theme(
            Legend=(margin=(10, 0, 0, 0),),
            Axis=(
                limits=(-2.5, 3.5, -4.5, 0.5),
                xticks=-1:1:3,
                yticks=-4:1:0,
            ),
        ),
    )

    kennicuttSchmidtLaw(
        [simulation_path],
        n_snapshots;
        quantity=:molecular,
        grid_size=30.0u"kpc",
        bin_size=0.8u"kpc",
        post_processing=GalaxyInspector.ppLeroy2008!,
        pp_kwargs=(; color=GalaxyInspector.WONG_ORANGE),
        fit=true,
        output_file=joinpath(output_path, "_ks_law/leroy2008_molecular.png"),
        trans_mode,
        filter_mode,
        sim_labels=[label],
        theme=Theme(
            Legend=(margin=(10, 0, 0, 0),),
            Axis=(
                limits=(-2.5, 3.5, -4.5, 0.5),
                xticks=-1:1:3,
                yticks=-4:1:0,
            ),
        ),
    )

    molecular_paths = [
        joinpath(output_path, "_ks_law/leroy2008_molecular.png"),
        joinpath(output_path, "_ks_law/sun2023_molecular.png"),
    ]

    GalaxyInspector.hvcatImages(
        1,
        molecular_paths;
        output_path=joinpath(output_path, "_ks_law/molecular.png"),
    )

    ################
    # Atomic KS law
    ################

    if logging
        println(log_file, "#"^50)
        println(log_file, "# Atomic KS law")
        println(log_file, "#"^50, "\n")
    end

    kennicuttSchmidtLaw(
        [simulation_path],
        n_snapshots;
        quantity=:atomic,
        grid_size=30.0u"kpc",
        bin_size=0.6u"kpc",
        post_processing=GalaxyInspector.ppBigiel2010!,
        pp_kwargs=(; galaxy=:all, quantity=:atomic, color=GalaxyInspector.WONG_ORANGE),
        fit=false,
        output_file=joinpath(output_path, "_ks_law/bigiel2010_atomic.png"),
        trans_mode,
        filter_mode,
        sim_labels=[label],
        theme=Theme(
            Legend=(margin=(10, 0, 0, 20),),
            Axis=(
                limits=(-1.8, 2.2, -4.5, 0.5),
                xticks=0:0.5:2,
                yticks=-4:1:0,
            ),
        ),
    )

    kennicuttSchmidtLaw(
        [simulation_path],
        n_snapshots;
        quantity=:atomic,
        grid_size=30.0u"kpc",
        bin_size=0.8u"kpc",
        post_processing=GalaxyInspector.ppLeroy2008!,
        pp_kwargs=(; quantity=:atomic, color=GalaxyInspector.WONG_ORANGE),
        fit=false,
        output_file=joinpath(output_path, "_ks_law/leroy2008_atomic.png"),
        trans_mode,
        filter_mode,
        sim_labels=[label],
        theme=Theme(
            Legend=(margin=(10, 0, 0, 20),),
            Axis=(
                limits=(-1.8, 2.2, -4.5, 0.5),
                xticks=0:0.5:2,
                yticks=-4:1:0,
            ),
        ),
    )

    atomic_paths = [
        joinpath(output_path, "_ks_law/leroy2008_atomic.png"),
        joinpath(output_path, "_ks_law/bigiel2010_atomic.png"),
    ]

    GalaxyInspector.hvcatImages(
        1,
        atomic_paths;
        output_path=joinpath(output_path, "_ks_law/atomic.png"),
    )

    #################
    # Neutral KS law
    #################

    if logging
        println(log_file, "#"^50)
        println(log_file, "# Neutral KS law")
        println(log_file, "#"^50, "\n")
    end

    kennicuttSchmidtLaw(
        [simulation_path],
        n_snapshots;
        quantity=:neutral,
        grid_size=30.0u"kpc",
        bin_size=0.6u"kpc",
        post_processing=GalaxyInspector.ppBigiel2010!,
        pp_kwargs=(; galaxy=:all, quantity=:neutral, color=GalaxyInspector.WONG_ORANGE),
        fit=false,
        output_file=joinpath(output_path, "_ks_law/bigiel2010_neutral.png"),
        trans_mode,
        filter_mode,
        sim_labels=[label],
        theme=Theme(
            Legend=(margin=(10, 0, 0, 20),),
            Axis=(
                limits=(0.4, 2.6, -4.5, 0.5),
                xticks=0.5:0.5:2.5,
                yticks=-4:1:0,
            ),
        ),
    )

    kennicuttSchmidtLaw(
        [simulation_path],
        n_snapshots;
        quantity=:neutral,
        grid_size=30.0u"kpc",
        bin_size=0.8u"kpc",
        post_processing=GalaxyInspector.ppLeroy2008!,
        pp_kwargs=(; quantity=:neutral, color=GalaxyInspector.WONG_ORANGE),
        fit=false,
        output_file=joinpath(output_path, "_ks_law/leroy2008_neutral.png"),
        trans_mode,
        filter_mode,
        sim_labels=[label],
        theme=Theme(
            Legend=(margin=(10, 0, 0, 20),),
            Axis=(
                limits=(0.4, 2.6, -4.5, 0.5),
                xticks=0.5:0.5:2.5,
                yticks=-4:1:0,
            ),
        ),
    )

    neutral_paths = [
        joinpath(output_path, "_ks_law/leroy2008_neutral.png"),
        joinpath(output_path, "_ks_law/bigiel2010_neutral.png"),
    ]

    GalaxyInspector.hvcatImages(
        1,
        neutral_paths;
        output_path=joinpath(output_path, "_ks_law/neutral.png"),
    )

    ###################
    # Total gas KS law
    ###################

    if logging
        println(log_file, "#"^50)
        println(log_file, "# Total gas KS law")
        println(log_file, "#"^50, "\n")
    end

    kennicuttSchmidtLaw(
        [simulation_path],
        n_snapshots;
        quantity=:gas,
        grid_size=30.0u"kpc",
        bin_size=0.6u"kpc",
        post_processing=GalaxyInspector.ppBigiel2010!,
        pp_kwargs=(; galaxy=:all, quantity=:neutral, color=GalaxyInspector.WONG_ORANGE),
        fit=false,
        output_file=joinpath(output_path, "_ks_law/bigiel2010_total.png"),
        trans_mode,
        filter_mode,
        sim_labels=[label],
        theme=Theme(
            Legend=(margin=(10, 0, 0, 20),),
            Axis=(
                limits=(0.4, 2.6, -4.5, 0.5),
                xticks=0.5:0.5:2.5,
                yticks=-4:1:0,
            ),
        ),
    )

    kennicuttSchmidtLaw(
        [simulation_path],
        n_snapshots;
        quantity=:gas,
        grid_size=30.0u"kpc",
        bin_size=0.8u"kpc",
        post_processing=GalaxyInspector.ppLeroy2008!,
        pp_kwargs=(; quantity=:neutral, color=GalaxyInspector.WONG_ORANGE),
        fit=false,
        output_file=joinpath(output_path, "_ks_law/leroy2008_total.png"),
        trans_mode,
        filter_mode,
        sim_labels=[label],
        theme=Theme(
            Legend=(margin=(10, 0, 0, 20),),
            Axis=(
                limits=(0.4, 2.6, -4.5, 0.5),
                xticks=0.5:0.5:2.5,
                yticks=-4:1:0,
            ),
        ),
    )

    total_paths = [
        joinpath(output_path, "_ks_law/leroy2008_total.png"),
        joinpath(output_path, "_ks_law/bigiel2010_total.png"),
    ]

    GalaxyInspector.hvcatImages(
        1,
        total_paths;
        output_path=joinpath(output_path, "_ks_law/total.png"),
    )

    ###############################
    # Final image with all KS laws
    ###############################

    gas_paths = joinpath.(
        output_path,
        "_ks_law",
        ["molecular.png", "atomic.png", "neutral.png", "total.png"],
    )

    GalaxyInspector.hvcatImages(
        4,
        gas_paths;
        output_path=joinpath(output_path, "ks_law.png"),
    )

    rm(joinpath(output_path, "_ks_law"); recursive=true, force=true)

    ################################################################################################
    # Stellar density maps of the last snapshot, face-on/edge-on projections
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(
            log_file,
            "# Stellar density maps of the last snapshot, face-on/edge-on projections",
        )
        println(log_file, "#"^100, "\n")
    end

    stellarDensityMaps(
        [simulation_path],
        n_snapshots;
        box_size=R3,
        output_path,
        trans_mode,
        filter_mode,
    )

    ################################################################################################
    # Gas components density map of the last snapshot, face-on/edge-on projections
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(
            log_file,
            "# Gas components density map of the last snapshot, face-on/edge-on projections",
        )
        println(log_file, "#"^100, "\n")
    end

    gasDensityMaps(
        [simulation_path],
        n_snapshots;
        box_size=R3,
        output_path,
        density_range=(2.0, NaN),
        trans_mode,
        filter_mode,
    )

    ################################################################################################
    # Face-on density maps for different quantities of the last snapshot
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# Face-on density maps for different quantities of the last snapshot")
        println(log_file, "#"^100, "\n")
    end

    temp_folder = joinpath(output_path, "_face_on_maps")

    box_size       = R3
    pixel_length   = box_size / 400.0
    half_box_size  = ustrip(u"kpc", box_size) / 2.0
    tick           = floor(half_box_size; sigdigits=1)
    size           = (910, 760)
    figure_padding = (20, 20, 30, 20)

    densityMap(
        [simulation_path],
        n_snapshots;
        box_size,
        pixel_length,
        l_unit=u"kpc",
        output_path=temp_folder,
        trans_mode,
        filter_mode,
        title="Gas density",
        colorbar=true,
        theme=Theme(;
            size,
            figure_padding,
            Colorbar=(label=L"\log_{10} \, \Sigma_\text{gas} \,\, [\mathrm{M_\odot \, kpc^{-2}}]",),
            Axis=(
                limits=(-half_box_size, half_box_size, -half_box_size, half_box_size),
                xticks=-tick:10:tick,
                yticks=-tick:10:tick,
            ),
        ),
    )

    densityMap(
        [simulation_path],
        n_snapshots;
        components=[:stellar],
        field_types=[:particles],
        box_size,
        pixel_length,
        l_unit=u"kpc",
        output_path=temp_folder,
        trans_mode,
        filter_mode,
        title="Stellar density",
        colorbar=true,
        theme=Theme(;
            size,
            figure_padding,
            Colorbar=(label=L"\log_{10} \, \Sigma_\star \,\, [\mathrm{M_\odot \, kpc^{-2}}]",),
            Axis=(
                limits=(-half_box_size, half_box_size, -half_box_size, half_box_size),
                xticks=-tick:10:tick,
                yticks=-tick:10:tick,
            ),
        ),

    )

    temperatureMap(
        [simulation_path],
        n_snapshots;
        field_type=:cells,
        box_size,
        pixel_length,
        output_path=temp_folder,
        trans_mode,
        filter_mode,
        title="Gas temperature",
        colorbar=true,
        theme=Theme(;
            size,
            figure_padding,
            Colorbar=(label=L"\log_{10} \, T \,\, [\mathrm{K}]",),
            Axis=(
                limits=(-half_box_size, half_box_size, -half_box_size, half_box_size),
                xticks=-tick:10:tick,
                yticks=-tick:10:tick,
            ),
        ),
    )

    metallicityMap(
        [simulation_path],
        n_snapshots;
        box_size,
        pixel_length,
        output_path=temp_folder,
        trans_mode,
        filter_mode,
        title="Gas metallicity",
        colorbar=true,
        theme=Theme(;
            size,
            figure_padding,
            Colorbar=(label=L"\log_{10} \, Z_\text{gas} \, [\mathrm{Z_\odot}]",),
            Axis=(
                limits=(-half_box_size, half_box_size, -half_box_size, half_box_size),
                xticks=-tick:10:tick,
                yticks=-tick:10:tick,
            ),
        ),
    )

    metallicityMap(
        [simulation_path],
        n_snapshots;
        components=[:stellar],
        field_types=[:particles],
        output_path=temp_folder,
        box_size,
        pixel_length,
        trans_mode,
        filter_mode,
        title="Stellar metallicity",
        colorbar=true,
        colorrange=(0.04, 0.8),
        theme=Theme(;
            size,
            figure_padding,
            Colorbar=(label=L"\log_{10} \, Z_\star \, [\mathrm{Z_\odot}]",),
            Axis=(
                limits=(-half_box_size, half_box_size, -half_box_size, half_box_size),
                xticks=-tick:10:tick,
                yticks=-tick:10:tick,
            ),
        ),
    )

    gasSFRMap(
        [simulation_path],
        n_snapshots;
        box_size,
        pixel_length,
        output_path=temp_folder,
        trans_mode,
        filter_mode,
        title="Gas SFR",
        colorbar=true,
        theme=Theme(;
            size,
            figure_padding,
            Colorbar=(label=L"\log_{10} \, \text{SFR}_\text{gas} \, [\mathrm{M_\odot \, yr^{-1}}]",),
            Axis=(
                limits=(-half_box_size, half_box_size, -half_box_size, half_box_size),
                xticks=-tick:10:tick,
                yticks=-tick:10:tick,
            ),
        ),
    )

    GalaxyInspector.hvcatImages(
        3,
        joinpath.(
            temp_folder,
            [
                "$(basename(simulation_path))_gas_xy_density_map_snap_$(n_snaps_str).png",
                "$(basename(simulation_path))_stellar_xy_density_map_snap_$(n_snaps_str).png",
                "$(basename(simulation_path))_xy_gas_sfr_map_snap_$(n_snaps_str).png",
                "$(basename(simulation_path))_xy_gas_temperature_map_snap_$(n_snaps_str).png",
                "$(basename(simulation_path))_gas_xy_all_metallicity_map_snap_$(n_snaps_str).png",
                "$(basename(simulation_path))_stellar_xy_all_metallicity_map_snap_$(n_snaps_str).png",
            ]
        );
        output_path=joinpath(output_path, "face_on_maps.png"),
    )

    rm(temp_folder, recursive=true, force=true)

    ################################################################################################
    # Evolution of the masses and of the fractions, for the different gas components
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(
            log_file,
            "# Evolution of the masses and of the fractions, for the different gas components",
        )
        println(log_file, "#"^100, "\n")
    end

    gasFractionsEvolution(
        [simulation_path];
        r_gas=R1,
        output_path,
        trans_mode,
        filter_mode,
    )

    ################################################################################################
    # Stellar and gas density video, face-on/edge-on projections
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# Stellar and gas density video, face-on/edge-on projections")
        println(log_file, "#"^100, "\n")
    end

    evolutionVideo(
        [simulation_path],
        :ode_molecular;
        filed_type=:cells,
        box_size=R3,
        output_path,
        density_range=(4.0, NaN),
        trans_mode,
        filter_mode,
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

    # Check if the simulation is cosmological
    cosmological = all(GalaxyInspector.isSimCosmological, simulation_paths)
    if cosmological
        trans_mode   = :stellar_subhalo
        filter_mode  = :subhalo
        extra_filter = dd -> GalaxyInspector.filterBySphere(dd, (0.0u"kpc", R1), :zero)
    else
        trans_mode   = :stellar_box
        filter_mode  = :sphere
        extra_filter = GalaxyInspector.filterNothing
    end

    # ################################################################################################
    # # Evolution of the gas components and the SFR
    # ################################################################################################

    # if logging
    #     println(log_file, "#"^100)
    #     println(log_file, "# Evolution of the gas components and the SFR")
    #     println(log_file, "#"^100, "\n")
    # end

    # quantities   = [:ionized_mass, :atomic_mass, :molecular_mass, :dust_mass, :sfr, :stellar_mass]
    # smooths      = [0, 0, 0, 0, n_snapshots ÷ 5, 0]
    # n_quantities = length(quantities)

    # jld2_paths = Vector{String}(undef, n_quantities * n_sims)

    # for (i, simulation) in pairs(simulation_paths)

    #     sim_name = basename(simulation)

    #     for (j, (quantity, smooth)) in enumerate(zip(quantities, smooths))

    #         if logging
    #             println(log_file, "#"^100)
    #             println(
    #                 log_file,
    #                 "# Computing the evolution of $(quantity) for $(sim_name)",
    #             )
    #             println(log_file, "#"^100, "\n")
    #         end

    #         y_plot_params = GalaxyInspector.plotParams(quantity)

    #         plotTimeSeries(
    #             [simulation],
    #             [lines!];
    #             output_path=temp_folder,
    #             filename="$(quantity)_$(sim_name)",
    #             da_functions=[GalaxyInspector.daEvolution],
    #             da_args=[(:physical_time, quantity)],
    #             da_kwargs=[(; filter_mode, extra_filter, smooth, scaling=identity)],
    #             x_unit=u"Gyr",
    #             y_unit=y_plot_params.unit,
    #             y_exp_factor=y_plot_params.exp_factor,
    #             save_figures=false,
    #             backup_results=true,
    #         )

    #         jld2_paths[(i - 1) * n_quantities + j] = joinpath(
    #             temp_folder,
    #             "$(quantity)_$(sim_name).jld2",
    #         )

    #     end

    # end

    # current_theme = merge(
    #     Theme(palette=(linestyle=[:solid],),),
    #     GalaxyInspector.DEFAULT_THEME,
    #     theme_latexfonts(),
    # )

    # colors = current_theme[:palette][:color][]

    # with_theme(current_theme) do

    #     f = Figure(size=(2750, 1400), figure_padding=(5, 10, 5, 10))

    #     for (i, quantity) in pairs(quantities)

    #         y_plot_params = GalaxyInspector.plotParams(quantity)

    #         ylabel = LaTeXString(
    #             replace(
    #                 y_plot_params.axis_label,
    #                 "auto_label" => GalaxyInspector.getLabel(
    #                     y_plot_params.var_name,
    #                     y_plot_params.exp_factor,
    #                     y_plot_params.unit,
    #                 ),
    #             ),
    #         )

    #         ax = CairoMakie.Axis(
    #             f[ceil(Int, i / 3), mod1(i, 3)];
    #             xlabel=L"t \, [\mathrm{Gyr}]",
    #             ylabel,
    #             aspect=AxisAspect(1.4),
    #         )

    #         for (simulation_path, label, color) in zip(simulation_paths, labels, colors)

    #             jld2_path = joinpath(temp_folder, "$(quantity)_$(basename(simulation_path)).jld2")

    #             jldopen(jld2_path, "r") do jld2_file

    #                 address = first(keys(jld2_file))

    #                 x, y = jld2_file[address]["simulation_001"]

    #                 lines!(ax, x, y; label, color, linewidth=3)

    #             end

    #         end

    #         if i == length(quantities)

    #             axislegend(ax, position=:rb, framevisible=false, nbanks=1)

    #         end


    #     end

    #     save(joinpath(output_path, "mass_evolution.png"), f)

    # end

    # rm(temp_folder; recursive=true)

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
    #         simulation_paths,
    #         n_snapshots,
    #         quantity;
    #         output_path=joinpath(output_path, "Molla2015"),
    #         filter_mode,
    #         sim_labels=labels,
    #         theme=Theme(
    #             size=(1500, 880),
    #             figure_padding=(10, 15, 5, 15),
    #             palette=(linestyle=[:solid],),
    #             Axis=(aspect=nothing,),
    #             Legend=(halign=:right, valign=:top, nbanks=1),
    #         ),
    #     )

    # end

    kennicuttSchmidtLaw(
        simulation_paths,
        n_snapshots;
        quantity=:gas,
        reduce_grid=:circular,
        grid_size=30.0u"kpc",
        bin_size=1.5u"kpc",
        post_processing=GalaxyInspector.ppKennicutt1998!,
        pp_kwargs=(; color=GalaxyInspector.WONG_ORANGE, extend=2.0),
        fit=false,
        output_file=joinpath(output_path, "ks_law_total.png"),
        trans_mode,
        filter_mode,
        sim_labels=labels,
        theme=Theme(
            figure_padding=(2, 20, 5, 15),
            Legend=(margin=(10, 0, 0, 20),),
            Axis=(
                limits=(1.0, 2.5, -3.0, 0.5),
                xticks=1.0:0.5:2.5,
                yticks=-3.0:0.5:0.5,
            ),
        ),
    )

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

    # ##################################################
    # # R⊙ tests
    # ##################################################

    # comparison(
    #     joinpath.("F:/simulations/current/", ["SFM_064_T10", "SFM_064_T09", "SFM_064_T11"]),
    #     joinpath(BASE_OUT_PATH, "comparison/Rsun_test"),
    #     LOGGING;
    #     labels=[
    #         L"R_\odot = 1.0 \times 10^{-17} \, \mathrm{cm^3 \, s^{-1}}",
    #         L"R_\odot = 3.5 \times 10^{-17} \, \mathrm{cm^3 \, s^{-1}}",
    #         L"R_\odot = 5.0 \times 10^{-17} \, \mathrm{cm^3 \, s^{-1}}",
    #     ],
    # )

    # ##################################################
    # # Zeff tests
    # ##################################################

    # comparison(
    #     joinpath.("F:/simulations/current/", ["SFM_064_T12", "SFM_064_T09", "SFM_064_T13"]),
    #     joinpath(BASE_OUT_PATH, "comparison/Zeff_test"),
    #     LOGGING;
    #     labels=[
    #         L"Z_\mathrm{eff} = 10^{-4} \, Z_\odot",
    #         L"Z_\mathrm{eff} = 10^{-3} \, Z_\odot",
    #         L"Z_\mathrm{eff} = 10^{-2} \, Z_\odot",
    #     ],
    # )

    # ##################################################
    # # Crho tests
    # ##################################################

    # comparison(
    #     joinpath.("F:/simulations/current/", ["SFM_064_T14", "SFM_064_T15", "SFM_064_T09"]),
    #     joinpath(BASE_OUT_PATH, "comparison/Crho_test"),
    #     LOGGING;
    #     labels=[L"C_\rho = 1", L"C_\rho = 50", L"C_\rho = 100"],
    # )

    # ##################################################
    # # taudd tests
    # ##################################################

    # comparison(
    #     joinpath.(
    #         "F:/simulations/current/",
    #         ["SFM_064_T16", "SFM_064_T09", "SFM_064_T17", "SFM_064_T18"],
    #     ),
    #     joinpath(BASE_OUT_PATH, "comparison/taudd_test"),
    #     LOGGING;
    #     labels=[
    #         L"\tau_\mathrm{dd} = 2.0 \, \mathrm{Gyr}",
    #         L"\tau_\mathrm{dd} = 2.3 \, \mathrm{Gyr}",
    #         L"\tau_\mathrm{dd} = 3.2 \, \mathrm{Gyr}",
    #         L"\tau_\mathrm{dd} = 20.0 \, \mathrm{Myr}",
    #     ],
    # )

    # ##################################################
    # # Blitz and STD tests
    # ##################################################

    # comparison(
    #     joinpath.(
    #         "F:/simulations/",
    #         ["2025-10-08/isolated/BLT_064", "2025-10-08/isolated/STD_064", "current/SFM_064_T09"],
    #     ),
    #     joinpath(BASE_OUT_PATH, "comparison/Blitz_test"),
    #     LOGGING;
    #     labels=["BLT", "STD", "SFM"],
    # )

    ##################################################
    # Papers tests
    ##################################################

    # comparison(
    #     joinpath.(
    #         "F:/simulations/",
    #         ["2025-10-08/isolated/SFM_064_T03", "current/SFM_064_T09"],
    #     ),
    #     joinpath(BASE_OUT_PATH, "comparison/paper_test"),
    #     LOGGING;
    #     labels=["Paper 1", "Paper 2"],
    # )

    ##################################################

    SIMULATIONS = joinpath.(
        "F:/simulations/",
        [
            "current/SFM_064_T09",
            "current/SFM_128_T19",
            "2025-10-08/cosmological/test_dust_13",
        ]
    )

    for simulation in SIMULATIONS
        basic_analysis(simulation, BASE_OUT_PATH, LOGGING)
    end

end
