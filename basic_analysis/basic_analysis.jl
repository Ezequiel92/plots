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
    InvertedIndices,
    JLD2,
    LaTeXStrings,
    LinearAlgebra,
    Logging,
    Measurements,
    NearestNeighbors,
    ProgressMeter,
    QuadGK,
    Rotations,
    StaticArrays,
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
    video::Bool=false,
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

    # Compute the number of snapshots
    n_snapshots = GalaxyInspector.countSnapshot(simulation_path)

    # Compute the number of snapshots as a string
    n_snaps_str = lpad(string(n_snapshots - 1), 3, "0")

    (
        !iszero(n_snapshots) ||
        throw(ArgumentError("basic_analysis: $(simulation_path) has no snapshots"))
    )

    # Check if the simulation is cosmological
    cosmological = GalaxyInspector.isSimCosmological(simulation_path)

    # Check if the simulation has our SF model
    sfm = GalaxyInspector.isSimSFM(simulation_path)

    # Characteristic radii
    R1 = 40.0u"kpc"
    R2 = 2.0u"kpc"
    R3 = 65.0u"kpc"

    # Choose the correct transformations and filters
    if cosmological
        trans_mode   = :stellar_subhalo
        filter_mode  = :subhalo
        extra_filter = dd -> GalaxyInspector.filterBySphere(dd, (0.0u"kpc", R1), :zero)
    else
        trans_mode   = :stellar_box
        filter_mode  = :sphere
        extra_filter = GalaxyInspector.filterNothing
    end

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
        range=(1.0e-6, 1.0),
        output_path,
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

        #TODO
        molla_quantities = [
            :stellar_area_density,
            :ode_molecular_stellar_area_density,
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
            output_path=joinpath(output_path, "molla_2015"),
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

    temp_folder = joinpath(output_path, "_ks_law")

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
        output_file=joinpath(temp_folder, "sun2023_molecular.png"),
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
        output_file=joinpath(temp_folder, "leroy2008_molecular.png"),
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
        joinpath(temp_folder, "leroy2008_molecular.png"),
        joinpath(temp_folder, "sun2023_molecular.png"),
    ]

    GalaxyInspector.hvcatImages(
        1,
        molecular_paths;
        output_path=joinpath(temp_folder, "molecular.png"),
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
        output_file=joinpath(temp_folder, "bigiel2010_atomic.png"),
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
        output_file=joinpath(temp_folder, "leroy2008_atomic.png"),
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
        joinpath(temp_folder, "leroy2008_atomic.png"),
        joinpath(temp_folder, "bigiel2010_atomic.png"),
    ]

    GalaxyInspector.hvcatImages(
        1,
        atomic_paths;
        output_path=joinpath(temp_folder, "atomic.png"),
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
        output_file=joinpath(temp_folder, "bigiel2010_neutral.png"),
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
        output_file=joinpath(temp_folder, "leroy2008_neutral.png"),
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
        joinpath(temp_folder, "leroy2008_neutral.png"),
        joinpath(temp_folder, "bigiel2010_neutral.png"),
    ]

    GalaxyInspector.hvcatImages(
        1,
        neutral_paths;
        output_path=joinpath(temp_folder, "neutral.png"),
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
        output_file=joinpath(temp_folder, "bigiel2010_total.png"),
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
        output_file=joinpath(temp_folder, "leroy2008_total.png"),
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
        joinpath(temp_folder, "leroy2008_total.png"),
        joinpath(temp_folder, "bigiel2010_total.png"),
    ]

    GalaxyInspector.hvcatImages(
        1,
        total_paths;
        output_path=joinpath(temp_folder, "total.png"),
    )

    ###############################
    # Final image with all KS laws
    ###############################

    GalaxyInspector.hvcatImages(
        4,
        joinpath.(temp_folder, ["molecular.png", "atomic.png", "neutral.png", "total.png"]);
        output_path=joinpath(output_path, "ks_law.png"),
    )

    rm(temp_folder; recursive=true, force=true)

    ################################################################################################
    # Stellar density maps of the last snapshot, face-on/edge-on projections
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# Stellar density maps of the last snapshot, face-on/edge-on projections")
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
    # Gas density maps of the last snapshot, face-on/edge-on projections
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# Gas density maps of the last snapshot, face-on/edge-on projections")
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
        colorrange=(0.0, 1.0),
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
    # Evolution of the masses and fractions, for the different gas components
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# Evolution of the masses and fractions, for the different gas components")
        println(log_file, "#"^100, "\n")
    end

    gasFractionsEvolution(
        [simulation_path];
        r_gas=R1,
        output_path,
        trans_mode,
        filter_mode,
    )

    if video

        ############################################################################################
        # Stellar and gas density video, face-on/edge-on projections
        ############################################################################################

        if logging
            println(log_file, "#"^100)
            println(log_file, "# Stellar and gas density video, face-on/edge-on projections")
            println(log_file, "#"^100, "\n")
        end

        # Choose the correct transformations and filters
        if cosmological
            trans_mode   = :stellar_subhalo
            filter_mode  = :subhalo
            extra_filter = dd -> GalaxyInspector.filterBySphere(dd, (0.0u"kpc", R1), :zero)
        else
            trans_mode   = :stellar_box
            filter_mode  = :all
            extra_filter = GalaxyInspector.filterNothing
        end

        #TODO
        evolutionVideo(
            [simulation_path],
            :ode_molecular_stellar;
            field_type=:cells,
            box_size=R3,
            output_path,
            density_range=(NaN, NaN),
            trans_mode,
            filter_mode,
            save_data=true,
        )

    end

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

    # Compute the number of simulations
    n_sims = length(simulation_paths)
    @assert n_sims == length(labels) "Number of simulations and labels must match."

    # Create the output folder
    mkpath(output_path)

    # Activate logging
    if logging
        log_file = open(joinpath(output_path, "logs.txt"), "w+")
        GalaxyInspector.setLogging!(logging; stream=log_file)
    end

    # Compute the number of snapshots
    n_snapshots = minimum(GalaxyInspector.countSnapshot.(simulation_paths))

    # Check if the simulation are cosmological
    cosmological = all(GalaxyInspector.isSimCosmological, simulation_paths)

    # Check if the simulations have out SF model
    sfm = all(GalaxyInspector.isSimSFM, simulation_paths)

    # Characteristic radii
    R1 = 40.0u"kpc"

    # Choose the correct transformations and filters
    if cosmological
        trans_mode   = :gas_subhalo
        filter_mode  = :subhalo
        extra_filter = dd -> GalaxyInspector.filterBySphere(dd, 0.0u"kpc", R1, :zero)
    else
        trans_mode   = :all_box
        filter_mode  = :sphere
        extra_filter = GalaxyInspector.filterNothing
    end

    ################################################################################################
    # Evolution of the gas and stellar components
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# Evolution of the gas and stellar components")
        println(log_file, "#"^100, "\n")
    end

    temp_folder = joinpath(output_path, "_mass_evolution")

    if sfm

        quantities = [
            :gas_mass,
            :ode_ionized_mass,
            :ode_atomic_mass,
            :ode_molecular_stellar_mass,
            :ode_metals_mass,
            :ode_dust_mass,
            :ode_neutral_mass,
            :ode_cold_mass,
            :stellar_mass,
            :sfr,
        ]

    else

        quantities = [
            :gas_mass,
            :ionized_mass,
            :br_atomic_mass,
            :br_molecular_mass,
            :neutral_mass,
            :Z_gas_mass,
            :stellar_mass,
            :sfr,
        ]

    end

    ff_request = Dict(:stellar => ["POS "])

    for quantity in quantities

        if logging
            println(log_file, "#"^100)
            println(log_file, "# Evolution of $(quantity)")
            println(log_file, "#"^100, "\n")
        end

        y_plot_params = GalaxyInspector.plotParams(quantity)

        plotTimeSeries(
            simulation_paths,
            [lines!];
            output_path=temp_folder,
            filename="$(quantity)",
            da_functions=[GalaxyInspector.daEvolution],
            da_args=[(:physical_time, quantity)],
            da_kwargs=[(; trans_mode, filter_mode, extra_filter, ff_request)],
            x_unit=u"Gyr",
            y_unit=y_plot_params.unit,
            y_exp_factor=y_plot_params.exp_factor,
            xaxis_label=L"t \, [\mathrm{Gyr}]",
            yaxis_var_name=y_plot_params.var_name,
            sim_labels=labels,
            theme=Theme(
                size=(1300, 880),
                figure_padding=(10, 15, 5, 15),
                fontsize=45,
                palette=(linestyle=[:solid],),
                Axis=(aspect=nothing,),
                Legend=(halign=:left, valign=:top, nbanks=1, margin=(10, 0, 0, 0),)
            )
        )

    end

    img_paths = joinpath.(temp_folder, ["$(quantity).png" for quantity in quantities])

    GalaxyInspector.hvcatImages(
        Int(length(img_paths) / 2),
        img_paths;
        output_path=joinpath(output_path, "mass_evolution.png"),
    )

    rm(temp_folder; recursive=true, force=true)

    ################################################################################################
    # Profiles for the last snapshot, comparison with Mollá et al. (2015)
    ################################################################################################

    # Choose the correct transformations and filters
    if cosmological
        trans_mode   = :stellar_subhalo
        filter_mode  = :subhalo
        extra_filter = dd -> GalaxyInspector.filterBySphere(dd, (0.0u"kpc", R1), :zero)
    else
        trans_mode   = :stellar_box
        filter_mode  = :sphere
        extra_filter = GalaxyInspector.filterNothing
    end

    #TODO
    if sfm
        quantities = [
        :stellar_area_density,
        :ode_molecular_stellar_area_density,
        :sfr_area_density,
        :ode_atomic_area_density,
        :O_stellar_abundance,
        :N_stellar_abundance,
        :C_stellar_abundance,
    ]
    else
        quantities = [
        :stellar_area_density,
        :br_molecular_area_density,
        :sfr_area_density,
        :br_atomic_area_density,
        :O_stellar_abundance,
        :N_stellar_abundance,
        :C_stellar_abundance,
    ]
    end

    for quantity in quantities

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
            output_path=joinpath(output_path, "molla_2015"),
            trans_mode,
            filter_mode,
            sim_labels=labels,
            theme=Theme(Legend=(halign=:right, valign=:top,)),
        )

    end

    ################################################################################################
    # Resolved Kennicutt–Schmidt law of the last snapshot, scatter with circular grid
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(
            log_file,
            "# Resolved Kennicutt–Schmidt law of the last snapshot, scatter with circular grid",
        )
        println(log_file, "#"^100, "\n")
    end

    temp_folder = joinpath(output_path, "_ks_law")

    ###################
    # Molecular KS law
    ###################

    if logging
        println(log_file, "#"^50)
        println(log_file, "# Molecular KS law")
        println(log_file, "#"^50, "\n")
    end

    kennicuttSchmidtLaw(
        simulation_paths,
        n_snapshots;
        quantity=:molecular,
        reduce_grid=:circular,
        grid_size=30.0u"kpc",
        bin_size=1.5u"kpc",
        post_processing=GalaxyInspector.ppSun2023!,
        pp_kwargs=(; color=GalaxyInspector.WONG_ORANGE),
        fit=false,
        output_file=joinpath(temp_folder, "sun2023_molecular.png"),
        trans_mode,
        filter_mode,
        sim_labels=labels,
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
        simulation_paths,
        n_snapshots;
        quantity=:molecular,
        reduce_grid=:circular,
        grid_size=30.0u"kpc",
        bin_size=1.5u"kpc",
        post_processing=GalaxyInspector.ppLeroy2008!,
        pp_kwargs=(; color=GalaxyInspector.WONG_ORANGE),
        fit=false,
        output_file=joinpath(temp_folder, "leroy2008_molecular.png"),
        trans_mode,
        filter_mode,
        sim_labels=labels,
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
        joinpath(temp_folder, "leroy2008_molecular.png"),
        joinpath(temp_folder, "sun2023_molecular.png"),
    ]

    GalaxyInspector.hvcatImages(
        1,
        molecular_paths;
        output_path=joinpath(temp_folder, "molecular.png"),
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
        simulation_paths,
        n_snapshots;
        quantity=:gas,
        reduce_grid=:circular,
        grid_size=30.0u"kpc",
        bin_size=1.5u"kpc",
        post_processing=GalaxyInspector.ppBigiel2010!,
        pp_kwargs=(; galaxy=:all, quantity=:neutral, color=GalaxyInspector.WONG_ORANGE),
        fit=false,
        output_file=joinpath(temp_folder, "bigiel2010_total.png"),
        trans_mode,
        filter_mode,
        sim_labels=labels,
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
        simulation_paths,
        n_snapshots;
        quantity=:gas,
        reduce_grid=:circular,
        grid_size=30.0u"kpc",
        bin_size=1.5u"kpc",
        post_processing=GalaxyInspector.ppLeroy2008!,
        pp_kwargs=(; quantity=:neutral, color=GalaxyInspector.WONG_ORANGE),
        fit=false,
        output_file=joinpath(temp_folder, "leroy2008_total.png"),
        trans_mode,
        filter_mode,
        sim_labels=labels,
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
        joinpath(temp_folder, "leroy2008_total.png"),
        joinpath(temp_folder, "bigiel2010_total.png"),
    ]

    GalaxyInspector.hvcatImages(
        1,
        total_paths;
        output_path=joinpath(temp_folder, "total.png"),
    )

    ###############################
    # Final image with all KS laws
    ###############################

    GalaxyInspector.hvcatImages(
        2,
        joinpath.(temp_folder, ["molecular.png", "total.png"]);
        output_path=joinpath(output_path, "ks_law.png"),
    )

    rm(temp_folder; recursive=true, force=true)

    ################################################################################################
    # Close files
    ################################################################################################

    logging && close(log_file)

    return nothing

end

function (@main)(ARGS)

    # If logging into a file will be enable
    LOGGING = true

    # If a video of the evolution will be made (vary slow)
    VIDEO = false

    # Output folder
    BASE_OUT_PATH = "./"

    ################################################################################################
    # Comparison
    ################################################################################################

    # #############
    # # Crho tests
    # #############

    # comparison(
    #     joinpath.("F:/simulations/current/", ["SFM_064_T14", "SFM_064_T15", "SFM_064_T09"]),
    #     joinpath(BASE_OUT_PATH, "comparison/Crho_test"),
    #     LOGGING;
    #     labels=[L"C_\rho = 1", L"C_\rho = 50", L"C_\rho = 100"],
    # )

    # #############
    # # R⊙ tests
    # #############

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

    # #############
    # # Zeff tests
    # #############

    # comparison(
    #     joinpath.("F:/simulations/current/", ["SFM_064_T12", "SFM_064_T09", "SFM_064_T13"]),
    #     joinpath(BASE_OUT_PATH, "comparison/Zeff_test"),
    #     LOGGING;
    #     labels=[
    #         L"Z_\text{eff} = 1.0 \times 10^{-4} \, Z_\odot",
    #         L"Z_\text{eff} = 1.0 \times 10^{-3} \, Z_\odot",
    #         L"Z_\text{eff} = 1.0 \times 10^{-2} \, Z_\odot",
    #     ],
    # )

    ################################################################################################
    # Basic analysis
    ################################################################################################

    SIMULATIONS = joinpath.(
        "F:/simulations/",
        [
            "current/SFM_01",
            "current/SFM_06",
            "current/SFM_06_CRHO50",
        ]
    )

    for simulation in SIMULATIONS
        basic_analysis(simulation, BASE_OUT_PATH, LOGGING; video=VIDEO)
    end

    # ###################
    # # Paper comparison
    # ###################

    # comparison(
    #     joinpath.("F:/simulations/current/", ["SFM_01", "SFM_06"]),
    #     joinpath(BASE_OUT_PATH, "comparison/paper_comparison"),
    #     LOGGING;
    #     labels=[
    #         "Paper I",
    #         "Paper II",
    #     ],
    # )

    # #####################
    # # Cluster comparison
    # #####################

    # comparison(
    #     joinpath.("F:/simulations/current/", ["SFM_06", "SFM_06_CRHO50"]),
    #     joinpath(BASE_OUT_PATH, "comparison/cluster_comparison"),
    #     LOGGING;
    #     labels=[
    #         L"C_\rho = 100",
    #         L"C_\rho = 50",
    #     ],
    # )

end
