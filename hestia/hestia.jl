cd(@__DIR__)

using Pkg

Pkg.activate("../../codes/GalaxyInspector/")
Pkg.instantiate()

using CairoMakie, JLD2, LaTeXStrings, Unitful, UnitfulAstro, Statistics, ProgressMeter

push!(LOAD_PATH, "../../codes/GalaxyInspector/src/")
using GalaxyInspector

function hestia(
    simulation_path::String,
    subhalo_idx::Dict{Symbol, Int},
    output_path::String,
    logging::Bool;
    simulation_label::AbstractString=basename(simulation_path),
)::Nothing

    output_folder = joinpath(output_path, simulation_label)

    # Create the output folder
    mkpath(output_folder)

    # Activate logging
    if logging
        log_file = open(joinpath(output_path, "logs.txt"), "w+")
        GalaxyInspector.setLogging!(logging; stream=log_file)
    end

    # Filter and transformation request
    trans_request = GalaxyInspector.mergeRequests(
        Dict(:stellar => ["POS ", "MASS", "VEL "]),
        Dict(
            :group    => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
            :subhalo  => ["S_LenType", "S_Pos", "S_Vel"],
        ),
    )

    # Choose the correct transformations and filters
    trans_mode_MW = (
        # Center at subhalo subhalo_idx[:MW]
        subhalo_idx[:MW] + 1,
        # Rotate by the PA of the stars in subhalo subhalo_idx[:MW]
        (:pa, :stellar, dd->GalaxyInspector.filterBySubhalo(dd, subhalo_idx[:MW] + 1)),
        # Transformation request
        trans_request,
    )
    trans_mode_M31 = (
        # Center at subhalo subhalo_idx[:MW]
        subhalo_idx[:M31] + 1,
        # Rotate by the PA of the stars in subhalo subhalo_idx[:M31]
        (:pa, :stellar, dd->GalaxyInspector.filterBySubhalo(dd, subhalo_idx[:M31] + 1)),
        # Transformation request
        trans_request,
    )

    filter_mode_MW = (
        # Filter function. Subhalo with absolute index subhalo_idx[:MW]
        dd -> GalaxyInspector.filterBySubhalo(dd, subhalo_idx[:MW] + 1),
        # Filter request, needed data already in the transformation request
        Dict{Symbol,Vector{String}}(),
    )
    filter_mode_M31 = (
        # Filter function. Subhalo with absolute index subhalo_idx[:M31]
        dd -> GalaxyInspector.filterBySubhalo(dd, subhalo_idx[:M31] + 1),
        # Filter request, needed data already in the transformation request
        Dict{Symbol,Vector{String}}(),
    )

    ################################################################################################
    # Report files
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# Report files")
        println(log_file, "#"^100, "\n")
    end

    simulationReport([simulation_path]; output_path=output_folder)

    for (trans_mode, filter_mode, galaxy_label) in zip(
        [trans_mode_MW, trans_mode_M31],
        [filter_mode_MW, filter_mode_M31],
        ["MW", "M31"]
    )

        ############################################################################################
        # SDSS mockup of the last snapshot, face-on/edge-on projections
        ############################################################################################

        if logging
            println(log_file, "#"^100)
            println(log_file, "# SDSS mockup of the last snapshot, face-on/edge-on projections")
            println(log_file, "#"^100, "\n")
        end

        base_request = GalaxyInspector.mergeRequests(
            GalaxyInspector.plotParams(:stellar_age).request,
            GalaxyInspector.plotParams(:stellar_metallicity).request,
            GalaxyInspector.plotParams(:stellar_mass).request,
            GalaxyInspector.plotParams(:neutral_mass).request,
            Dict(:stellar=>["POS "], :gas=>["POS ", "RHO "]),
        )

        translation, rotation, trans_request = GalaxyInspector.selectTransformation(
            trans_mode,
            base_request,
        )
        filter_function, request = GalaxyInspector.selectFilter(filter_mode, trans_request)

        # Set up the grid
        grid = GalaxyInspector.CubicGrid(80.0u"kpc", 800)

        temp_folder = joinpath(output_folder, "_sdss_mockups")

        projection_planes = [:xy, :xz]

        plotSnapshot(
            [simulation_path, simulation_path],
            request,
            [heatmap!];
            output_path=temp_folder,
            base_filename="sdss_mockup",
            show_progress=true,
            slice=1,
            transform_box=true,
            translation,
            rotation,
            filter_function,
            da_functions=[GalaxyInspector.daSDSSMockup],
            da_args=[(grid,)],
            da_kwargs=[
                (;
                    projection_plane=projection_planes[1],
                    extinction=:neutral,
                    mmap_path=temp_folder,
                ),
                (;
                    projection_plane=projection_planes[2],
                    extinction=:neutral,
                    mmap_path=temp_folder,
                ),
            ],
            save_figures=false,
            backup_results=false,
            backup_raw_results=true,
            sim_labels=string.(projection_planes),
        )

        jld2_path = joinpath(temp_folder, "sdss_mockup_raw.jld2")

        l_unit = u"kpc"
        limit  = ustrip(l_unit, 80.0u"kpc" / 2.0)
        x_axis = extrema(ustrip.(l_unit, grid.x_axis))
        y_axis = extrema(ustrip.(l_unit, grid.y_axis))

        iterator = zip(projection_planes, [limit, limit / 2.0], [40, 20], [860, 500])

        for (projection_plane, ylimit, ytick, ysize) in iterator

            pp_string = string(projection_plane)
            xlabel    = GalaxyInspector.getLabel(pp_string[1:1], 0, l_unit)
            ylabel    = GalaxyInspector.getLabel(pp_string[2:2], 0, l_unit)

            current_theme = merge(
                Theme(
                    size=(880, ysize),
                    figure_padding=(5, 20, 0, 0),
                    Axis=(
                        aspect=DataAspect(),
                        limits=(-limit, limit, -ylimit, ylimit),
                        xticks=-40:10:40,
                        yticks=-ytick:10:ytick,
                    ),
                    Text=(color=:white, fontsize=30),
                ),
                GalaxyInspector.DEFAULT_THEME,
                theme_latexfonts(),
            )

            with_theme(current_theme) do

                jldopen(jld2_path, "r") do jld2_file

                    for snap in keys(jld2_file)

                        dataset_key = "$(snap)/$(projection_plane)"

                        f = Figure()

                        ax = CairoMakie.Axis(f[1, 1]; xlabel, ylabel)

                        image!(ax, x_axis, y_axis, jld2_file[dataset_key])

                        text!(
                            ax,
                            0.5,
                            0.98;
                            text=simulation_label * " - " * galaxy_label,
                            align=(:center, :top),
                            space=:relative,
                        )

                        output_file = "$(basename(simulation_path))_$(galaxy_label)_$(snap)_sdss_mockup_$(pp_string).png"

                        Makie.save(joinpath(output_folder, output_file), f)

                    end

                end

            end

        end

        rm(temp_folder; recursive=true)

    end

    ################################################################################################
    # Close files
    ################################################################################################

    logging && close(log_file)

    return nothing

end

function @main(args)

    # If logging into a file will be enable
    LOGGING = true

    # Output folder
    BASE_OUT_PATH = "./"

    # Set the correct internal length unit for Hestia
    GalaxyInspector.DEFAULT_L_UNIT[] = 1.0u"Mpc"

    ################################################################################################
    # Basic analysis
    ################################################################################################

    SIMULATIONS = joinpath.(
        "/run/media/elozano/DISK4T/hestia/",
        [
            # "17_11",
            "09_18",
            # "37_11",
        ]
    )

    # Subhalo number for the Milky-way and M31 at redshift 0
    SUBHALO_IDX = [
        # Dict(:MW => 1, :M31 => 0),
        Dict(:MW => 3911, :M31 => 2608),
        # Dict(:MW => 920, :M31 => 0),
    ]

    SIMULATION_LABELS = basename.(SIMULATIONS)

    iterator = zip(SIMULATIONS, SUBHALO_IDX, SIMULATION_LABELS)

    for (simulation_path, subhalo_idx, simulation_label) in iterator

        hestia(simulation_path, subhalo_idx, BASE_OUT_PATH, LOGGING; simulation_label)

    end

end
