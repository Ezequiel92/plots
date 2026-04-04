cd(@__DIR__)

using Pkg

Pkg.activate("../../codes/GalaxyInspector/")
Pkg.instantiate()

using CairoMakie, JLD2, LaTeXStrings, Unitful, UnitfulAstro, Statistics, ProgressMeter

push!(LOAD_PATH, "../../codes/GalaxyInspector/src/")
using GalaxyInspector

function lozano_2026(
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

    # Characteristic radii
    R1 = 40.0u"kpc"

    # Choose the correct transformations and filters
    trans_mode   = :gas_subhalo
    filter_mode  = :subhalo
    extra_filter = dd -> GalaxyInspector.filterBySphere(dd, 0.0u"kpc", R1, :zero)

    ################################################################################################
    # Efficiency per free-fall time histograms of the last snapshot
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# Efficiency per free-fall time histograms of the last snapshot")
        println(log_file, "#"^100, "\n")
    end

    plot_params = GalaxyInspector.plotParams(:gas_eff)

    translation, rotation, trans_request = GalaxyInspector.selectTransformation(
        trans_mode,
        plot_params.request,
    )
    filter_function, request = GalaxyInspector.selectFilter(filter_mode, trans_request)

    grid = GalaxyInspector.LinearGrid(1.0e-6, 1.0, 100; log=true)

    xaxis_label = L"\log_{10} \," * GalaxyInspector.getLabel(
        plot_params.axis_label,
        0,
        plot_params.unit,
    )

    GalaxyInspector.plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        output_path,
        base_filename="eff_histogram",
        slice=128,
        transform_box=true,
        translation,
        rotation,
        filter_function,
        da_functions=[GalaxyInspector.daHistogram],
        da_args=[(:gas_eff, grid)],
        da_kwargs=[(; log=true, filter_function=extra_filter)],
        post_processing=GalaxyInspector.ppVerticalFlags!,
        pp_args=([-2.0],),
        pp_kwargs=(; colors=[:black]),
        xaxis_label,
        xaxis_var_name=plot_params.var_name,
        yaxis_var_name=L"\mathrm{Normalized \,\, counts}",
        theme=Theme(
            palette=(linestyle=[:solid],),
            Legend=(nbanks=1, halign=:left, valign=:top, margin=(10, 0, 0, 5), labelsize=25),
            Lines=(linewidth=2,),
        ),
        sim_labels=labels,
    )

    ################################################################################################
    # Efficiency per free-fall time histograms of the last snapshot
    ################################################################################################

    if logging
        println(log_file, "#"^100)
        println(log_file, "# Efficiency per free-fall time histograms of the last snapshot")
        println(log_file, "#"^100, "\n")
    end

    molecularFractionEvolution(
        simulation_paths;
        slice=50:128,
        x_quantity=:redshift,
        measurements=false,
        output_path,
        trans_mode,
        filter_mode,
        sim_labels=labels,
        theme=Theme(Axis=(xreversed=true,), Legend=(labelsize=25,)),
    )

    ################################################################################
    # Gas efficiency per free-fall vs time, gas number density, and gas metallicity
    ################################################################################

    if logging
        println(log_file, "#"^100)
        println(
            log_file,
            "# Gas efficiency per free-fall vs time, gas number density, and gas metallicity",
        )
        println(log_file, "#"^100, "\n")
    end

    with_theme(merge(
        Theme(size=(1700, 670), palette=(linestyle=[:solid],), Lines=(linewidth=3,)),
        theme_latexfonts(),
        GalaxyInspector.DEFAULT_THEME,
    )) do

        f = Figure()

        t_limit = 0.2u"Gyr"

        ϵff_plot_params = GalaxyInspector.plotParams(:gas_eff)
        n_plot_params   = GalaxyInspector.plotParams(:gas_number_density)
        Z_plot_params   = GalaxyInspector.plotParams(:gas_metallicity)

        translation, rotation, trans_request = GalaxyInspector.selectTransformation(
            trans_mode,
            GalaxyInspector.mergeRequests(
                n_plot_params.request,
                Z_plot_params.request,
                ϵff_plot_params.request,
            ),
        )
        filter_function, request = GalaxyInspector.selectFilter(filter_mode, trans_request)

        ax_01 = CairoMakie.Axis(
            f[1, 1],
            xlabel=L"t \,\, [\mathrm{Gyr}]",
            ylabel=LaTeXString(
                replace(
                    L"$\log_{10} \, $auto_label",
                    "auto_label" => GalaxyInspector.getLabel(
                        ϵff_plot_params.var_name,
                        ϵff_plot_params.exp_factor,
                        ϵff_plot_params.unit,
                    ),
                ),
            ),
            xticks=0:2:14,
            yticks=-3.0:0.5:-0.5,
            limits=(nothing, nothing, -3.2, -0.4),
        )

        ax_02 = CairoMakie.Axis(
            f[1, 2];
            xlabel=LaTeXString(
                replace(
                    L"$\log_{10} \, $auto_label",
                    "auto_label" => GalaxyInspector.getLabel(
                        n_plot_params.var_name,
                        n_plot_params.exp_factor,
                        n_plot_params.unit,
                    ),
                ),
            ),
            ylabelvisible=false,
            yticklabelsvisible=false,
            limits=(nothing, nothing, -3.2, -0.4),
        )

        ax_03 = CairoMakie.Axis(
            f[1, 3];
            xlabel=LaTeXString(
                replace(
                    L"$\log_{10} \, $auto_label",
                    "auto_label" => GalaxyInspector.getLabel(
                        Z_plot_params.var_name,
                        Z_plot_params.exp_factor,
                        Z_plot_params.unit,
                    ),
                ),
            ),
            ylabelvisible=false,
            yticklabelsvisible=false,
            limits=(nothing, nothing, -3.2, -0.4),
        )

        # Start at `t_limit`
        initial_snap = GalaxyInspector.findClosestSnapshot(simulation_paths[1], t_limit)

        # Snapshot range
        slice = initial_snap:128

        for (simulation, label) in zip(simulation_paths, labels)

            #########################################################################
            # Gas efficiency per free-fall vs time (from t_limit up to t = 13.8 Gyr)
            #########################################################################

            # Median of ϵff (50th percentile)
            ys      = Vector{Float64}(undef, length(slice))
            # ϵff 25th percentile
            ys_low  = Vector{Float64}(undef, length(slice))
            # ϵff 75th percentile
            ys_high = Vector{Float64}(undef, length(slice))
            # Physical time
            times   = Vector{Float64}(undef, length(slice))

            prog_bar = Progress(
                length(slice),
                dt=0.5,
                desc="Analyzing and plotting the data... ",
                color=:blue,
                barglyphs=BarGlyphs("|#  |"),
            )

            for (i, snap_n) in pairs(slice)

                data_dict = makeDataDict(simulation, snap_n, request)

                # Filter and transform the values
                GalaxyInspector.translateData!(data_dict, translation...)
                GalaxyInspector.rotateData!(data_dict, rotation...)
                GalaxyInspector.filterData!(data_dict; filter_function)
                GalaxyInspector.filterData!(data_dict; filter_function=extra_filter)

                ϵffs = GalaxyInspector.scatterQty(data_dict, :gas_eff)

                # Ignore the ϵffs that are zero and NaN
                filter!(x->!iszero(x) && !isnan(x), ϵffs)

                # Median of ϵff (50th percentile)
                ys[i]      = median(ϵffs)
                # ϵff 25th percentile
                ys_low[i]  = quantile(ϵffs, 0.25)
                # ϵff 75th percentile
                ys_high[i] = quantile(ϵffs, 0.75)
                # Physical time
                times[i]   = ustrip(u"Gyr", data_dict[:snap_data].physical_time)

                next!(prog_bar)

            end

            # Go a to a log scale for the y axis
            ys      = log10.(ys)
            ys_low  = log10.(ys_low)
            ys_high = log10.(ys_high)

            # Smooth the values
            _, ϵffs = GalaxyInspector.smoothWindow(times, ys, 40)
            _, ϵffs_low = GalaxyInspector.smoothWindow(times, ys_low, 40)
            times, ϵffs_high = GalaxyInspector.smoothWindow(times, ys_high, 40)

            lines!(ax_01, times, ϵffs; label)
            band!(ax_01, times, ϵffs_low, ϵffs_high)

            #################################################
            # Gas efficiency per free-fall vs number density
            #################################################

            temp_folder = joinpath(output_path, "_eff_vs_number_density")

            plotSnapshot(
                [simulation],
                request,
                [scatter!];
                pf_kwargs=[(;)],
                output_path=temp_folder,
                base_filename="eff_vs_number_density",
                slice=128,
                filter_function,
                da_functions=[GalaxyInspector.daScatterGalaxy],
                da_args=[(:gas_number_density, :gas_eff)],
                da_kwargs=[(; filter_function=extra_filter)],
                transform_box=true,
                translation,
                rotation,
                x_unit=n_plot_params.unit,
                y_unit=ϵff_plot_params.unit,
                save_figures=false,
                backup_results=true,
            )

            jld2_path = joinpath(temp_folder, "eff_vs_number_density.jld2")

            jldopen(jld2_path, "r") do jld2_file

                n_gas, ϵ_ff = jld2_file["snap_127"][basename(simulation)]

                # Remove point where n or ϵff are zero or NaN
                zero_filter = (x, y) -> iszero(x) || iszero(y) || isnan(x) || isnan(y)
                delete_idxs = zero_filter.(n_gas, ϵ_ff)
                deleteat!(n_gas, delete_idxs)
                deleteat!(ϵ_ff, delete_idxs)

                # Set up a logarithmic grid for the x axis
                grid = GalaxyInspector.LinearGrid(extrema(n_gas)..., 40; log=true)

                # Set up a logarithmic scale for the y axis
                log_ϵ_ff = log10.(ϵ_ff)

                binned_ϵffs = GalaxyInspector.listHistogram1D(n_gas, log_ϵ_ff, grid)

                # Gas number density
                x_ngas = grid.x_axis
                # Median of ϵff (50th percentile)
                ys_ϵff = Vector{Float64}(undef, grid.n_bins)
                # ϵff 25th percentile
                ys_low_ϵff = Vector{Float64}(undef, grid.n_bins)
                # ϵff 75th percentile
                ys_high_ϵff = Vector{Float64}(undef, grid.n_bins)

                for i in eachindex(binned_ϵffs)

                    if isempty(binned_ϵffs[i])

                        ys_ϵff[i]      = NaN
                        ys_low_ϵff[i]  = NaN
                        ys_high_ϵff[i] = NaN

                    else

                        ys_ϵff[i]      = median(binned_ϵffs[i])
                        ys_low_ϵff[i]  = quantile(binned_ϵffs[i], 0.25)
                        ys_high_ϵff[i] = quantile(binned_ϵffs[i], 0.75)

                    end

                end

                lines!(ax_02, log10.(x_ngas), ys_ϵff)

                band!(ax_02, log10.(x_ngas), ys_low_ϵff, ys_high_ϵff)

            end

            rm(temp_folder; recursive=true)

            ##################################################
            # Gas efficiency per free-fall vs gas metallicity
            ##################################################

            temp_folder = joinpath(output_path, "_eff_vs_Z")

            plotSnapshot(
                [simulation],
                request,
                [scatter!];
                pf_kwargs=[(;)],
                output_path=temp_folder,
                base_filename="eff_vs_Z",
                slice=128,
                filter_function,
                da_functions=[GalaxyInspector.daScatterGalaxy],
                da_args=[(:gas_metallicity, :gas_eff)],
                da_kwargs=[(; filter_function=extra_filter)],
                transform_box=true,
                translation,
                rotation,
                x_unit=Z_plot_params.unit,
                y_unit=ϵff_plot_params.unit,
                save_figures=false,
                backup_results=true,
            )

            jld2_path = joinpath(temp_folder, "eff_vs_Z.jld2")

            jldopen(jld2_path, "r") do jld2_file

                Z_gas, ϵ_ff = jld2_file["snap_127"][basename(simulation)]

                # Remove point where Z or ϵff are zero or NaN
                zero_filter = (x, y) -> iszero(x) || iszero(y) || isnan(x) || isnan(y)
                delete_idxs = zero_filter.(Z_gas, ϵ_ff)
                deleteat!(Z_gas, delete_idxs)
                deleteat!(ϵ_ff, delete_idxs)

                # Set up a logarithmic grid for the x axis
                grid = GalaxyInspector.LinearGrid(extrema(Z_gas)..., 40; log=true)

                # Set up a logarithmic scale for the y axis
                log_ϵ_ff = log10.(ϵ_ff)

                binned_ϵffs = GalaxyInspector.listHistogram1D(Z_gas, log_ϵ_ff, grid)

                # Gas metallicity
                x_Zgas      = grid.x_axis
                # Median of ϵff (50th percentile)
                ys_ϵff      = Vector{Float64}(undef, grid.n_bins)
                # ϵff 25th percentile
                ys_low_ϵff  = Vector{Float64}(undef, grid.n_bins)
                # ϵff 75th percentile
                ys_high_ϵff = Vector{Float64}(undef, grid.n_bins)

                for i in eachindex(binned_ϵffs)

                    if isempty(binned_ϵffs[i])

                        ys_ϵff[i]      = NaN
                        ys_low_ϵff[i]  = NaN
                        ys_high_ϵff[i] = NaN

                    else

                        ys_ϵff[i]      = median(binned_ϵffs[i])
                        ys_low_ϵff[i]  = quantile(binned_ϵffs[i], 0.25)
                        ys_high_ϵff[i] = quantile(binned_ϵffs[i], 0.75)

                    end

                end

                lines!(ax_03, log10.(x_Zgas), ys_ϵff)

                band!(ax_03, log10.(x_Zgas), ys_low_ϵff, ys_high_ϵff)

            end

            rm(temp_folder; recursive=true)

        end

        axislegend(ax_01, position=:lt, framevisible=false, nbanks=1)

        Makie.save(joinpath(output_path, "eff_vs_Z_n_t.png"), f)

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

    ################################################################################################
    # Basic analysis
    ################################################################################################

    SIMULATIONS = joinpath.(
        "/run/media/elozano/DISK4T/simulations/",
        [
            "SFM_01",  # Paper I
            "SFM_02",  # + Metals and dust
            "SFM_03",  # + Shielding
            "SFM_04",  # + Photon attenuation
            "SFM_05",  # + UVB
            "SFM_06",  # + LWB = Paper II
        ]
    )

    LABELS = [
        "Paper I",
        "+ Z & D",
        "+ S",
        "+ PA",
        "+ UVB",
        "+ LWB = Paper II",
    ]

    lozano_2026(SIMULATIONS, BASE_OUT_PATH, LOGGING; labels=LABELS)

end
