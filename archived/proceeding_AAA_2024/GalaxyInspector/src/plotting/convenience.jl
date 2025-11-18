####################################################################################################
# Opinionated convenience functions
####################################################################################################

"""
    snapshotReport(
        simulation_paths::Vector{String},
        slices::Vector{Int};
        <keyword arguments>
    )::Nothing

Write a text file with information about a given snapshot.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. One text file will be printed for each simulation.
  - `slices::Vector{Int}`: Selects which snapshots to plot for each simulation, starts at 1 and is independent of the number in the file name. If every snapshot is present, the relation is `slice_n` = filename_number + 1.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be considered in the "filtered" section of the report. The options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `halo_idx::Int=1`: Index of the target halo (FoF group) for the corresponding section. Starts at 1.
  - `subhalo_rel_idx::Int=1`: Index of the target subhalo (subfind), relative to the target halo, for the corresponding section. Starts at 1.
"""
function snapshotReport(
    simulation_paths::Vector{String},
    slices::Vector{Int};
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    halo_idx::Int=1,
    subhalo_rel_idx::Int=1,
)::Nothing

    for (i, simulation_path) in pairs(simulation_paths)

        ############################################################################################
        # Load the relevant values and check for missing files
        ############################################################################################

        # Make a dataframe with the following columns:
        #   - 1. DataFrame index
        #   - 2. Number in the file name
        #   - 3. Scale factor
        #   - 4. Redshift
        #   - 5. Physical time
        #   - 6. Lookback time
        #   - 7. Snapshot path
        #   - 8. Group catalog path
        simulation_table = makeSimulationTable(simulation_path)

        # Select the target slice for the `i`-th simulation
        slice_n = ring(slices, i)

        # Get the number in the filename
        snap_number = safeSelect(simulation_table[!, :numbers], slice_n)

        # Check that after slicing there is one snapshot left
        (
            !isempty(snap_number) ||
            throw(ArgumentError("snapshotReport: There are no snapshot number \
            $(slice_n), the contents of $(simulation_path) are: \n$(simulation_table)"))
        )

        # Find the target row
        snapshot_row = filter(:numbers => ==(lpad(snap_number, 3, "0")), simulation_table)

        # Read the path to the target snapshot
        snapshot_path = snapshot_row[1, :snapshot_paths]

        # Check that the snapshot path is not missing
        snapshot_filename = "\"$(SNAP_BASENAME)_$(lpad(snap_number, 3, "0"))\""
        (
            !ismissing(snapshot_path) ||
            throw(ArgumentError("snapshotReport: The snapshot $(snapshot_filename) is missing \
            in $(simulation_path)"))
        )

        # Read the path to the target group catalog file
        groupcat_path = snapshot_row[1, :groupcat_paths]

        # Check if the simulation is cosmological
        cosmological = isCosmological(snapshot_path)

        # Read the physical time since the Big Bang
        physical_time = round(ustrip(u"Gyr", snapshot_row[1, :physical_times]), digits=2)

        # Select the ordinal index of the target snapshot
        o_idx = snapshot_row[1, :ids]

        # Compute the number of snapshots in the folder
        snapshot_length = count(!ismissing, simulation_table[!, :snapshot_paths])

        # Compute the number of group catalog files in the folder
        groupcat_length = count(!ismissing, simulation_table[!, :groupcat_paths])

        # Read the header of the target snapshot
        snapshot_header = readSnapHeader(snapshot_path)

        ############################################################################################
        # Print the report header
        ############################################################################################

        # Create the output file
        filename = "$(SNAP_BASENAME)_$(lpad(snap_number, 3, "0"))_of_$(basename(simulation_path))"
        file = open(joinpath(mkpath(output_path), "report_for_$(filename).txt"), "w")

        println(file, "#"^100)
        println(file, "\nSimulation name:  $(basename(simulation_path))")
        println(file, "Physical time:    $(physical_time) Gyr")

        if cosmological
            # For cosmological simulations print the scale factor and the redshift
            scale_factor = round(snapshot_row[1, :scale_factors], digits=3)
            redshift = round(snapshot_row[1, :redshifts], digits=3)

            println(file, "Scale factor:     $(scale_factor)")
            println(file, "Redshift:         $(redshift)")
            println(file, "Cosmological:     Yes")
        else
            println(file, "Cosmological:     No")
        end

        if !PHYSICAL_UNITS && cosmological
            println(file, "Report units:     Comoving\n")
        else
            println(file, "Report units:     Physical\n")
        end

        println(file, "#"^100)
        println(file, "\nSnapshot number:  $(o_idx) (of $(snapshot_length))")
        println(file, "Snapshot path:    $(snapshot_path)\n")

        if !ismissing(groupcat_path)

            println(file, "Subfind number:   $(o_idx) (of $(groupcat_length))")
            println(file, "Subfind path:     $(groupcat_path)\n")

        end

        println(file, "#"^100)

        ############################################################################################
        # Read the data in the snapshot and group catalog file
        ############################################################################################

        # Select one snapshot file
        if isfile(snapshot_path)
            file_path = snapshot_path
        else
            file_path = minimum(glob("$(SNAP_BASENAME)_*.*.hdf5", snapshot_path))
        end

        physical_components = [:gas, :halo, :stars, :black_hole]

        # Detect which of the main physical components are present in the snapshot
        component_list = h5open(file_path, "r") do snapshot
            filter(
                in(physical_components),
                [get(PARTICLE_TYPE, key, nothing) for key in keys(snapshot)],
            )
        end

        # Select the filter function, translation, rotation, and request dictionary
        filter_function, translation, rotation, request = selectFilter(
            filter_mode,
            mergeRequests(
                Dict(component => ["POS ", "MASS", "VEL "] for component in component_list),
                Dict(
                    :gas => [
                        "NHP ",
                        "NH  ",
                        "PRES",
                        "FRAC",
                        "ODIT",
                        "TAUS",
                        "RHO ",
                        "RHOC",
                        "PARZ",
                        "ETAD",
                        "ETAI",
                        "PARR",
                        "ID  ",
                        "SFFL",
                    ],
                    :stars => [
                        "ODIT",
                        "PARA",
                        "TAUS",
                        "RHOC",
                        "PARZ",
                        "ETAD",
                        "ETAI",
                        "PARR",
                        "FRAC",
                        "GMAS",
                        "GSFR",
                        "GPRE",
                        "ID  ",
                    ],
                ),
            ),
        )

        # Read the snapshot data
        if !in(filter_mode, [:all, :sphere])

            # Check that the group catalog data is available
            if !ismissing(groupcat_path) && isSubfindActive(groupcat_path)

                # Create the data dictionary
                data_dict = makeDataDict(simulation_path, slice_n, request)

            else

                throw(ArgumentError("snapshotReport: You asked for a filter base on the \
                halos/subhalos or you gave a personalized filter, but I could not find a valid \
                group catalog file"))

            end

        else

            # Create the data dictionary
            data_dict = readSnapshot(snapshot_path, request)

        end

        ############################################################################################
        # Print the global properties of the simulation
        ############################################################################################

        println(file, "\nGlobal properties (full simulation box):")

        ##############################################################
        # Print the total number of cells/particles of each component
        ##############################################################

        println(file, "\n\tCell/particle number (full simulation box):\n")

        total_count = 0
        for component in component_list

            count = getindex(getfield(snapshot_header, :num_total), PARTICLE_INDEX[component] + 1)

            total_count += count

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(25 - length(title))

            println(file, "\t\t$(title)$(count)")

        end

        println(file, "\n\t\tTotal count:             $(total_count)\n")

        #########################################
        # Print the total mass of each component
        #########################################

        println(file, "\tMasses (full simulation box):\n")

        total_mass = 0.0u"Msun"
        for component in component_list

            mass = sum(data_dict[component]["MASS"]; init=0.0u"Msun")

            total_mass += mass

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(25 - length(title))

            println(file, "\t\t$(title)$(round(typeof(1.0u"Msun"), mass, sigdigits=3))")

        end

        println(
            file,
            "\n\t\tTotal mass:              $(round(typeof(1.0u"Msun"), total_mass, sigdigits=3))\n",
        )

        ########################################
        # Print the mass of each hydrogen phase
        ########################################

        if :gas in component_list

            println(file, "\tHydrogen masses (full simulation box):\n")

            gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

            hii_mass = sum(computeMass(data_dict, :ionized); init=0.0u"Msun")
            hii_percent = round((hii_mass / gas_mass) * 100, sigdigits=3)

            hi_mass = sum(computeMass(data_dict, :atomic); init=0.0u"Msun")
            hi_percent = round((hi_mass / gas_mass) * 100, sigdigits=3)

            h2_mass = sum(computeMass(data_dict, :molecular); init=0.0u"Msun")
            h2_percent = round((h2_mass / gas_mass) * 100, sigdigits=3)

            h2p_mass = sum(computeMass(data_dict, :br_molecular); init=0.0u"Msun")
            h2p_percent = round((h2p_mass / gas_mass) * 100, sigdigits=3)

            hn_mass = sum(computeMass(data_dict, :neutral); init=0.0u"Msun")
            hn_percent = round((hn_mass / gas_mass) * 100, sigdigits=3)

            if !iszero(hii_mass)
                title = "Ionized mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), hii_mass, sigdigits=3)) \
                    ($(hii_percent)% of total gas mass)",
                )
            end

            if !iszero(hi_mass)
                title = "Atomic mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), hi_mass, sigdigits=3)) \
                    ($(hi_percent)% of total gas mass)",
                )
            end

            if !iszero(h2_mass)
                title = "Molecular mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), h2_mass, sigdigits=3)) \
                    ($(h2_percent)% of total gas mass)",
                )
            end

            if !iszero(h2p_mass)
                title = "Molecular mass (BR):"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), h2p_mass, sigdigits=3)) \
                    ($(h2p_percent)% of total gas mass)",
                )
            end

            if !iszero(hn_mass)
                title = "Neutral mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), hn_mass, sigdigits=3)) \
                    ($(hn_percent)% of total gas mass)\n",
                )
            end

        end

        ############################################################################################
        # Filter the simulation box
        ############################################################################################

        filterData!(data_dict; filter_function)

        println(file, "#"^100)
        println(file, "Filtered box with:")

        if filter_mode isa Symbol
            println(file, "\n\tFilter mode: $(filter_mode)")
        else
            println(file, "\n\tFilter function: $(String(Symbol(filter_function)))")
        end

        println(file, "\tTranslation: $(translation)")
        println(file, "\tRotation: $(rotation)")
        println(file, "#"^100)

        ############################################################################################
        # Print the global properties of the simulation after filtering
        ############################################################################################

        println(file, "\nGlobal properties (filtered box):\n")

        ########################################################
        # Print the number of cells/particles of each component
        ########################################################

        println(file, "\tCell/particle number (filtered box):\n")

        total_count = 0
        for component in component_list

            count = length(data_dict[component]["MASS"])

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(25 - length(title))

            total_count += count

            println(file, "\t\t$(title)$(count)")

        end

        println(file, "\n\t\tTotal count:             $(total_count)\n")

        ###################################
        # Print the mass of each component
        ###################################

        println(file, "\tMasses (filtered box):\n")

        total_mass = 0.0u"Msun"
        for component in component_list

            mass = sum(data_dict[component]["MASS"]; init=0.0u"Msun")

            total_mass += mass

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(25 - length(title))

            println(file, "\t\t$(title)$(round(typeof(1.0u"Msun"), mass, sigdigits=3))")

        end

        println(
            file,
            "\n\t\tTotal mass:              $(round(typeof(1.0u"Msun"), total_mass, sigdigits=3))\n",
        )

        ########################################
        # Print the mass of each hydrogen phase
        ########################################

        if :gas in component_list

            println(file, "\tHydrogen masses (filtered box):\n")

            gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

            hii_mass = sum(computeMass(data_dict, :ionized); init=0.0u"Msun")
            hii_percent = round((hii_mass / gas_mass) * 100, sigdigits=3)

            hi_mass = sum(computeMass(data_dict, :atomic); init=0.0u"Msun")
            hi_percent = round((hi_mass / gas_mass) * 100, sigdigits=3)

            h2_mass = sum(computeMass(data_dict, :molecular); init=0.0u"Msun")
            h2_percent = round((h2_mass / gas_mass) * 100, sigdigits=3)

            h2p_mass = sum(computeMass(data_dict, :br_molecular); init=0.0u"Msun")
            h2p_percent = round((h2p_mass / gas_mass) * 100, sigdigits=3)

            hn_mass = sum(computeMass(data_dict, :neutral); init=0.0u"Msun")
            hn_percent = round((hn_mass / gas_mass) * 100, sigdigits=3)

            if !iszero(hii_mass)
                title = "Ionized mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), hii_mass, sigdigits=3)) \
                    ($(hii_percent)% of total gas mass)",
                )
            end

            if !iszero(hi_mass)
                title = "Atomic mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), hi_mass, sigdigits=3)) \
                    ($(hi_percent)% of total gas mass)",
                )
            end

            if !iszero(h2_mass)
                title = "Molecular mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), h2_mass, sigdigits=3)) \
                    ($(h2_percent)% of total gas mass)",
                )
            end

            if !iszero(h2p_mass)
                title = "Molecular mass (BR):"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), h2p_mass, sigdigits=3)) \
                    ($(h2p_percent)% of total gas mass)",
                )
            end

            if !iszero(hn_mass)
                title = "Neutral mass:"
                title *= " "^(25 - length(title))
                println(
                    file,
                    "\t\t$(title)$(round(typeof(1.0u"Msun"), hn_mass, sigdigits=3)) \
                    ($(hn_percent)% of total gas mass)\n",
                )
            end

        end

        #############################################
        # Print the center of mass of each component
        #############################################

        println(file, "\tCenter of mass (filtered box):\n")

        global_cm = computeGlobalCenterOfMass(data_dict)

        for component in component_list

            cm = computeCenterOfMass(data_dict[component]["POS "], data_dict[component]["MASS"])

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(25 - length(title))

            println(file, "\t\t$(title)$(round.(ustrip.(u"Mpc", cm), sigdigits=6)) $(u"Mpc")")
            println(file, "\t\tDistance to global CM:   $(sqrt(sum((global_cm - cm).^2)))\n")

        end

        global_cm = round.(ustrip.(u"Mpc", global_cm), sigdigits=6)
        println(file, "\t\tGlobal center of mass:   $(global_cm) $(u"Mpc")\n")

        #################################################################
        # Print the fraction of gas cells that have enter the SF routine
        #################################################################

        if !isempty(data_dict[:gas]["SFFL"])

            println(
                file,
                "\tFraction of gas cells that have enter the SF routine (filtered box):\n"
            )

            gas_masses = computeMass(data_dict, :gas)

            total_number = length(gas_masses)
            stellar_gas_number = count(isone, data_dict[:gas]["SFFL"])
            fraction = (stellar_gas_number / total_number) * 100

            idxs = findall(isone, data_dict[:gas]["SFFL"])
            stellar_gas_mass = sum(gas_masses[idxs])
            mass_fraction = (stellar_gas_mass / sum(gas_masses)) * 100

            println(file, "\t\t$(round(fraction, sigdigits=3))% of the cells")
            println(file, "\t\t$(round(mass_fraction, sigdigits=3))% of the mass\n")

        end

        ###############################################
        # Print the properties of the star forming gas
        ###############################################

        if any(
            !isempty,
            [
                data_dict[:stars]["ODIT"],
                data_dict[:stars]["PARA"],
                data_dict[:stars]["TAUS"],
                data_dict[:stars]["RHOC"],
                data_dict[:stars]["PARZ"],
                data_dict[:stars]["ETAD"],
                data_dict[:stars]["ETAI"],
                data_dict[:stars]["PARR"],
                data_dict[:stars]["FRAC"],
                data_dict[:stars]["GMAS"],
                data_dict[:stars]["GSFR"],
                data_dict[:stars]["GPRE"],
            ],
        )

            println(file, "\tProperties of the gas that has formed stars (filtered box):\n")

        end

        if !isempty(data_dict[:stars]["ODIT"])

            odit = ustrip.(u"Myr", data_dict[:stars]["ODIT"])

            println(file, "\t\tTotal integration time:\n")
            println(file, "\t\t\tMean:    $(round(mean(odit), sigdigits=4)) Myr")
            println(file, "\t\t\tMedian:  $(round(median(odit), sigdigits=4)) Myr")
            println(file, "\t\t\tMode:    $(round(mode(odit)[1], sigdigits=4)) Myr")
            println(file, "\t\t\tMinimum: $(round(minimum(odit), sigdigits=4)) Myr")
            println(file, "\t\t\tMaximum: $(round(maximum(odit), sigdigits=4)) Myr\n")

            odit_50 = computeMassQty(odit, data_dict[:stars]["MASS"]; percent=50.0)
            odit_90 = computeMassQty(odit, data_dict[:stars]["MASS"]; percent=90.0)
            odit_95 = computeMassQty(odit, data_dict[:stars]["MASS"]; percent=95.0)

            println(file, "\t\t\tTotal integration time enclosing X% of the stellar mass:\n")
            println(file, "\t\t\t\t$(round(odit_50, sigdigits=4)) Myr (50%)")
            println(file, "\t\t\t\t$(round(odit_90, sigdigits=4)) Myr (90%)")
            println(file, "\t\t\t\t$(round(odit_95, sigdigits=4)) Myr (95%)\n")

            odit_percent_low = computeMassFraction(
                odit,
                data_dict[:stars]["MASS"],
                (0.0, 100.0), # Myr
            ) * 100

            println(file, "\t\t\tFraction of stellar mass with a total integration time < 100 Myr:\n")
            println(file, "\t\t\t\t$(round(odit_percent_low, sigdigits=4))%\n")

            odit_percent_high = computeMassFraction(
                odit,
                data_dict[:stars]["MASS"],
                (100.0, Inf), # Myr
            ) * 100

            println(file, "\t\t\tFraction of stellar mass with a total integration time > 100 Myr:\n")
            println(file, "\t\t\t\t$(round(odit_percent_high, sigdigits=4))%\n")

        end

        if !isempty(data_dict[:stars]["PARA"])

            para = data_dict[:stars]["PARA"]

            println(file, "\t\tScale factor:\n")
            println(file, "\t\t\tMean:    $(round(mean(para), sigdigits=4))")
            println(file, "\t\t\tMedian:  $(round(median(para), sigdigits=4))")
            println(file, "\t\t\tMode:    $(round(mode(para)[1], sigdigits=4))")
            println(file, "\t\t\tMinimum: $(round(minimum(para), sigdigits=4))")
            println(file, "\t\t\tMaximum: $(round(maximum(para), sigdigits=4))\n")

            para_50 = computeMassQty(para, data_dict[:stars]["MASS"]; percent=50.0)
            para_90 = computeMassQty(para, data_dict[:stars]["MASS"]; percent=90.0)
            para_95 = computeMassQty(para, data_dict[:stars]["MASS"]; percent=95.0)

            println(file, "\t\t\tScale factor enclosing X% of the stellar mass:\n")
            println(file, "\t\t\t\t$(round(para_50, sigdigits=4)) (50%)")
            println(file, "\t\t\t\t$(round(para_90, sigdigits=4)) (90%)")
            println(file, "\t\t\t\t$(round(para_95, sigdigits=4)) (95%)\n")

        end

        if !isempty(data_dict[:stars]["TAUS"])

            τS = ustrip.(u"Myr", data_dict[:stars]["TAUS"])

            println(file, "\t\tStar formation time parameter (τS):\n")
            println(file, "\t\t\tMean:    $(round(mean(τS), sigdigits=4)) Myr")
            println(file, "\t\t\tMedian:  $(round(median(τS), sigdigits=4)) Myr")
            println(file, "\t\t\tMode:    $(round(mode(τS)[1], sigdigits=4)) Myr")
            println(file, "\t\t\tMinimum: $(round(minimum(τS), sigdigits=4)) Myr")
            println(file, "\t\t\tMaximum: $(round(maximum(τS), sigdigits=4)) Myr\n")

            τS_50 = computeMassQty(τS, data_dict[:stars]["MASS"]; percent=50.0)
            τS_90 = computeMassQty(τS, data_dict[:stars]["MASS"]; percent=90.0)
            τS_95 = computeMassQty(τS, data_dict[:stars]["MASS"]; percent=95.0)

            println(file, "\t\t\tτS enclosing X% of the stellar mass:\n")
            println(file, "\t\t\t\t$(round(τS_50, sigdigits=4)) Myr (50%)")
            println(file, "\t\t\t\t$(round(τS_90, sigdigits=4)) Myr (90%)")
            println(file, "\t\t\t\t$(round(τS_95, sigdigits=4)) Myr (95%)\n")

        end

        if !isempty(data_dict[:stars]["RHOC"])

            rhoc = ustrip.(u"cm^-3", data_dict[:stars]["RHOC"])

            println(file, "\t\tCell density:\n")
            println(file, "\t\t\tMean:    $(round(mean(rhoc), sigdigits=4)) cm^-3")
            println(file, "\t\t\tMedian:  $(round(median(rhoc), sigdigits=4)) cm^-3")
            println(file, "\t\t\tMode:    $(round(mode(rhoc)[1], sigdigits=4)) cm^-3")
            println(file, "\t\t\tMinimum: $(round(minimum(rhoc), sigdigits=4)) cm^-3")
            println(file, "\t\t\tMaximum: $(round(maximum(rhoc), sigdigits=4)) cm^-3\n")

            rhoc_50 = computeMassQty(rhoc, data_dict[:stars]["MASS"]; percent=50.0)
            rhoc_90 = computeMassQty(rhoc, data_dict[:stars]["MASS"]; percent=90.0)
            rhoc_95 = computeMassQty(rhoc, data_dict[:stars]["MASS"]; percent=95.0)

            println(file, "\t\t\tCell density enclosing X% of the stellar mass:\n")
            println(file, "\t\t\t\t$(round(rhoc_50, sigdigits=4)) cm^-3 (50%)")
            println(file, "\t\t\t\t$(round(rhoc_90, sigdigits=4)) cm^-3 (90%)")
            println(file, "\t\t\t\t$(round(rhoc_95, sigdigits=4)) cm^-3 (95%)\n")

        end

        if !isempty(data_dict[:stars]["PARZ"])

            parz = data_dict[:stars]["PARZ"] ./ SOLAR_METALLICITY

            println(file, "\t\tMetallicity:\n")
            println(file, "\t\t\tMean:    $(round(mean(parz), sigdigits=4)) Z⊙")
            println(file, "\t\t\tMedian:  $(round(median(parz), sigdigits=4)) Z⊙")
            println(file, "\t\t\tMode:    $(round(mode(parz)[1], sigdigits=4)) Z⊙")
            println(file, "\t\t\tMinimum: $(round(minimum(parz), sigdigits=4)) Z⊙")
            println(file, "\t\t\tMaximum: $(round(maximum(parz), sigdigits=4)) Z⊙\n")

            parz_50 = computeMassQty(parz, data_dict[:stars]["MASS"]; percent=50.0)
            parz_90 = computeMassQty(parz, data_dict[:stars]["MASS"]; percent=90.0)
            parz_95 = computeMassQty(parz, data_dict[:stars]["MASS"]; percent=95.0)

            println(file, "\t\t\tMetallicity enclosing X% of the stellar mass:\n")
            println(file, "\t\t\t\t$(round(parz_50, sigdigits=4)) Z⊙ (50%)")
            println(file, "\t\t\t\t$(round(parz_90, sigdigits=4)) Z⊙ (90%)")
            println(file, "\t\t\t\t$(round(parz_95, sigdigits=4)) Z⊙ (95%)\n")

        end

        if !isempty(data_dict[:stars]["ETAD"])

            ηd = data_dict[:stars]["ETAD"]

            println(file, "\t\tPhotodissociation parameter (ηd):\n")
            println(file, "\t\t\tMean:    $(round(mean(ηd), sigdigits=4))")
            println(file, "\t\t\tMedian:  $(round(median(ηd), sigdigits=4))")
            println(file, "\t\t\tMode:    $(round(mode(ηd)[1], sigdigits=4))")
            println(file, "\t\t\tMinimum: $(round(minimum(ηd), sigdigits=4))")
            println(file, "\t\t\tMaximum: $(round(maximum(ηd), sigdigits=4))\n")

            ηd_50 = computeMassQty(ηd, data_dict[:stars]["MASS"]; percent=50.0)
            ηd_90 = computeMassQty(ηd, data_dict[:stars]["MASS"]; percent=90.0)
            ηd_95 = computeMassQty(ηd, data_dict[:stars]["MASS"]; percent=95.0)

            println(file, "\t\t\tηd enclosing X% of the stellar mass:\n")
            println(file, "\t\t\t\t$(round(ηd_50, sigdigits=4)) (50%)")
            println(file, "\t\t\t\t$(round(ηd_90, sigdigits=4)) (90%)")
            println(file, "\t\t\t\t$(round(ηd_95, sigdigits=4)) (95%)\n")

        end

        if !isempty(data_dict[:stars]["ETAI"])

            ηi = data_dict[:stars]["ETAI"]

            println(file, "\t\tPhotoionization parameter (ηi):\n")
            println(file, "\t\t\tMean:    $(round(mean(ηi), sigdigits=4))")
            println(file, "\t\t\tMedian:  $(round(median(ηi), sigdigits=4))")
            println(file, "\t\t\tMode:    $(round(mode(ηi)[1], sigdigits=4))")
            println(file, "\t\t\tMinimum: $(round(minimum(ηi), sigdigits=4))")
            println(file, "\t\t\tMaximum: $(round(maximum(ηi), sigdigits=4))\n")

            ηi_50 = computeMassQty(ηi, data_dict[:stars]["MASS"]; percent=50.0)
            ηi_90 = computeMassQty(ηi, data_dict[:stars]["MASS"]; percent=90.0)
            ηi_95 = computeMassQty(ηi, data_dict[:stars]["MASS"]; percent=95.0)

            println(file, "\t\t\tηi enclosing X% of the stellar mass:\n")
            println(file, "\t\t\t\t$(round(ηi_50, sigdigits=4)) (50%)")
            println(file, "\t\t\t\t$(round(ηi_90, sigdigits=4)) (90%)")
            println(file, "\t\t\t\t$(round(ηi_95, sigdigits=4)) (95%)\n")

        end

        if !isempty(data_dict[:stars]["PARR"])

            R = data_dict[:stars]["PARR"]

            println(file, "\t\tMass recycling parameter (R):\n")
            println(file, "\t\t\tMean:    $(round(mean(R), sigdigits=4))")
            println(file, "\t\t\tMedian:  $(round(median(R), sigdigits=4))")
            println(file, "\t\t\tMode:    $(round(mode(R)[1], sigdigits=4))")
            println(file, "\t\t\tMinimum: $(round(minimum(R), sigdigits=4))")
            println(file, "\t\t\tMaximum: $(round(maximum(R), sigdigits=4))\n")

            R_50 = computeMassQty(R, data_dict[:stars]["MASS"]; percent=50.0)
            R_90 = computeMassQty(R, data_dict[:stars]["MASS"]; percent=90.0)
            R_95 = computeMassQty(R, data_dict[:stars]["MASS"]; percent=95.0)

            println(file, "\t\t\tR enclosing X% of the stellar mass:\n")
            println(file, "\t\t\t\t$(round(R_50, sigdigits=4)) (50%)")
            println(file, "\t\t\t\t$(round(R_90, sigdigits=4)) (90%)")
            println(file, "\t\t\t\t$(round(R_95, sigdigits=4)) (95%)\n")

        end

        if !isempty(data_dict[:stars]["GMAS"])

            gmas = ustrip.(u"Msun", data_dict[:stars]["GMAS"]) ./ exp10(4.0)

            println(file, "\t\tParent gas mass:\n")
            println(file, "\t\t\tMean:    $(round(mean(gmas), sigdigits=4)) × 10⁴ M⊙")
            println(file, "\t\t\tMedian:  $(round(median(gmas), sigdigits=4)) × 10⁴ M⊙")
            println(file, "\t\t\tMode:    $(round(mode(gmas)[1], sigdigits=4)) × 10⁴ M⊙")
            println(file, "\t\t\tMinimum: $(round(minimum(gmas), sigdigits=4)) × 10⁴ M⊙")
            println(file, "\t\t\tMaximum: $(round(maximum(gmas), sigdigits=4)) × 10⁴ M⊙\n")

            gmas_50 = computeMassQty(gmas, data_dict[:stars]["MASS"]; percent=50.0)
            gmas_90 = computeMassQty(gmas, data_dict[:stars]["MASS"]; percent=90.0)
            gmas_95 = computeMassQty(gmas, data_dict[:stars]["MASS"]; percent=95.0)

            println(file, "\t\t\tParent gas mass enclosing X% of the stellar mass:\n")
            println(file, "\t\t\t\t$(round(gmas_50, sigdigits=4)) × 10⁴ M⊙ (50%)")
            println(file, "\t\t\t\t$(round(gmas_90, sigdigits=4)) × 10⁴ M⊙ (90%)")
            println(file, "\t\t\t\t$(round(gmas_95, sigdigits=4)) × 10⁴ M⊙ (95%)\n")

        end

        if !isempty(data_dict[:stars]["GSFR"])

            gsfr = ustrip.(u"Msun*yr^-1", data_dict[:stars]["GSFR"])

            println(file, "\t\tParent SFR:\n")
            println(file, "\t\t\tMean:    $(round(mean(gsfr), sigdigits=4)) M⊙ yr^-1")
            println(file, "\t\t\tMedian:  $(round(median(gsfr), sigdigits=4)) M⊙ yr^-1")
            println(file, "\t\t\tMode:    $(round(mode(gsfr)[1], sigdigits=4)) M⊙ yr^-1")
            println(file, "\t\t\tMinimum: $(round(minimum(gsfr), sigdigits=4)) M⊙ yr^-1")
            println(file, "\t\t\tMaximum: $(round(maximum(gsfr), sigdigits=4)) M⊙ yr^-1\n")

            gsfr_50 = computeMassQty(gsfr, data_dict[:stars]["MASS"]; percent=50.0)
            gsfr_90 = computeMassQty(gsfr, data_dict[:stars]["MASS"]; percent=90.0)
            gsfr_95 = computeMassQty(gsfr, data_dict[:stars]["MASS"]; percent=95.0)

            println(file, "\t\t\tParent SFR enclosing X% of the stellar mass:\n")
            println(file, "\t\t\t\t$(round(gsfr_50, sigdigits=4)) M⊙ yr^-1 (50%)")
            println(file, "\t\t\t\t$(round(gsfr_90, sigdigits=4)) M⊙ yr^-1 (90%)")
            println(file, "\t\t\t\t$(round(gsfr_95, sigdigits=4)) M⊙ yr^-1 (95%)\n")

        end

        if !isempty(data_dict[:stars]["GPRE"])

            gpre = ustrip.(u"dyn*cm^-2", data_dict[:stars]["GPRE"])

            println(file, "\t\tParent gas pressure:\n")
            println(file, "\t\t\tMean:    $(round(mean(gpre), sigdigits=4)) dyn cm^-2")
            println(file, "\t\t\tMedian:  $(round(median(gpre), sigdigits=4)) dyn cm^-2")
            println(file, "\t\t\tMode:    $(round(mode(gpre)[1], sigdigits=4)) dyn cm^-2")
            println(file, "\t\t\tMinimum: $(round(minimum(gpre), sigdigits=4)) dyn cm^-2")
            println(file, "\t\t\tMaximum: $(round(maximum(gpre), sigdigits=4)) dyn cm^-2\n")

            gpre_50 = computeMassQty(gpre, data_dict[:stars]["MASS"]; percent=50.0)
            gpre_90 = computeMassQty(gpre, data_dict[:stars]["MASS"]; percent=90.0)
            gpre_95 = computeMassQty(gpre, data_dict[:stars]["MASS"]; percent=95.0)

            println(file, "\t\t\tParent gas pressure enclosing X% of the stellar mass:\n")
            println(file, "\t\t\t\t$(round(gpre_50, sigdigits=4)) dyn cm^-2 (50%)")
            println(file, "\t\t\t\t$(round(gpre_90, sigdigits=4)) dyn cm^-2 (90%)")
            println(file, "\t\t\t\t$(round(gpre_95, sigdigits=4)) dyn cm^-2 (95%)\n")

        end

        ############################################################
        # Print the values of the ODEs parameters for the gas cells
        ############################################################

        quantities = ["ODIT", "TAUS", "RHOC", "PARZ", "ETAD", "ETAI", "PARR"]
        names = [ "integration time", "τS", "cell density", "metallicity", "ηd", "ηi", "R"]
        units = [
            u"Myr",
            u"Myr",
            u"cm^-3",
            Unitful.NoUnits,
            Unitful.NoUnits,
            Unitful.NoUnits,
            Unitful.NoUnits,
        ]

        if any(isBlockPresent.(:gas, quantities, snapshot_path))
            println(file, "\tODE parameters for the gas cells (filtered box):\n")
        end

        for (quantity, unit, name) in zip(quantities, units, names)

            if isBlockPresent(:gas, quantity, snapshot_path)

                values = filter(!isnan, ustrip.(unit, data_dict[:gas][quantity]))

                isempty(values) && continue

                println(file, "\t\tMean $name:    $(round(mean(values), sigdigits=4)) $unit")
                println(file, "\t\tMedian $name:  $(round(median(values), sigdigits=4)) $unit")
                println(file, "\t\tMode $name:    $(round(mode(values)[1], sigdigits=4)) $unit")
                println(file, "\t\tMinimum $name: $(round(minimum(values), sigdigits=4)) $unit")
                println(file, "\t\tMaximum $name: $(round(maximum(values), sigdigits=4)) $unit\n")

            end

        end

        ###############################
        # Translate the simulation box
        ###############################

        translateData!(data_dict, translation)

        ##########################################################
        # Print the normalized angular momentum of each component
        ##########################################################

        println(file, "\tNormalized angular momentum (filtered box):\n")

        for component in component_list

            L = computeTotalAngularMomentum(
                data_dict[component]["POS "],
                data_dict[component]["VEL "],
                data_dict[component]["MASS"];
            )

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(25 - length(title))

            println(file, "\t\t$(title)$(round.(L, sigdigits=3))")

        end

        global_L = round.(computeGlobalAngularMomentum(data_dict), sigdigits=3)

        println(file, "\n\t\tGlobal angular momentum: $(global_L)\n")

        #############################################
        # Print the spin parameter of each component
        #############################################

        println(file, "\tSpin parameter (R = $(DISK_R) - filtered box):\n")

        for component in component_list

            λ = computeSpinParameter(
                data_dict[component]["POS "],
                data_dict[component]["VEL "],
                data_dict[component]["MASS"],
            )

            title = "$(PARTICLE_NAMES[component]):"
            title *= " "^(25 - length(title))

            println(file, "\t\t$(title)$(round.(λ, sigdigits=3))")

        end

        global_λ = round.(computeGlobalSpinParameter(data_dict), sigdigits=3)

        println(file, "\n\t\tTotal spin parameter:    $(global_λ)\n")

        ############################
        # Rotate the simulation box
        ############################

        rotateData!(data_dict, rotation)

        ###################################################################################
        # Print the total height of a cylinder, of infinite radius, containing 90% and 95%
        # of the stellar mass
        ###################################################################################

        if !isempty(data_dict[:stars]["MASS"])

            mass_height_90 = computeMassHeight(
                data_dict[:stars]["POS "],
                data_dict[:stars]["MASS"];
                percent=90.0,
            )

            mass_height_95 = computeMassHeight(
                data_dict[:stars]["POS "],
                data_dict[:stars]["MASS"];
                percent=95.0,
            )

            println(
                file,
                "\tTotal height of a cylinder, of infinite radius, containing X% of the stellar \
                mass (filtered box):\n"
            )
            println(file, "\t\t$(round(ustrip(u"kpc", mass_height_90), sigdigits=4)) $(u"kpc") (90%)")
            println(file, "\t\t$(round(ustrip(u"kpc", mass_height_95), sigdigits=4)) $(u"kpc") (95%)\n")

        end

        ############################################################################################
        # Print the properties of the target halo and subhalo
        ############################################################################################

        if !ismissing(groupcat_path) && isSubfindActive(groupcat_path)

            # Check that the requested halo index is within bounds
            n_groups_total = readGroupCatHeader(groupcat_path).n_groups_total
            (
                0 < halo_idx <= n_groups_total ||
                throw(ArgumentError("snapshotReport: There is only $(n_groups_total) FoF \
                groups in $(simulation_path), so `halo_idx` = $(halo_idx) is outside of bounds"))
            )

            request = Dict(
                :subhalo => [
                    "S_Mass",
                    "S_MassType",
                    "S_LenType",
                    "S_CM",
                    "S_Pos",
                    "S_Vel",
                    "S_HalfmassRad",

                ],
                :group => [
                    "G_Mass",
                    "G_MassType",
                    "G_M_Crit200",
                    "G_LenType",
                    "G_Nsubs",
                    "G_CM",
                    "G_Pos",
                    "G_Vel",
                    "G_R_Crit200",
                ],
            )

            # Read the necessary data
            gc_data = readGroupCatalog(groupcat_path, snapshot_path, request)

            # Check that the requested subhalo index is within bounds
            g_n_subs = gc_data[:group]["G_Nsubs"]
            n_subfinds = g_n_subs[halo_idx]
            (
                subhalo_rel_idx <= n_subfinds ||
                throw(ArgumentError("snapshotReport: There is only $(n_subfinds) subhalos \
                for the FoF group $(halo_idx) in $(simulation_path), so `subhalo_rel_idx` \
                = $(subhalo_rel_idx) is outside of bounds"))
            )

            # Compute the number of subhalos and particles up to the last halo before `halo_idx`
            if isone(halo_idx)
                n_subs_floor = 0
            else
                n_subs_floor = sum(g_n_subs[1:(halo_idx - 1)]; init=0)
            end

            # Compute the subhalo absolute index
            subhalo_abs_idx = n_subs_floor + subhalo_rel_idx

            # Load the necessary data
            s_mass          = gc_data[:subhalo]["S_Mass"][subhalo_abs_idx]
            s_mass_type     = gc_data[:subhalo]["S_MassType"][:, subhalo_abs_idx]
            s_len_type      = gc_data[:subhalo]["S_LenType"][:, subhalo_abs_idx]
            s_cm            = gc_data[:subhalo]["S_CM"][:, subhalo_abs_idx]
            s_pos           = gc_data[:subhalo]["S_Pos"][:, subhalo_abs_idx]
            s_vel           = gc_data[:subhalo]["S_Vel"][:, subhalo_abs_idx]
            s_half_mass_rad = gc_data[:subhalo]["S_HalfmassRad"][subhalo_abs_idx]
            g_mass          = gc_data[:group]["G_Mass"][halo_idx]
            g_mass_type     = gc_data[:group]["G_MassType"][:, halo_idx]
            g_m_crit_200    = gc_data[:group]["G_M_Crit200"][halo_idx]
            g_len_type      = gc_data[:group]["G_LenType"][:, halo_idx]
            g_n_subs        = gc_data[:group]["G_Nsubs"][halo_idx]
            g_cm            = gc_data[:group]["G_CM"][:, halo_idx]
            g_pos           = gc_data[:group]["G_Pos"][:, halo_idx]
            g_vel           = gc_data[:group]["G_Vel"][:, halo_idx]
            g_r_crit_200    = gc_data[:group]["G_R_Crit200"][halo_idx]

            #####################################
            # Print the fraction of insitu stars
            #####################################

            if snapshot_length >= 2

                insitu_idx = filterByBirthPlace(
                    data_dict,
                    :exsitu;
                    halo_idx,
                    subhalo_rel_idx,
                )[:stars]

                iMs = sum(data_dict[:stars]["MASS"][insitu_idx]; init=0.0u"Msun")
                tMs = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun")
                insitu_fraction = round(uconvert.(Unitful.NoUnits, (iMs / tMs) * 100); sigdigits=2)

                println(file, "#"^100)
                println(file, "\nFraction of insitu stars (filtered box): $(insitu_fraction)%\n")

            end

            if :gas in component_list

                ############################################################
                # Indices of cells and particles within the disc radius and
                # between the disc radius and the virial radius
                ############################################################

                disc_idxs = filterWithinSphere(data_dict, (0.0u"kpc", DISK_R), :zero)
                halo_idxs = filterWithinSphere(data_dict, (DISK_R, g_r_crit_200), :zero)

                stellar_masses     = computeMass(data_dict, :stars)
                gas_masses         = computeMass(data_dict, :gas)
                ionized_masses     = computeMass(data_dict, :ionized)
                atomic_masses      = computeMass(data_dict, :atomic)
                molecular_masses   = computeMass(data_dict, :molecular)
                molecular_p_masses = computeMass(data_dict, :br_molecular)
                neutral_masses     = computeMass(data_dict, :neutral)

                stellar_mass_inside  = stellar_masses[disc_idxs[:stars]]
                stellar_mass_outside = stellar_masses[halo_idxs[:stars]]

                gas_mass_inside  = gas_masses[disc_idxs[:gas]]
                gas_mass_outside = gas_masses[halo_idxs[:gas]]

                if !isempty(ionized_masses)
                    ionized_mass_inside  = ionized_masses[disc_idxs[:gas]]
                    ionized_mass_outside = ionized_masses[halo_idxs[:gas]]
                end

                if !isempty(atomic_masses)
                    atomic_mass_inside  = atomic_masses[disc_idxs[:gas]]
                    atomic_mass_outside = atomic_masses[halo_idxs[:gas]]
                end

                if !isempty(molecular_masses)
                    molecular_mass_inside  = molecular_masses[disc_idxs[:gas]]
                    molecular_mass_outside = molecular_masses[halo_idxs[:gas]]
                end

                if !isempty(molecular_p_masses)
                    molecular_p_mass_inside  = molecular_p_masses[disc_idxs[:gas]]
                    molecular_p_mass_outside = molecular_p_masses[halo_idxs[:gas]]
                end

                if !isempty(neutral_masses)
                    neutral_mass_inside  = neutral_masses[disc_idxs[:gas]]
                    neutral_mass_outside = neutral_masses[halo_idxs[:gas]]
                end

                println(file, "Characteristic radii (filtered box):\n")

                #######################################################
                # Print the radius containing 90% and 95% of the mass,
                # withing de disc (r < `DISK_R`)
                #######################################################

                ########
                # Stars
                ########

                mass_radius_90 = computeMassRadius(
                    data_dict[:stars]["POS "][:, disc_idxs[:stars]],
                    stellar_mass_inside;
                    percent=90.0,
                )

                mass_radius_95 = computeMassRadius(
                    data_dict[:stars]["POS "][:, disc_idxs[:stars]],
                    stellar_mass_inside;
                    percent=95.0,
                )

                println(file, "\tRadius containing X% of the stellar mass (r < $(DISK_R)):\n")
                println(
                    file,
                    "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                )
                println(
                    file,
                    "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                )

                ############
                # Total gas
                ############

                mass_radius_90 = computeMassRadius(
                    data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                    gas_mass_inside;
                    percent=90.0,
                )

                mass_radius_95 = computeMassRadius(
                    data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                    gas_mass_inside;
                    percent=95.0,
                )

                println(
                    file,
                    "\tRadius containing X% of the total gas mass (r < $(DISK_R)):\n",
                )
                println(
                    file,
                    "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                )
                println(
                    file,
                    "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                )

                ##############
                # Ionized gas
                ##############

                if !isempty(ionized_masses)

                    mass_radius_90 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        ionized_mass_inside;
                        percent=90.0,
                    )

                    mass_radius_95 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        ionized_mass_inside;
                        percent=95.0,
                    )

                    println(
                        file,
                        "\tRadius containing X% of the ionized gas mass (r < $(DISK_R)):\n",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                    )

                end

                #############
                # Atomic gas
                #############

                if !isempty(atomic_masses)

                    mass_radius_90 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        atomic_mass_inside;
                        percent=90.0,
                    )

                    mass_radius_95 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        atomic_mass_inside;
                        percent=95.0,
                    )

                    println(
                        file,
                        "\tRadius containing X% of the atomic gas mass (r < $(DISK_R)):\n",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                    )

                end

                ################
                # Molecular gas
                ################

                if !isempty(molecular_masses)

                    mass_radius_90 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        molecular_mass_inside;
                        percent=90.0,
                    )

                    mass_radius_95 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        molecular_mass_inside;
                        percent=95.0,
                    )

                    println(
                        file,
                        "\tRadius containing X% of the molecular gas mass (r < $(DISK_R)):\n",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                    )

                end

                #####################
                # Molecular gas (BR)
                #####################

                if !isempty(molecular_p_masses)

                    mass_radius_90 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        molecular_p_mass_inside;
                        percent=90.0,
                    )

                    mass_radius_95 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        molecular_p_mass_inside;
                        percent=95.0,
                    )

                    println(
                        file,
                        "\tRadius containing X% of the molecular gas mass (BR) (r < $(DISK_R)):\n",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                    )

                end

                ##############
                # Neutral gas
                ##############

                if !isempty(neutral_masses)

                    mass_radius_90 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        neutral_mass_inside;
                        percent=90.0,
                    )

                    mass_radius_95 = computeMassRadius(
                        data_dict[:gas]["POS "][:, disc_idxs[:gas]],
                        neutral_mass_inside;
                        percent=95.0,
                    )

                    println(
                        file,
                        "\tRadius containing X% of the neutral gas mass (r < $(DISK_R)):\n",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_90), sigdigits=4)) $(u"kpc") (90%)",
                    )
                    println(
                        file,
                        "\t\t$(round(ustrip(u"kpc", mass_radius_95), sigdigits=4)) $(u"kpc") (95%)\n",
                    )

                end

                #####################################################
                # Print the masses withing de disc (r < DISK_R) and
                # outside the disc (DISK_R < r < R200)
                #####################################################

                println(file, "#"^100)
                println(file, "\nCharacteristic fractions and masses (filtered box):\n")

                println(file, "\t", "#"^20)
                println(file, "\tR200: $(round(typeof(1.0u"kpc"), g_r_crit_200, sigdigits=4))")
                println(file, "\t", "#"^20, "\n")

                ########
                # Stars
                ########

                total_stellar_mass_inside  = sum(stellar_mass_inside; init=0.0u"Msun")
                total_stellar_mass_outside = sum(stellar_mass_outside; init=0.0u"Msun")
                total_stellar_mass         = total_stellar_mass_inside + total_stellar_mass_outside

                s_inside_percent  = (total_stellar_mass_inside / total_stellar_mass) * 100.0
                s_outside_percent = (total_stellar_mass_outside / total_stellar_mass) * 100.0

                println(file, "\t", "#"^40)
                println(file, "\tStars:")
                println(file, "\t", "#"^40, "\n")

                println(file, "\tStellar mass inside the disc (r < $(DISK_R)):\n")
                println(
                    file,
                    "\t\t$(round(typeof(1.0u"Msun"), total_stellar_mass_inside, sigdigits=3)) \
                    ($(round(s_inside_percent, sigdigits=3))% of the total stellar mass)\n",
                )

                println(file, "\tStellar mass outside the disc ($(DISK_R) < r < R200):\n")
                println(
                    file,
                    "\t\t$(round(typeof(1.0u"Msun"), total_stellar_mass_outside, sigdigits=3)) \
                    ($(round(s_outside_percent, sigdigits=3))% of the total stellar mass)\n",
                )

                ############
                # Total gas
                ############

                total_gas_mass_inside  = sum(gas_mass_inside; init=0.0u"Msun")
                total_gas_mass_outside = sum(gas_mass_outside; init=0.0u"Msun")
                total_gas_mass         = total_gas_mass_inside + total_gas_mass_outside

                g_inside_percent  = (total_gas_mass_inside / total_gas_mass) * 100.0
                g_outside_percent = (total_gas_mass_outside / total_gas_mass) * 100.0

                println(file, "\t", "#"^40)
                println(file, "\tTotal gas:")
                println(file, "\t", "#"^40, "\n")

                println(file, "\tGas mass inside the disc (r < $(DISK_R)):\n")
                println(
                    file,
                    "\t\t$(round(typeof(1.0u"Msun"), total_gas_mass_inside, sigdigits=3)) \
                    ($(round(g_inside_percent, sigdigits=3))% of the gas mass)\n",
                )

                println(file, "\tGas mass outside the disc ($(DISK_R) < r < R200):\n")
                println(
                    file,
                    "\t\t$(round(typeof(1.0u"Msun"), total_gas_mass_outside, sigdigits=3)) \
                    ($(round(g_outside_percent, sigdigits=3))% of the gas mass)\n",
                )

                ##############
                # Ionized gas
                ##############

                if !isempty(ionized_masses)

                    total_ion_mass_inside  = sum(ionized_mass_inside; init=0.0u"Msun")
                    total_ion_mass_outside = sum(ionized_mass_outside; init=0.0u"Msun")

                    i_inside_percent  = (total_ion_mass_inside  / total_gas_mass_inside) * 100.0
                    i_outside_percent = (total_ion_mass_outside / total_gas_mass_outside) * 100.0

                    println(file, "\t", "#"^40)
                    println(file, "\tIonized gas:")
                    println(file, "\t", "#"^40, "\n")

                    println(file, "\tIonized mass inside the disc (r < $(DISK_R)):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_ion_mass_inside, sigdigits=3)) \
                        ($(round(i_inside_percent, sigdigits=3))% of the gas mass)\n",
                    )

                    println(file, "\tIonized mass outside the disc ($(DISK_R) < r < R200):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_ion_mass_outside, sigdigits=3)) \
                        ($(round(i_outside_percent, sigdigits=3))% of the gas mass)\n",
                    )

                end

                ##############
                # Neutral gas
                ##############

                if !isempty(neutral_masses)

                    total_neu_mass_inside  = sum(neutral_mass_inside; init=0.0u"Msun")
                    total_neu_mass_outside = sum(neutral_mass_outside; init=0.0u"Msun")

                    m_inside_percent  = (total_neu_mass_inside  / total_gas_mass_inside) * 100.0
                    m_outside_percent = (total_neu_mass_outside / total_gas_mass_outside) * 100.0

                    println(file, "\t", "#"^40)
                    println(file, "\tNeutral gas:")
                    println(file, "\t", "#"^40, "\n")

                    println(file, "\tNeutral mass inside the disc (r < $(DISK_R)):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_neu_mass_inside, sigdigits=3)) \
                        ($(round(m_inside_percent, sigdigits=3))% of the gas mass)\n",
                    )

                    println(file, "\tNeutral mass outside the disc ($(DISK_R) < r < R200):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_neu_mass_outside, sigdigits=3)) \
                        ($(round(m_outside_percent, sigdigits=3))% of the gas mass)\n",
                    )

                end

                #############
                # Atomic gas
                #############

                if !isempty(atomic_masses)

                    total_ato_mass_inside  = sum(atomic_mass_inside; init=0.0u"Msun")
                    total_ato_mass_outside = sum(atomic_mass_outside; init=0.0u"Msun")

                    a_inside_percent  = (total_ato_mass_inside  / total_gas_mass_inside) * 100.0
                    a_outside_percent = (total_ato_mass_outside / total_gas_mass_outside) * 100.0

                    println(file, "\t", "#"^40)
                    println(file, "\tAtomic gas:")
                    println(file, "\t", "#"^40, "\n")

                    println(file, "\tAtomic mass inside the disc (r < $(DISK_R)):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_ato_mass_inside, sigdigits=3)) \
                        ($(round(a_inside_percent, sigdigits=3))% of the gas mass)\n",
                    )

                    println(file, "\tAtomic mass outside the disc ($(DISK_R) < r < R200):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_ato_mass_outside, sigdigits=3)) \
                        ($(round(a_outside_percent, sigdigits=3))% of the gas mass)\n",
                    )

                end

                ################
                # Molecular gas
                ################

                if !isempty(molecular_masses)

                    total_mol_mass_inside  = sum(molecular_mass_inside; init=0.0u"Msun")
                    total_mol_mass_outside = sum(molecular_mass_outside; init=0.0u"Msun")

                    m_inside_percent  = (total_mol_mass_inside  / total_gas_mass_inside) * 100.0
                    m_outside_percent = (total_mol_mass_outside / total_gas_mass_outside) * 100.0

                    println(file, "\t", "#"^40)
                    println(file, "\tMolecular gas:")
                    println(file, "\t", "#"^40, "\n")

                    println(file, "\tMolecular mass inside the disc (r < $(DISK_R)):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_mol_mass_inside, sigdigits=3)) \
                        ($(round(m_inside_percent, sigdigits=3))% of the gas mass)\n",
                    )

                    println(file, "\tMolecular mass outside the disc ($(DISK_R) < r < R200):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_mol_mass_outside, sigdigits=3)) \
                        ($(round(m_outside_percent, sigdigits=3))% of the gas mass)\n",
                    )

                end

                #####################
                # Molecular gas (BR)
                #####################

                if !isempty(molecular_p_masses)

                    total_mol_p_mass_inside  = sum(molecular_p_mass_inside; init=0.0u"Msun")
                    total_mol_p_mass_outside = sum(molecular_p_mass_outside; init=0.0u"Msun")

                    m_p_inside_percent  = (total_mol_p_mass_inside  / total_gas_mass_inside) * 100.0
                    m_p_outside_percent = (total_mol_p_mass_outside / total_gas_mass_outside) * 100.0

                    println(file, "\t", "#"^40)
                    println(file, "\tMolecular gas (BR recipe):")
                    println(file, "\t", "#"^40, "\n")

                    println(file, "\tMolecular mass (BR recipe) inside the disc (r < $(DISK_R)):\n")
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_mol_p_mass_inside, sigdigits=3)) \
                        ($(round(m_p_inside_percent, sigdigits=3))% of the gas mass)\n",
                    )

                    println(
                        file,
                        "\tMolecular mass (BR recipe) outside the disc ($(DISK_R) < r < R200):\n",
                    )
                    println(
                        file,
                        "\t\t$(round(typeof(1.0u"Msun"), total_mol_p_mass_outside, sigdigits=3)) \
                        ($(round(m_p_outside_percent, sigdigits=3))% of the gas mass)\n",
                    )

                end

            end

            ########################################################################################
            # Halo and subhalo global properties
            ########################################################################################

            println(file, "#"^100)
            println(file, "\nHalo and subhalo global properties:")

            ############################
            # Print the halo properties
            ############################

            println(file, "\n", "#"^71)
            println(file, "NOTE: Stellar particle counts include wind particles from here on out!")
            println(file, "#"^71)

            # Compute halo number
            halo_n = lpad(halo_idx - 1, 3, "0")

            println(file, "\nHalo $(halo_n) properties:\n")

            ########################################################################################

            println(file, "\tCell/particle number (in halo $(halo_n)):\n")
            for (i, len) in pairs(g_len_type)

                component = PARTICLE_NAMES[INDEX_PARTICLE[i - 1]]
                println(file, "\t\t$(component):$(" "^(22 - length(component))) $(len)")

            end

            ########################################################################################

            println(file, "\n\tNumber of subhalos (in halo $(halo_n)):\n\n\t\t$(g_n_subs)\n")

            ########################################################################################

            println(file, "\tMasses (in halo $(halo_n)):\n")
            for (i, mass) in pairs(g_mass_type)

                symbol_name = INDEX_PARTICLE[i - 1]
                component_label = PARTICLE_NAMES[symbol_name]

                if symbol_name == :stars
                    component_label = "Stellar/Wind particles"
                elseif symbol_name == :tracer
                    mass = TRACER_MASS * internalUnits("MASS", snapshot_path)
                end

                println(
                    file,
                    "\t\t$(component_label):$(" "^(22 - length(component_label))) \
                    $(round(typeof(1.0u"Msun"), mass, sigdigits=3))",
                )

            end

            println(
                file,
                "\n\t\tTotal mass:             $(round.(typeof(1.0u"Msun"), g_mass, sigdigits=3))",
            )

            ########################################################################################

            println(
                file,
                "\n\tCenter of mass (in halo $(halo_n)): \
                \n\n\t\t$(round.(ustrip.(u"Mpc", g_cm), sigdigits=6)) $(u"Mpc")\n",
            )

            println(
                file,
                "\tPosition of the particle with the minimum gravitational potential energy \
                (in halo $(halo_n)): \n\n\t\t$(round.(ustrip.(u"Mpc", g_pos), sigdigits=6)) $(u"Mpc")\n",
            )

            separation = sqrt(sum((g_cm - g_pos) .^ 2))
            println(
                file,
                "\tSeparation between the minimum potencial and the global CM (in halo $(halo_n)): \
                \n\n\t\t$(round(typeof(1.0u"kpc"), separation, sigdigits=6))\n",
            )

            ########################################################################################

            vel_cm = round.(Float64.(ustrip.(u"km*s^-1", g_vel)), sigdigits=6)
            println(
                file,
                "\tVelocity of the center of mass (in halo $(halo_n)):\n\n\t\t$(vel_cm) \
                $(u"km*s^-1")\n",
            )

            ########################################################################################

            println(
                file,
                "\tTotal mass enclosed in a sphere with a mean density 200 times the critical \
                density (in halo $(halo_n)): \
                \n\n\t\t$(round(typeof(1.0u"Msun"), g_m_crit_200, sigdigits=3))\n",
            )

            ########################################################################################

            println(
                file,
                "\tRadius of a sphere with a mean density 200 times the critical density \
                (in halo $(halo_n)): \n\n\t\t$(round(typeof(1.0u"kpc"), g_r_crit_200, sigdigits=4))\n",
            )

            ###############################
            # Print the subhalo properties
            ###############################

            # Compute subhalo number
            subhalo_n = lpad(subhalo_rel_idx - 1, 3, "0")

            println(file, "#"^100)
            println(
                file,
                "\nSubhalo $(subhalo_n) (of halo $(halo_n)) properties:\n",
            )

            ########################################################################################

            println(file, "\tCell/particle number (in subhalo $(subhalo_n)):\n")
            for (i, len) in pairs(s_len_type)

                component = PARTICLE_NAMES[INDEX_PARTICLE[i - 1]]
                println(file, "\t\t$(component):$(" "^(22 - length(component))) $(len)")

            end

            ########################################################################################

            println(file, "\n\tMasses (in subhalo $(subhalo_n)):\n")
            for (i, mass) in pairs(s_mass_type)

                symbol_name = INDEX_PARTICLE[i - 1]
                component_label = PARTICLE_NAMES[symbol_name]

                if symbol_name == :stars
                    component_label = "Stellar/Wind particles"
                elseif symbol_name == :tracer
                    mass = TRACER_MASS * internalUnits("MASS", snapshot_path)
                end

                println(
                    file,
                    "\t\t$(component_label):$(" "^(22 - length(component_label))) \
                    $(round(typeof(1.0u"Msun"), mass, sigdigits=3))",
                )

            end

            println(
                file,
                "\n\t\tTotal mass:             $(round.(typeof(1.0u"Msun"), s_mass, sigdigits=3))",
            )

            ########################################################################################

            println(
                file,
                "\n\tCenter of mass (in subhalo $(subhalo_n)): \
                \n\n\t\t$(round.(ustrip.(u"Mpc", s_cm), sigdigits=6)) $(u"Mpc")\n",
            )

            println(
                file,
                "\tPosition of the particle with the minimum gravitational potential energy \
                (in subhalo $(subhalo_n)): \n\n\t\t$(round.(ustrip.(u"Mpc", s_pos), sigdigits=6)) $(u"Mpc")\n",
            )

            separation = sqrt(sum((s_cm - s_pos) .^ 2))
            println(
                file,
                "\tSeparation between the minimum potencial and the global CM (in subhalo $(subhalo_n)): \
                \n\n\t\t$(round(typeof(1.0u"kpc"), separation, sigdigits=6))\n",
            )

            ########################################################################################

            vel_cm = round.(Float64.(ustrip.(u"km*s^-1", s_vel)), sigdigits=6)
            println(
                file,
                "\tVelocity of the center of mass (in subhalo $(subhalo_n)): \
                \n\n\t\t$(vel_cm) $(u"km*s^-1")\n",
            )

            ########################################################################################

            println(
                file,
                "\tRadius containing half of the total mass (in subhalo $(subhalo_n)): \
                \n\n\t\t$(round(typeof(1.0u"kpc"), s_half_mass_rad, sigdigits=4))",
            )

        end

        close(file)

    end

    return nothing

end

"""
    simulationReport(
        simulation_paths::Vector{String};
        <keyword arguments>
    )::Nothing

Write a text file with information about a given simulation

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. One text file will be printed for each simulation.
  - `output_path::String="./"`: Path to the output folder.
"""
function simulationReport(
    simulation_paths::Vector{String};
    output_path::String="./",
)::Nothing

    for simulation_path in simulation_paths

        # Make a dataframe with the following columns:
        #   - 1. DataFrame index
        #   - 2. Number in the file name
        #   - 3. Scale factor
        #   - 4. Redshift
        #   - 5. Physical time
        #   - 6. Lookback time
        #   - 7. Snapshot path
        #   - 8. Group catalog path
        simulation_table = makeSimulationTable(simulation_path)

        # Compute the number of snapshots in the folder
        snapshot_n = count(!ismissing, simulation_table[!, :snapshot_paths])

        # Check that there is at least one snapshot
        (
            !iszero(snapshot_n) ||
            throw(ArgumentError("simulationReport: There are no snapshots in $(simulation_path)"))
        )

        # Compute the number of group catalog files in the folder
        groupcat_n = count(!ismissing, simulation_table[!, :groupcat_paths])

        # Check if the simulation is cosmological
        cosmological = isCosmological(first(skipmissing(simulation_table[!, :snapshot_paths])))

        ############################################################################################
        # Print the report header
        ############################################################################################

        # Create the output file
        file = open(
            joinpath(mkpath(output_path), "report_for_$(basename(simulation_path)).txt"),
            "w",
        )

        println(file, "#"^100)
        println(file, "\nSimulation name:  $(basename(simulation_path))")

        if cosmological
            println(file, "Cosmological:     Yes")
        else
            println(file, "Cosmological:     No")
        end

        if !PHYSICAL_UNITS && cosmological
            println(file, "Report units:     Comoving\n")
        else
            println(file, "Report units:     Physical\n")
        end

        println(file, "#"^100)
        println(file, "\nNumber of snapshots:       $(snapshot_n)")
        println(file, "Number of group catalogs:  $(groupcat_n)\n")

        # Print the simulation time ranges
        if snapshot_n > 1

            min_pt, max_pt = round.(
                ustrip.(u"Gyr", extrema(simulation_table[!, :physical_times])),
                digits=2,
            )

            println(file, "Physical time range:       $(min_pt) - $(max_pt) Gyr")

            if cosmological

                # For cosmological simulations print the scale factor and the redshift
                min_a, max_a = round.(extrema(simulation_table[!, :scale_factors]), digits=3)
                min_z, max_z = round.(extrema(simulation_table[!, :redshifts]), digits=3)

                println(file, "Scale factor range:        $(min_a) - $(max_a)")
                println(file, "Redshift range:            $(max_z) - $(min_z)")

            end

            println(file,)

        else

            pt = round(ustrip(u"Gyr", simulation_table[1, :physical_times]), digits=2)

            println(file, "Physical time:             $(pt) Gyr")

            if cosmological

                # For cosmological simulations print the scale factor and the redshift
                a = round(simulation_table[1, :scale_factors], digits=3)
                z = round(simulation_table[1, :redshifts], digits=3)

                println(file, "Scale factor:              $(a)")
                println(file, "Redshift:                  $(z)")

            end

            println(file,)

        end

        # Set flags to print only the first instance of each condition
        first_star_flag            = false
        first_subhalo_flag         = false
        first_star_in_subhalo_flag = false

        for snapshot_row in eachrow(simulation_table)

            snapshot_path = snapshot_row[:snapshot_paths]

            # Skip this row if there is no snapshot
            !ismissing(snapshot_path) || continue

            # Read the snapshot header
            snapshot_header = readSnapHeader(snapshot_path)

            # Read the group catalog path
            groupcat_path = snapshot_row[:groupcat_paths]

            # Read the different time ticks
            physical_time = round(ustrip(u"Gyr", snapshot_row[:physical_times]), digits=2)
            if cosmological
                # For cosmological simulations read the scale factor and the redshift
                scale_factor = round(snapshot_row[:scale_factors], digits=3)
                redshift = round(snapshot_row[:redshifts], digits=3)
            end

            # Read how many stars there are in this snapshot
            star_number = countStars(snapshot_path)

            # Check if there is subfind information in the group catalog file
            subfind_active = !ismissing(groupcat_path) && isSubfindActive(groupcat_path)

            if subfind_active

                # Read the group catalog header
                groupcat_header = readGroupCatHeader(groupcat_path)

                # Read the number of halos
                n_groups_total = readGroupCatHeader(groupcat_path).n_groups_total

                # Make the subfind request
                request = Dict(:subhalo => ["S_LenType", "S_CM", "S_Pos"])

                # Read the necessary data
                gc_data = readGroupCatalog(groupcat_path, snapshot_path, request)

                # Load the necessary data
                s_cm       = gc_data[:subhalo]["S_CM"][:, 1]
                s_pos      = gc_data[:subhalo]["S_Pos"][:, 1]
                s_len_type = gc_data[:subhalo]["S_LenType"][:, 1]

                # Select the filter function and request dictionary
                filter_function, _, _, request = selectFilter(:subhalo, Dict(:stars => ["MASS"]))

                # Create a metadata dictionary
                metadata = Dict(
                    :snap_data => Snapshot(
                        snapshot_path,
                        1,
                        1,
                        0.0u"yr",
                        0.0u"yr",
                        0.0,
                        0.0,
                        snapshot_header,
                    ),
                    :gc_data => GroupCatalog(groupcat_path, groupcat_header),
                )

                # Create the data dictionary
                data_dict = merge(
                    metadata,
                    readSnapshot(snapshot_path, request),
                    readGroupCatalog(groupcat_path, snapshot_path, request),
                )

                filterData!(data_dict; filter_function)

                # Compute the number of stars in the main subhalo
                stellar_n_subhalo = length(data_dict[:stars]["MASS"])

            end

            ########################################################################################
            # First stars
            ########################################################################################

            if star_number > 0 && !first_star_flag

                println(file, "#"^100)
                println(file, "\nFirst snapshot with star formation:")

                println(file, "\n\tSnapshot:         $(basename(snapshot_path))")
                println(file, "\tPhysical time:    $(physical_time) Gyr")

                if cosmological
                    # For cosmological simulations print the scale factor and the redshift
                    println(file, "\tScale factor:     $(scale_factor)")
                    println(file, "\tRedshift:         $(redshift)")
                end

                println(file, "\tNumber of stars:  $(star_number)")

                if subfind_active

                    println(file, "\tNumber of halos:  $(n_groups_total)")

                    println(file, "\n\tMain subhalo properties:")

                    println(file, "\n\t\tNumber of stellar particles:  $(stellar_n_subhalo)\n")

                    println(
                        file,
                        "\t\tCenter of mass:\n\n\t\t\t\
                        $(round.(ustrip.(u"Mpc", s_cm), sigdigits=6)) $(u"Mpc")\n",
                    )

                    println(
                        file,
                        "\t\tPosition of the particle with the minimum gravitational potential \
                        energy: \n\n\t\t\t\
                        $(round.(ustrip.(u"Mpc", s_pos), sigdigits=6)) $(u"Mpc")\n",
                    )

                    separation = sqrt(sum((s_cm - s_pos) .^ 2))
                    println(
                        file,
                        "\t\tSeparation between the minimum potencial and the global CM: \
                        \n\n\t\t\t$(round(typeof(1.0u"kpc"), separation, sigdigits=6))\n",
                    )

                else
                    println(file, "\n" * "#"^51)
                    println(file, "There is no subfind information for this snapshot!")
                    println(file, "#"^51 * "\n")
                end

                first_star_flag = true

            end

            ########################################################################################
            # First subhalos
            ########################################################################################

            if subfind_active && !first_subhalo_flag

                # Read the number of halos
                n_groups_total = readGroupCatHeader(groupcat_path).n_groups_total

                println(file, "#"^100)
                println(file, "\nFirst snapshot with subfind information:")

                println(file, "\n\tSnapshot:         $(basename(snapshot_path))")

                println(file, "\tPhysical time:    $(physical_time) Gyr")

                if cosmological
                    # For cosmological simulations print the scale factor and the redshift
                    println(file, "\tScale factor:     $(scale_factor)")
                    println(file, "\tRedshift:         $(redshift)")
                end

                println(file, "\tNumber of halos:  $(n_groups_total)")

                println(file, "\n\tMain subhalo properties:")

                println(
                    file,
                    "\n\t\tCenter of mass:\n\n\t\t\t$(round.(ustrip.(u"Mpc", s_cm), sigdigits=6)) \
                    $(u"Mpc")\n",
                )

                println(
                    file,
                    "\t\tPosition of the particle with the minimum gravitational potential energy: \
                    \n\n\t\t\t$(round.(ustrip.(u"Mpc", s_pos), sigdigits=6)) $(u"Mpc")\n",
                )

                separation = sqrt(sum((s_cm - s_pos) .^ 2))
                println(
                    file,
                    "\t\tSeparation between the minimum potencial and the global CM: \
                    \n\n\t\t\t$(round(typeof(1.0u"kpc"), separation, sigdigits=6))\n",
                )

                println(file, "\t\t", "#"^54)
                println(file, "\t\tNOTE: Stellar particle counts include wind particles!")
                println(file, "\t\t", "#"^54)

                println(file, "\n\t\tCell/particle number:\n")

                for (i, len) in pairs(s_len_type)

                    component = PARTICLE_NAMES[INDEX_PARTICLE[i - 1]]
                    println(file, "\t\t\t$(component):$(" "^(22 - length(component))) $(len)")

                end

                println(file)

                first_subhalo_flag = true

            end

            ########################################################################################
            # First stars in the main subhalo
            ########################################################################################

            if subfind_active && !first_star_in_subhalo_flag && stellar_n_subhalo > 0

                println(file, "#"^100)
                println(file, "\nFirst snapshot with star formation in the main subhalo:")

                println(file, "\n\tSnapshot:         $(basename(snapshot_path))")
                println(file, "\tPhysical time:    $(physical_time) Gyr")

                if cosmological
                    # For cosmological simulations print the scale factor and the redshift
                    println(file, "\tScale factor:     $(scale_factor)")
                    println(file, "\tRedshift:         $(redshift)")
                end

                println(file, "\tNumber of stars:  $(stellar_n_subhalo)\n")

                first_star_in_subhalo_flag = true

            end

            # End the loop when all the conditions have been met
            if first_star_flag && first_subhalo_flag && first_star_in_subhalo_flag
                break
            end

        end

        close(file)

    end

    return nothing

end

"""
    sfrTXT(
        simulation_paths::Vector{String},
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a time series of the data in the `sfr.txt` file.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor.
      + `:redshift`      -> Redshift.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass` -> Stellar mass.
      + `:sfr`          -> The star formation rate.
  - `smooth::Int=0`: The result will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `output_path::String="./"`: Path to the output folder.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function sfrTXT(
    simulation_paths::Vector{String},
    x_quantity::Symbol,
    y_quantity::Symbol;
    smooth::Int=0,
    output_path::String="./",
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    plotTimeSeries(
        simulation_paths,
        [lines!];
        pf_kwargs=[(;)],
        # `plotTimeSeries` configuration
        output_path,
        filename="$(y_quantity)_vs_$(x_quantity)",
        output_format=".png",
        show_progress=true,
        # Data manipulation options
        slice=(:),
        da_functions=[daSFRtxt],
        da_args=[(x_quantity, y_quantity)],
        da_kwargs=[(; smooth)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        # Plotting options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
    )

    return nothing

end

"""
    cpuTXT(
        simulation_paths::Vector{String},
        target::String,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a time series of the data in the `cpu.txt` file.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `target::String`: Target process.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:time_step`              -> Time step.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:clock_time_s`           -> Clock time duration of the time step in seconds.
      + `:clock_time_percent`     -> Clock time duration of the time step as a percentage.
      + `:tot_clock_time_s`       -> Total clock time in seconds.
      + `:tot_clock_time_percent` -> Total clock time as a percentage.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:time_step`              -> Time step.
      + `:physical_time`          -> Physical time since the Big Bang.
      + `:clock_time_s`           -> Clock time duration of the time step in seconds.
      + `:clock_time_percent`     -> Clock time duration of the time step as a percentage.
      + `:tot_clock_time_s`       -> Total clock time in seconds.
      + `:tot_clock_time_percent` -> Total clock time as a percentage.
  - `smooth::Int=0`: The result will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `yscale::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `x_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the x coordinates fit within `x_trim`.
  - `y_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the y coordinates fit within `y_trim`. This option does not affect histograms.
  - `output_path::String="./"`: Path to the output folder.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function cpuTXT(
    simulation_paths::Vector{String},
    target::String,
    x_quantity::Symbol,
    y_quantity::Symbol;
    smooth::Int=0,
    yscale::Function=identity,
    x_trim::NTuple{2,<:Real}=(-Inf, Inf),
    y_trim::NTuple{2,<:Real}=(-Inf, Inf),
    output_path::String="./",
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    safe_str_target = replace(target, "/" => "-", "_" => "-")

    plotTimeSeries(
        simulation_paths,
        [lines!];
        pf_kwargs=[(;)],
        # `plotTimeSeries` configuration
        output_path,
        filename="$(y_quantity)_vs_$(x_quantity)_for_$(safe_str_target)",
        output_format=".png",
        show_progress=true,
        # Data manipulation options
        slice=(:),
        da_functions=[daCPUtxt],
        da_args=[(target, x_quantity, y_quantity)],
        da_kwargs=[(; smooth)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim,
        y_trim,
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=yscale,
        # Plotting options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title=L"\mathrm{Process: \,\, %$(safe_str_target)}",
    )

    return nothing

end

"""
    stellarBirthHalos(
        simulation_path::String,
        slice_n::Int;
        <keyword arguments>
    )::Nothing

Write, to a pair of CSV files, in which halo and subhalo every star in snapshot `slice_n` was born.

# Arguments

  - `simulation_paths::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects the target snapshot. Starts at 1 and is independent of the number in the file name. If every snapshot is present, the relation is `slice_n` = filename_number + 1.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
"""
function stellarBirthHalos(
    simulation_path::String,
    slice_n::Int;
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
)::Nothing

    # Select the filter function and request dictionary
    filter_function, _, _, request = selectFilter(filter_mode, Dict(:stars=>["ID  "]))

    # Read the relevant data of the snapshot
    data_dict = makeDataDict(
        simulation_path,
        slice_n,
        request,
    )

    # Filter the data
    filterData!(data_dict; filter_function)

    # Find the birth place of every star
    birth_halo, birth_subhalo = locateStellarBirthPlace(data_dict)

    # Write the results to CSV files
    CSV.write(
        joinpath(output_path, "stellar_birth_halos.gz"),
        Tables.table(birth_halo);
        newline=',',
        writeheader=false,
        compress=true,
    )
    CSV.write(
        joinpath(output_path, "stellar_birth_subhalos.gz"),
        Tables.table(birth_subhalo);
        newline=',',
        writeheader=false,
        compress=true,
    )

    return nothing

end

"""
    densityMap(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D histogram of the density.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored. If set to 0, an animation using every snapshots will be made.
  - `quantities::Vector{Symbol}=[:gas_mass]`: Quantities for which the density will be calculated. The options are:

      + `:stellar_mass`      -> Stellar mass.
      + `:gas_mass`          -> Gas mass.
      + `:hydrogen_mass`     -> Hydrogen mass.
      + `:dm_mass`           -> Dark matter mass.
      + `:bh_mass`           -> Black hole mass.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`      -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`  -> Stellar gas mass (according to out SF model).
  - `types::Vector{Symbol}=[:cells]`: List of component types for the density fields, each element can be either `:particles` or Voronoi `:cells`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Pixel (bin of the 2D histogram) side length.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection, averaging the value of neighboring pixels. It has to divide the size of `grid` exactly.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the values to be plotted.
  - `da_ff::Function=filterNothing`: Filter function for the data analysis function. It must be a function with the signature:

    `da_ff(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `da_ff` filter function.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function densityMap(
    simulation_paths::Vector{String},
    slice::IndexType;
    quantities::Vector{Symbol}=[:gas_mass],
    types::Vector{Symbol}=[:cells],
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    reduce_factor::Int=1,
    theme::Attributes=Theme(),
    title::Union{Symbol,<:AbstractString}="",
    annotation::AbstractString="",
    colorbar::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid
    grid = CubicGrid(box_size, resolution)

    pf_kwargs = isnothing(colorrange) ? [(;)] : [(; colorrange)]

    for (i, quantity) in pairs(quantities)

        filter_function, translation, rotation, request = selectFilter(
            filter_mode,
            mergeRequests(plotParams(quantity).request, ff_request),
        )

        for simulation_path in simulation_paths

            # Get the simulation name as a string
            sim_name = basename(simulation_path)

            for projection_plane in projection_planes

                # Construct the file name
                base_filename = "$(sim_name)_$(quantity)_$(projection_plane)_density_map"

                plotSnapshot(
                    [simulation_path],
                    request,
                    [heatmap!];
                    pf_kwargs,
                    # `plotSnapshot` configuration
                    output_path,
                    base_filename,
                    output_format=".png",
                    show_progress=true,
                    # Data manipulation options
                    slice=iszero(slice) ? (:) : slice,
                    filter_function,
                    da_functions=[daDensity2DProjection],
                    da_args=[(grid, quantity, ring(types, i))],
                    da_kwargs=[(; reduce_factor, projection_plane, filter_function=da_ff)],
                    post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                    pp_args=(annotation,),
                    pp_kwargs=(; color=:white),
                    transform_box=true,
                    translation,
                    rotation,
                    smooth=0,
                    x_unit=u"kpc",
                    y_unit=u"kpc",
                    x_exp_factor=0,
                    y_exp_factor=0,
                    x_trim=(-Inf, Inf),
                    y_trim=(-Inf, Inf),
                    x_edges=false,
                    y_edges=false,
                    x_func=identity,
                    y_func=identity,
                    # Axes options
                    xaxis_label="auto_label",
                    yaxis_label="auto_label",
                    xaxis_var_name=string(projection_plane)[1:1],
                    yaxis_var_name=string(projection_plane)[2:2],
                    xaxis_scale_func=identity,
                    yaxis_scale_func=identity,
                    # Plotting options
                    save_figures=!iszero(slice),
                    backup_results=iszero(slice),
                    theme=merge(
                        theme,
                        Theme(
                            size=colorbar ? (880, 760) : (880, 880),
                            figure_padding=(5, 20, 20, 10),
                            Axis=(limits=(-limit, limit, -limit, limit),),
                        ),
                    ),
                    sim_labels=nothing,
                    title,
                    colorbar,
                    # Animation options
                    animation=iszero(slice),
                    animation_filename="$(base_filename).mp4",
                    framerate=5,
                )

            end

        end

    end

    return nothing

end

"""
    gasSFRMap(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D map of the gas SFR.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored. If set to 0, an animation using every snapshots will be made.
  - `types::Symbol=:cells`: Gas type for the SFR fields. It can be either `:particles` or Voronoi `:cells`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Pixel (bin of the 2D histogram) side length.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the values to be plotted.
  - `da_ff::Function=filterNothing`: Filter function for the data analysis function. It must be a function with the signature:

    `da_ff(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `da_ff` filter function.
"""
function gasSFRMap(
    simulation_paths::Vector{String},
    slice::IndexType;
    type::Symbol=:cells,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    theme::Attributes=Theme(),
    title::Union{Symbol,<:AbstractString}="",
    annotation::AbstractString="",
    colorbar::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid
    grid = CubicGrid(box_size, resolution)

    pf_kwargs = isnothing(colorrange) ? [(;)] : [(; colorrange)]

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(
            plotParams(:gas_sfr).request,
            plotParams(:gas_mass_density).request,
            ff_request,
        ),
    )

    for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        for projection_plane in projection_planes

            # Construct the file name
            base_filename = "$(sim_name)_$(projection_plane)_gas_sfr_map"

            plotSnapshot(
                [simulation_path],
                request,
                [heatmap!];
                pf_kwargs,
                # `plotSnapshot` configuration
                output_path,
                base_filename,
                output_format=".png",
                show_progress=true,
                # Data manipulation options
                slice=iszero(slice) ? (:) : slice,
                filter_function,
                da_functions=[daGasSFR2DProjection],
                da_args=[(grid, type)],
                da_kwargs=[(; projection_plane, filter_function=da_ff)],
                post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                pp_args=(annotation,),
                pp_kwargs=(; color=:white),
                transform_box=true,
                translation,
                rotation,
                smooth=0,
                x_unit=u"kpc",
                y_unit=u"kpc",
                x_exp_factor=0,
                y_exp_factor=0,
                x_trim=(-Inf, Inf),
                y_trim=(-Inf, Inf),
                x_edges=false,
                y_edges=false,
                x_func=identity,
                y_func=identity,
                # Axes options
                xaxis_label="auto_label",
                yaxis_label="auto_label",
                xaxis_var_name=string(projection_plane)[1:1],
                yaxis_var_name=string(projection_plane)[2:2],
                xaxis_scale_func=identity,
                yaxis_scale_func=identity,
                # Plotting options
                save_figures=!iszero(slice),
                backup_results=iszero(slice),
                theme=merge(
                    theme,
                    Theme(
                        size=colorbar ? (880, 760) : (880, 880),
                        figure_padding=(5, 20, 20, 10),
                        Axis=(limits=(-limit, limit, -limit, limit),),
                    ),
                ),
                sim_labels=nothing,
                title,
                colorbar,
                # Animation options
                animation=iszero(slice),
                animation_filename="$(base_filename).mp4",
                framerate=5,
            )

        end

    end

    return nothing

end

"""
    densityMapVelField(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D histogram of the density, with the velocity field.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored. If set to 0, an animation using every snapshots will be made.
  - `quantities::Vector{Symbol}=[:gas_mass]`: Quantities for which the density will be calculated. The options are:

      + `:stellar_mass`      -> Stellar mass.
      + `:gas_mass`          -> Gas mass.
      + `:hydrogen_mass`     -> Hydrogen mass.
      + `:dm_mass`           -> Dark matter mass.
      + `:bh_mass`           -> Black hole mass.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`      -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`  -> Stellar gas mass (according to out SF model).
  - `types::Vector{Symbol}=[:cells]`: List of component types for the density fields, each element can be either `:particles` or Voronoi `:cells`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Pixel (bin of the 2D histogram) side length.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the values to be plotted.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function densityMapVelField(
    simulation_paths::Vector{String},
    slice::IndexType;
    quantities::Vector{Symbol}=[:gas_mass],
    types::Vector{Symbol}=[:cells],
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    theme::Attributes=Theme(),
    title::Union{Symbol,<:AbstractString}="",
    annotation::AbstractString="",
    colorbar::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid for the heatmap
    grid_hm = CubicGrid(box_size, resolution)

    # Set up the grid for the velocity field
    grid_vf = SquareGrid(box_size, 25)

    pf_kwargs = isnothing(colorrange) ? [(;), (;)] : [(; colorrange), (;)]

    for (i, quantity) in pairs(quantities)

        filter_function, translation, rotation, request = selectFilter(
            filter_mode,
            plotParams(quantity).request,
        )

        if quantity ∈ [
            :gas_mass,
            :hydrogen_mass,
            :molecular_mass,
            :br_molecular_mass,
            :atomic_mass,
            :ionized_mass,
            :neutral_mass,
        ]
            component = :gas
        elseif quantity == :stellar_mass
            component = :stars
        elseif quantity == :dm_mass
            component = :halo
        elseif quantity == :bh_mass
            component = :black_hole
        else
            throw(ArgumentError("densityMapVelField: `quantities` contains :$(quantity), \
            which is not a valid symbol. See the documentation for valid options."))
        end

        for simulation_path in simulation_paths

            # Get the simulation name as a string
            sim_name = basename(simulation_path)

            for projection_plane in projection_planes

                # Construct the file name
                base_filename = "$(sim_name)_$(quantity)_$(projection_plane)_density_map"

                plotSnapshot(
                    [simulation_path, simulation_path],
                    request,
                    [heatmap!, arrows!];
                    pf_kwargs,
                    # `plotSnapshot` configuration
                    output_path,
                    base_filename,
                    output_format=".png",
                    show_progress=true,
                    # Data manipulation options
                    slice=iszero(slice) ? (:) : slice,
                    filter_function,
                    da_functions=[daDensity2DProjection, daVelocityField],
                    da_args=[(grid_hm, quantity, ring(types, i)), (grid_vf, component)],
                    da_kwargs=[(; projection_plane), (; projection_plane)],
                    post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                    pp_args=(annotation,),
                    pp_kwargs=(; color=:white),
                    transform_box=true,
                    translation,
                    rotation,
                    smooth=0,
                    x_unit=u"kpc",
                    y_unit=u"kpc",
                    x_exp_factor=0,
                    y_exp_factor=0,
                    x_trim=(-Inf, Inf),
                    y_trim=(-Inf, Inf),
                    x_edges=false,
                    y_edges=false,
                    x_func=identity,
                    y_func=identity,
                    # Axes options
                    xaxis_label="auto_label",
                    yaxis_label="auto_label",
                    xaxis_var_name=string(projection_plane)[1:1],
                    yaxis_var_name=string(projection_plane)[2:2],
                    xaxis_scale_func=identity,
                    yaxis_scale_func=identity,
                    # Plotting options
                    save_figures=!iszero(slice),
                    backup_results=iszero(slice),
                    theme=merge(
                        theme,
                        Theme(
                            size=colorbar ? (880, 680) : (880, 880),
                            figure_padding=(1, 50, 1, 1),
                            Axis=(limits=(-limit, limit, -limit, limit),),
                            Colorbar=(
                                label=L"\mathrm{log}_{10} \Sigma \,\, [\mathrm{M_\odot \, kpc^{-2}}]",
                            ),
                        ),
                    ),
                    sim_labels=nothing,
                    title,
                    colorbar,
                    # Animation options
                    animation=iszero(slice),
                    animation_filename="$(base_filename).mp4",
                    framerate=5,
                )

            end

        end

    end

    return nothing

end

"""
    metallicityMap(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D histogram of the metallicity.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored. If set to 0, an animation using every snapshots will be made.
  - `components::Vector{Symbol}=[:gas]`: Target component. It can be either `:stars` or `:gas`.
  - `types::Vector{Symbol}=[:cells]`: List of component types for the metallicity fields, each element can be either `:particles` or Voronoi `:cells`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Pixel (bin of the 2D histogram) side length.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the values to be plotted.
  - `da_ff::Function=filterNothing`: Filter function for the data analysis function. It must be a function with the signature:

    `da_ff(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `da_ff` filter function.
"""
function metallicityMap(
    simulation_paths::Vector{String},
    slice::IndexType;
    components::Vector{Symbol}=[:gas],
    types::Vector{Symbol}=[:cells],
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    theme::Attributes=Theme(),
    title::Union{Symbol,<:AbstractString}="",
    annotation::AbstractString="",
    colorbar::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid
    grid = CubicGrid(box_size, resolution)

    pf_kwargs = isnothing(colorrange) ? [(;)] : [(; colorrange)]

    for (i, component) in pairs(components)

        if component == :gas

            filter_function, translation, rotation, request = selectFilter(
                filter_mode,
                mergeRequests(
                    plotParams(:gas_metallicity).request,
                    plotParams(:gas_mass_density).request,
                    ff_request,
                ),
            )

        elseif component == :stars

            filter_function, translation, rotation, request = selectFilter(
                filter_mode,
                mergeRequests(plotParams(:stellar_metallicity).request, ff_request),
            )

        else

            throw(ArgumentError("metallicityMap: I don't recognize the component \
            :$(component)"))

        end

        for simulation_path in simulation_paths

            # Get the simulation name as a string
            sim_name = basename(simulation_path)

            for projection_plane in projection_planes

                # Construct the file name
                base_filename = "$(sim_name)_$(component)_$(projection_plane)_metallicity_map"

                plotSnapshot(
                    [simulation_path],
                    request,
                    [heatmap!];
                    pf_kwargs,
                    # `plotSnapshot` configuration
                    output_path,
                    base_filename,
                    output_format=".png",
                    show_progress=true,
                    # Data manipulation options
                    slice=iszero(slice) ? (:) : slice,
                    filter_function,
                    da_functions=[daMetallicity2DProjection],
                    da_args=[(grid, component, ring(types, i))],
                    da_kwargs=[(; element=:all, projection_plane, filter_function=da_ff)],
                    post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                    pp_args=(annotation,),
                    pp_kwargs=(; color=:white),
                    transform_box=true,
                    translation,
                    rotation,
                    smooth=0,
                    x_unit=u"kpc",
                    y_unit=u"kpc",
                    x_exp_factor=0,
                    y_exp_factor=0,
                    x_trim=(-Inf, Inf),
                    y_trim=(-Inf, Inf),
                    x_edges=false,
                    y_edges=false,
                    x_func=identity,
                    y_func=identity,
                    # Axes options
                    xaxis_label="auto_label",
                    yaxis_label="auto_label",
                    xaxis_var_name=string(projection_plane)[1:1],
                    yaxis_var_name=string(projection_plane)[2:2],
                    xaxis_scale_func=identity,
                    yaxis_scale_func=identity,
                    # Plotting options
                    save_figures=!iszero(slice),
                    backup_results=iszero(slice),
                    theme=merge(
                        theme,
                        Theme(
                            size=colorbar ? (880, 760) : (880, 880),
                            figure_padding=(5, 20, 20, 10),
                            Axis=(limits=(-limit, limit, -limit, limit),),
                        ),
                    ),
                    sim_labels=nothing,
                    title,
                    colorbar,
                    # Animation options
                    animation=iszero(slice),
                    animation_filename="$(base_filename).mp4",
                    framerate=5,
                )

            end

        end

    end

    return nothing

end

"""
    temperatureMap(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot a 2D histogram of the temperature.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored. If set to 0, an animation using every snapshots will be made.
  - `type::Symbol=:cells`: If the gas will be assumed to be in `:particles` or in Voronoi `:cells`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `projection_planes::Vector{Symbol}=[:xy]`: Projection planes. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
  - `box_size::Unitful.Length=100u"kpc"`: Physical side length of the plot window.
  - `pixel_length::Unitful.Length=0.1u"kpc"`: Pixel (bin of the 2D histogram) side length.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `annotation::AbstractString=""`: Text to be added into the top left corner of the plot. If left empty, nothing is printed.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing`: Sets the start and end points of the colormap. Use `nothing` to use the extrema of the values to be plotted.
"""
function temperatureMap(
    simulation_paths::Vector{String},
    slice::IndexType;
    type::Symbol=:cells,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    projection_planes::Vector{Symbol}=[:xy],
    box_size::Unitful.Length=100u"kpc",
    pixel_length::Unitful.Length=0.1u"kpc",
    theme::Attributes=Theme(),
    title::Union{Symbol,<:AbstractString}="",
    annotation::AbstractString="",
    colorbar::Bool=false,
    colorrange::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
)::Nothing

    # Compute the axes limits, to avoid white padding around the heatmap grid
    limit = ustrip(u"kpc", box_size / 2.0)

    # Compute number of pixel per side
    resolution = round(Int, box_size / pixel_length)

    # Set up the grid
    grid = CubicGrid(box_size, resolution)

    pf_kwargs = isnothing(colorrange) ? [(;)] : [(; colorrange)]

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(plotParams(:temperature).request, plotParams(:gas_mass_density).request),
    )

    for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        for projection_plane in projection_planes

            # Construct the file name
            base_filename = "$(sim_name)_$(projection_plane)_temperature_map"

            plotSnapshot(
                [simulation_path],
                request,
                [heatmap!];
                pf_kwargs,
                # `plotSnapshot` configuration
                output_path,
                base_filename,
                output_format=".png",
                show_progress=true,
                # Data manipulation options
                slice=iszero(slice) ? (:) : slice,
                filter_function,
                da_functions=[daTemperature2DProjection],
                da_args=[(grid, type)],
                da_kwargs=[(; projection_plane)],
                post_processing=isempty(annotation) ? getNothing : ppAnnotation!,
                pp_args=(annotation,),
                pp_kwargs=(; color=:white),
                transform_box=true,
                translation,
                rotation,
                smooth=0,
                x_unit=u"kpc",
                y_unit=u"kpc",
                x_exp_factor=0,
                y_exp_factor=0,
                x_trim=(-Inf, Inf),
                y_trim=(-Inf, Inf),
                x_edges=false,
                y_edges=false,
                x_func=identity,
                y_func=identity,
                # Axes options
                xaxis_label="auto_label",
                yaxis_label="auto_label",
                xaxis_var_name=string(projection_plane)[1:1],
                yaxis_var_name=string(projection_plane)[2:2],
                xaxis_scale_func=identity,
                yaxis_scale_func=identity,
                # Plotting options
                save_figures=!iszero(slice),
                backup_results=iszero(slice),
                theme=merge(
                    theme,
                    Theme(
                        size=colorbar ? (880, 760) : (880, 880),
                        figure_padding=(5, 20, 20, 10),
                        Axis=(limits=(-limit, limit, -limit, limit),),
                    ),
                ),
                sim_labels=nothing,
                title,
                colorbar,
                # Animation options
                animation=iszero(slice),
                animation_filename="$(base_filename).mp4",
                framerate=5,
            )

        end

    end

    return nothing

end

"""
    scatterPlot(
        simulation_paths::Vector{String},
        slice::IndexType,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot two quantities as a scatter plot, one marker for every cell/particle.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`                -> Stellar mass.
      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:dm_mass`                     -> Dark matter mass.
      + `:bh_mass`                     -> Black hole mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`            -> Stellar gas mass (according to out SF model).
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:ionized_neutral_fraction`    -> Fraction of ionized gas to neutral gas.
      + `:stellar_gas_fraction`        -> Stellar gas fraction (according to out SF model).
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_td`                      -> Total gas depletion time.
      + `:molecular_td`                -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time.
      + `:br_molecular_td`             -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_td`                   -> Atomic hydrogen (``\\mathrm{HI}``) depletion time.
      + `:ionized_td`                  -> Ionized hydrogen (``\\mathrm{HII}``) depletion time.
      + `:neutral_td`                  -> Neutral hydrogen (``\\mathrm{HI + H_2}``) depletion time.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`         -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`         -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`     -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`          -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`         -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`              -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`         -> Stellar circularity.
      + `:stellar_vcirc`               -> Stellar circular velocity.
      + `:stellar_vradial`             -> Stellar radial speed.
      + `:stellar_vtangential`         -> Stellar tangential speed.
      + `:stellar_vzstar`              -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                 -> Stellar age.
      + `:sfr`                         -> The star formation rate.
      + `:ssfr`                        -> The specific star formation rate.
      + `:observational_sfr`           -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`          -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:stellar_eff`                 -> The star formation efficiency per free-fall time for the gas that has turn into stars.
      + `:gas_eff`                     -> The star formation efficiency per free-fall time for the gas.
      + `:molecular_eff`               -> The star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas.
      + `:br_molecular_eff`            -> The star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_eff`                  -> The star formation efficiency per free-fall time for the atomic hydrogen (``\\mathrm{HI}``) gas.
      + `:ionized_eff`                 -> The star formation efficiency per free-fall time for the ionized hydrogen (``\\mathrm{HII}``) gas.
      + `:neutral_eff`                 -> The star formation efficiency per free-fall time for the neutral hydrogen (``\\mathrm{HI + H_2}``) gas.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
      + `:ode_gas_it`                  -> Integration time.
      + `:ode_gas_accu_it`             -> Accumulated integration time.
      + `:ode_gas_tau_s`               -> Star formation time scale, ``\\tau_\\mathrm{S}``.
      + `:ode_gas_eta_d`               -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``.
      + `:ode_gas_eta_i`               -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``.
      + `:ode_gas_r`                   -> Mass recycling parameter, ``R``.
      + `:ode_gas_cold_mf`             -> Cold gas mass fraction.
      + `:ode_stellar_it`              -> Integration time, for the gas that form the stars.
      + `:ode_stellar_accu_it`         -> Accumulated integration time, for the gas that form the stars.
      + `:ode_stellar_tau_s`           -> Star formation time scale, ``\\tau_\\mathrm{S}``, for the gas that form the stars.
      + `:ode_stellar_eta_d`           -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``, for the gas that form the stars.
      + `:ode_stellar_eta_i`           -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``, for the gas that form the stars.
      + `:ode_stellar_r`               -> Mass recycling parameter, ``R``, for the gas that form the stars.
      + `:ode_stellar_cold_mf`         -> Cold gas mass fraction, for the gas that form the stars.
      + `:ode_stellar_gas_rho`         -> Gas mass density, for the gas that form the stars.
      + `:ode_stellar_gas_Z`           -> Gas metallicity, for the gas that form the stars (solar units).
      + `:ode_stellar_gas_mass`        -> Cell mass, for the gas that form the stars.
      + `:ode_stellar_gas_sfr`         -> SFR associated to the gas particles/cells within the code, for the gas that form the stars.
      + `:ode_stellar_gas_P`           -> Gas pressure, for the gas that form the stars.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`                -> Stellar mass.
      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:dm_mass`                     -> Dark matter mass.
      + `:bh_mass`                     -> Black hole mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`            -> Stellar gas mass (according to out SF model).
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:ionized_neutral_fraction`    -> Fraction of ionized gas to neutral gas.
      + `:stellar_gas_fraction`        -> Stellar gas fraction (according to out SF model).
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_td`                      -> Total gas depletion time.
      + `:molecular_td`                -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time.
      + `:br_molecular_td`             -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_td`                   -> Atomic hydrogen (``\\mathrm{HI}``) depletion time.
      + `:ionized_td`                  -> Ionized hydrogen (``\\mathrm{HII}``) depletion time.
      + `:neutral_td`                  -> Neutral hydrogen (``\\mathrm{HI + H_2}``) depletion time.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`         -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`         -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`     -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`          -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`         -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`              -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`         -> Stellar circularity.
      + `:stellar_vcirc`               -> Stellar circular velocity.
      + `:stellar_vradial`             -> Stellar radial speed.
      + `:stellar_vtangential`         -> Stellar tangential speed.
      + `:stellar_vzstar`              -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                 -> Stellar age.
      + `:sfr`                         -> The star formation rate.
      + `:ssfr`                        -> The specific star formation rate.
      + `:observational_sfr`           -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`          -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:stellar_eff`                 -> The star formation efficiency per free-fall time for the gas that has turn into stars.
      + `:gas_eff`                     -> The star formation efficiency per free-fall time for the gas.
      + `:molecular_eff`               -> The star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas.
      + `:br_molecular_eff`            -> The star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_eff`                  -> The star formation efficiency per free-fall time for the atomic hydrogen (``\\mathrm{HI}``) gas.
      + `:ionized_eff`                 -> The star formation efficiency per free-fall time for the ionized hydrogen (``\\mathrm{HII}``) gas.
      + `:neutral_eff`                 -> The star formation efficiency per free-fall time for the neutral hydrogen (``\\mathrm{HI + H_2}``) gas.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
      + `:ode_gas_it`                  -> Integration time.
      + `:ode_gas_accu_it`             -> Accumulated integration time.
      + `:ode_gas_tau_s`               -> Star formation time scale, ``\\tau_\\mathrm{S}``.
      + `:ode_gas_eta_d`               -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``.
      + `:ode_gas_eta_i`               -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``.
      + `:ode_gas_r`                   -> Mass recycling parameter, ``R``.
      + `:ode_gas_cold_mf`             -> Cold gas mass fraction.
      + `:ode_stellar_it`              -> Integration time, for the gas that form the stars.
      + `:ode_stellar_accu_it`         -> Accumulated integration time, for the gas that form the stars.
      + `:ode_stellar_tau_s`           -> Star formation time scale, ``\\tau_\\mathrm{S}``, for the gas that form the stars.
      + `:ode_stellar_eta_d`           -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``, for the gas that form the stars.
      + `:ode_stellar_eta_i`           -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``, for the gas that form the stars.
      + `:ode_stellar_r`               -> Mass recycling parameter, ``R``, for the gas that form the stars.
      + `:ode_stellar_cold_mf`         -> Cold gas mass fraction, for the gas that form the stars.
      + `:ode_stellar_gas_rho`         -> Gas mass density, for the gas that form the stars.
      + `:ode_stellar_gas_Z`           -> Gas metallicity, for the gas that form the stars (solar units).
      + `:ode_stellar_gas_mass`        -> Cell mass, for the gas that form the stars.
      + `:ode_stellar_gas_sfr`         -> SFR associated to the gas particles/cells within the code, for the gas that form the stars.
      + `:ode_stellar_gas_P`           -> Gas pressure, for the gas that form the stars.
  - `xlog::Bool=false`: If true, sets everything so the x axis is log10(`x_quantity`).
  - `ylog::Bool=false`: If true, sets everything so the y axis is log10(`y_quantity`).
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `da_ff::Function=filterNothing`: A function with the signature:

    `da_ff(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `da_ff` filter function.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function scatterPlot(
    simulation_paths::Vector{String},
    slice::IndexType,
    x_quantity::Symbol,
    y_quantity::Symbol;
    xlog::Bool=false,
    ylog::Bool=false,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(x_plot_params.request, y_plot_params.request, ff_request),
    )

    # Set arguments for a log x axis
    if xlog
        x_log        = x_plot_params.unit
        x_unit       = Unitful.NoUnits
        unit_label   = getUnitLabel(0, x_plot_params.unit; latex=true)
        if isempty(unit_label)
            xaxis_label  = L"$\log_{10} \, $auto_label"
        else
            xaxis_label  = L"$\log_{10} \, $auto_label [%$(unit_label)]"
        end
        x_exp_factor = 0
    else
        x_log        = nothing
        x_unit       = x_plot_params.unit
        xaxis_label  = x_plot_params.axis_label
        x_exp_factor = x_plot_params.exp_factor
    end

    # Set arguments for a log y axis
    if ylog
        y_log        = y_plot_params.unit
        y_unit       = Unitful.NoUnits
        unit_label   = getUnitLabel(0, y_plot_params.unit; latex=true)
        if isempty(unit_label)
            yaxis_label  = L"$\log_{10} \, $auto_label"
        else
            yaxis_label  = L"$\log_{10} \, $auto_label [%$(unit_label)]"
        end
        y_exp_factor = 0
    else
        y_log        = nothing
        y_unit       = y_plot_params.unit
        yaxis_label  = y_plot_params.axis_label
        y_exp_factor = y_plot_params.exp_factor
    end

    for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        plotSnapshot(
            [simulation_path],
            request,
            [scatter!];
            pf_kwargs=[(; markersize=2)],
            # `plotSnapshot` configuration
            output_path,
            base_filename="$(sim_name)_$(y_quantity)_vs_$(x_quantity)",
            output_format=".png",
            show_progress=true,
            # Data manipulation options
            slice=slice,
            filter_function,
            da_functions=[daScatterGalaxy],
            da_args=[(x_quantity, y_quantity)],
            da_kwargs=[(; x_log, y_log, filter_function=da_ff)],
            post_processing=getNothing,
            pp_args=(),
            pp_kwargs=(;),
            transform_box=true,
            translation,
            rotation,
            smooth=0,
            x_unit,
            y_unit,
            x_exp_factor,
            y_exp_factor,
            x_trim=(-Inf, Inf),
            y_trim=(-Inf, Inf),
            x_edges=false,
            y_edges=false,
            x_func=identity,
            y_func=identity,
            # Axes options
            xaxis_label,
            yaxis_label,
            xaxis_var_name=x_plot_params.var_name,
            yaxis_var_name=y_plot_params.var_name,
            xaxis_scale_func=identity,
            yaxis_scale_func=identity,
            # Plotting and animation options
            save_figures=true,
            backup_results=false,
            theme,
            sim_labels=nothing,
            title="",
            colorbar=false,
            # Animation options
            animation=false,
            animation_filename="animation.mp4",
            framerate=10,
        )

    end

    return nothing

end

"""
    scatterDensityMap(
        simulation_paths::Vector{String},
        slice::IndexType,
        x_quantity::Symbol,
        y_quantity::Symbol,
        z_quantity::Symbol,
        z_unit::Unitful.Units;
        <keyword arguments>
    )::Nothing

Plot two quantities as a density scatter plot (2D histogram), weighted by `z_quantity`.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`                -> Stellar mass.
      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:dm_mass`                     -> Dark matter mass.
      + `:bh_mass`                     -> Black hole mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`            -> Stellar gas mass (according to out SF model).
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:ionized_neutral_fraction`    -> Fraction of ionized gas to neutral gas.
      + `:stellar_gas_fraction`        -> Stellar gas fraction (according to out SF model).
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_td`                      -> Total gas depletion time.
      + `:molecular_td`                -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time.
      + `:br_molecular_td`             -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_td`                   -> Atomic hydrogen (``\\mathrm{HI}``) depletion time.
      + `:ionized_td`                  -> Ionized hydrogen (``\\mathrm{HII}``) depletion time.
      + `:neutral_td`                  -> Neutral hydrogen (``\\mathrm{HI + H_2}``) depletion time.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`         -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`         -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`     -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`          -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`         -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`              -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`         -> Stellar circularity.
      + `:stellar_vcirc`               -> Stellar circular velocity.
      + `:stellar_vradial`             -> Stellar radial speed.
      + `:stellar_vtangential`         -> Stellar tangential speed.
      + `:stellar_vzstar`              -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                 -> Stellar age.
      + `:sfr`                         -> The star formation rate.
      + `:ssfr`                        -> The specific star formation rate.
      + `:observational_sfr`           -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`          -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:stellar_eff`                 -> The star formation efficiency per free-fall time for the gas that has turn into stars.
      + `:gas_eff`                     -> The star formation efficiency per free-fall time for the gas.
      + `:molecular_eff`               -> The star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas.
      + `:br_molecular_eff`            -> The star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_eff`                  -> The star formation efficiency per free-fall time for the atomic hydrogen (``\\mathrm{HI}``) gas.
      + `:ionized_eff`                 -> The star formation efficiency per free-fall time for the ionized hydrogen (``\\mathrm{HII}``) gas.
      + `:neutral_eff`                 -> The star formation efficiency per free-fall time for the neutral hydrogen (``\\mathrm{HI + H_2}``) gas.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
      + `:ode_gas_it`                  -> Integration time.
      + `:ode_gas_accu_it`             -> Accumulated integration time.
      + `:ode_gas_tau_s`               -> Star formation time scale, ``\\tau_\\mathrm{S}``.
      + `:ode_gas_eta_d`               -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``.
      + `:ode_gas_eta_i`               -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``.
      + `:ode_gas_r`                   -> Mass recycling parameter, ``R``.
      + `:ode_gas_cold_mf`             -> Cold gas mass fraction.
      + `:ode_stellar_it`              -> Integration time, for the gas that form the stars.
      + `:ode_stellar_accu_it`         -> Accumulated integration time, for the gas that form the stars.
      + `:ode_stellar_tau_s`           -> Star formation time scale, ``\\tau_\\mathrm{S}``, for the gas that form the stars.
      + `:ode_stellar_eta_d`           -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``, for the gas that form the stars.
      + `:ode_stellar_eta_i`           -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``, for the gas that form the stars.
      + `:ode_stellar_r`               -> Mass recycling parameter, ``R``, for the gas that form the stars.
      + `:ode_stellar_cold_mf`         -> Cold gas mass fraction, for the gas that form the stars.
      + `:ode_stellar_gas_rho`         -> Gas mass density, for the gas that form the stars.
      + `:ode_stellar_gas_Z`           -> Gas metallicity, for the gas that form the stars (solar units).
      + `:ode_stellar_gas_mass`        -> Cell mass, for the gas that form the stars.
      + `:ode_stellar_gas_sfr`         -> SFR associated to the gas particles/cells within the code, for the gas that form the stars.
      + `:ode_stellar_gas_P`           -> Gas pressure, for the gas that form the stars.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`                -> Stellar mass.
      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:dm_mass`                     -> Dark matter mass.
      + `:bh_mass`                     -> Black hole mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`            -> Stellar gas mass (according to out SF model).
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:ionized_neutral_fraction`    -> Fraction of ionized gas to neutral gas.
      + `:stellar_gas_fraction`        -> Stellar gas fraction (according to out SF model).
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_td`                      -> Total gas depletion time.
      + `:molecular_td`                -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time.
      + `:br_molecular_td`             -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_td`                   -> Atomic hydrogen (``\\mathrm{HI}``) depletion time.
      + `:ionized_td`                  -> Ionized hydrogen (``\\mathrm{HII}``) depletion time.
      + `:neutral_td`                  -> Neutral hydrogen (``\\mathrm{HI + H_2}``) depletion time.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`         -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`         -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`     -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`          -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`         -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`              -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`         -> Stellar circularity.
      + `:stellar_vcirc`               -> Stellar circular velocity.
      + `:stellar_vradial`             -> Stellar radial speed.
      + `:stellar_vtangential`         -> Stellar tangential speed.
      + `:stellar_vzstar`              -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                 -> Stellar age.
      + `:sfr`                         -> The star formation rate.
      + `:ssfr`                        -> The specific star formation rate.
      + `:observational_sfr`           -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`          -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:stellar_eff`                 -> The star formation efficiency per free-fall time for the gas that has turn into stars.
      + `:gas_eff`                     -> The star formation efficiency per free-fall time for the gas.
      + `:molecular_eff`               -> The star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas.
      + `:br_molecular_eff`            -> The star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_eff`                  -> The star formation efficiency per free-fall time for the atomic hydrogen (``\\mathrm{HI}``) gas.
      + `:ionized_eff`                 -> The star formation efficiency per free-fall time for the ionized hydrogen (``\\mathrm{HII}``) gas.
      + `:neutral_eff`                 -> The star formation efficiency per free-fall time for the neutral hydrogen (``\\mathrm{HI + H_2}``) gas.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
      + `:ode_gas_it`                  -> Integration time.
      + `:ode_gas_accu_it`             -> Accumulated integration time.
      + `:ode_gas_tau_s`               -> Star formation time scale, ``\\tau_\\mathrm{S}``.
      + `:ode_gas_eta_d`               -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``.
      + `:ode_gas_eta_i`               -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``.
      + `:ode_gas_r`                   -> Mass recycling parameter, ``R``.
      + `:ode_gas_cold_mf`             -> Cold gas mass fraction.
      + `:ode_stellar_it`              -> Integration time, for the gas that form the stars.
      + `:ode_stellar_accu_it`         -> Accumulated integration time, for the gas that form the stars.
      + `:ode_stellar_tau_s`           -> Star formation time scale, ``\\tau_\\mathrm{S}``, for the gas that form the stars.
      + `:ode_stellar_eta_d`           -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``, for the gas that form the stars.
      + `:ode_stellar_eta_i`           -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``, for the gas that form the stars.
      + `:ode_stellar_r`               -> Mass recycling parameter, ``R``, for the gas that form the stars.
      + `:ode_stellar_cold_mf`         -> Cold gas mass fraction, for the gas that form the stars.
      + `:ode_stellar_gas_rho`         -> Gas mass density, for the gas that form the stars.
      + `:ode_stellar_gas_Z`           -> Gas metallicity, for the gas that form the stars (solar units).
      + `:ode_stellar_gas_mass`        -> Cell mass, for the gas that form the stars.
      + `:ode_stellar_gas_sfr`         -> SFR associated to the gas particles/cells within the code, for the gas that form the stars.
      + `:ode_stellar_gas_P`           -> Gas pressure, for the gas that form the stars.
  - `z_quantity::Symbol`: Quantity for the z axis (weights). The options are:

      + `:stellar_mass`                -> Stellar mass.
      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:dm_mass`                     -> Dark matter mass.
      + `:bh_mass`                     -> Black hole mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`            -> Stellar gas mass (according to out SF model).
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:ionized_neutral_fraction`    -> Fraction of ionized gas to neutral gas.
      + `:stellar_gas_fraction`        -> Stellar gas fraction (according to out SF model).
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_td`                      -> Total gas depletion time.
      + `:molecular_td`                -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time.
      + `:br_molecular_td`             -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_td`                   -> Atomic hydrogen (``\\mathrm{HI}``) depletion time.
      + `:ionized_td`                  -> Ionized hydrogen (``\\mathrm{HII}``) depletion time.
      + `:neutral_td`                  -> Neutral hydrogen (``\\mathrm{HI + H_2}``) depletion time.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`         -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`         -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`     -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`          -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`         -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`              -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`         -> Stellar circularity.
      + `:stellar_vcirc`               -> Stellar circular velocity.
      + `:stellar_vradial`             -> Stellar radial speed.
      + `:stellar_vtangential`         -> Stellar tangential speed.
      + `:stellar_vzstar`              -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                 -> Stellar age.
      + `:sfr`                         -> The star formation rate.
      + `:ssfr`                        -> The specific star formation rate.
      + `:observational_sfr`           -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`          -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:stellar_eff`                 -> The star formation efficiency per free-fall time for the gas that has turn into stars.
      + `:gas_eff`                     -> The star formation efficiency per free-fall time for the gas.
      + `:molecular_eff`               -> The star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas.
      + `:br_molecular_eff`            -> The star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_eff`                  -> The star formation efficiency per free-fall time for the atomic hydrogen (``\\mathrm{HI}``) gas.
      + `:ionized_eff`                 -> The star formation efficiency per free-fall time for the ionized hydrogen (``\\mathrm{HII}``) gas.
      + `:neutral_eff`                 -> The star formation efficiency per free-fall time for the neutral hydrogen (``\\mathrm{HI + H_2}``) gas.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
      + `:ode_gas_it`                  -> Integration time.
      + `:ode_gas_accu_it`             -> Accumulated integration time.
      + `:ode_gas_tau_s`               -> Star formation time scale, ``\\tau_\\mathrm{S}``.
      + `:ode_gas_eta_d`               -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``.
      + `:ode_gas_eta_i`               -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``.
      + `:ode_gas_r`                   -> Mass recycling parameter, ``R``.
      + `:ode_gas_cold_mf`             -> Cold gas mass fraction.
      + `:ode_stellar_it`              -> Integration time, for the gas that form the stars.
      + `:ode_stellar_accu_it`         -> Accumulated integration time, for the gas that form the stars.
      + `:ode_stellar_tau_s`           -> Star formation time scale, ``\\tau_\\mathrm{S}``, for the gas that form the stars.
      + `:ode_stellar_eta_d`           -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``, for the gas that form the stars.
      + `:ode_stellar_eta_i`           -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``, for the gas that form the stars.
      + `:ode_stellar_r`               -> Mass recycling parameter, ``R``, for the gas that form the stars.
      + `:ode_stellar_cold_mf`         -> Cold gas mass fraction, for the gas that form the stars.
      + `:ode_stellar_gas_rho`         -> Gas mass density, for the gas that form the stars.
      + `:ode_stellar_gas_Z`           -> Gas metallicity, for the gas that form the stars (solar units).
      + `:ode_stellar_gas_mass`        -> Cell mass, for the gas that form the stars.
      + `:ode_stellar_gas_sfr`         -> SFR associated to the gas particles/cells within the code, for the gas that form the stars.
      + `:ode_stellar_gas_P`           -> Gas pressure, for the gas that form the stars.
  - `z_unit::Unitful.Units`: Target unit for the z axis.
  - `x_range::Union{NTuple{2,<:Number},Nothing}=nothing`: x axis range. If set to `nothing`, the extrema of the values will be used.
  - `y_range::Union{NTuple{2,<:Number},Nothing}=nothing`: y axis range. If set to `nothing`, the extrema of the values will be used.
  - `xlog::Bool=false`: If true, sets everything so the x axis is log10(`x_quantity`).
  - `ylog::Bool=false`: If true, sets everything so the y axis is log10(`y_quantity`).
  - `total::Bool=true`: If the sum (default) or the mean of `z_quantity` will be used as the value of each pixel.
  - `n_bins::Int=100`: Number of bins per side of the square grid.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `da_ff::Function=filterNothing`: A function with the signature:

    `da_ff(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `da_ff` filter function.
  - `title::Union{Symbol,<:AbstractString}=""`: Title for the figure. If left empty, no title is printed. It can also be set to one of the following options:

      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
      + `:scale_factor`  -> Scale factor (only relevant for cosmological simulations).
      + `:redshift`      -> Redshift (only relevant for cosmological simulations).
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function scatterDensityMap(
    simulation_paths::Vector{String},
    slice::IndexType,
    x_quantity::Symbol,
    y_quantity::Symbol,
    z_quantity::Symbol,
    z_unit::Unitful.Units;
    x_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    y_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    xlog::Bool=false,
    ylog::Bool=false,
    total::Bool=true,
    n_bins::Int=100,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    title::Union{Symbol,<:AbstractString}="",
    colorbar::Bool=false,
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)
    z_plot_params = plotParams(z_quantity)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(
            x_plot_params.request,
            y_plot_params.request,
            z_plot_params.request,
            ff_request,
        ),
    )

    n_sims = length(simulation_paths)

    # Set arguments for a log x axis
    if xlog
        x_log        = x_plot_params.unit
        x_unit       = Unitful.NoUnits
        unit_label   = getUnitLabel(0, x_plot_params.unit; latex=true)
        if isempty(unit_label)
            xaxis_label  = L"$\log_{10} \, $auto_label"
        else
            xaxis_label  = L"$\log_{10} \, $auto_label [%$(unit_label)]"
        end
        x_exp_factor = 0
    else
        x_log        = nothing
        x_unit       = x_plot_params.unit
        xaxis_label  = x_plot_params.axis_label
        x_exp_factor = x_plot_params.exp_factor
    end

    # Set arguments for a log y axis
    if ylog
        y_log        = y_plot_params.unit
        y_unit       = Unitful.NoUnits
        unit_label   = getUnitLabel(0, y_plot_params.unit; latex=true)
        if isempty(unit_label)
            yaxis_label  = L"$\log_{10} \, $auto_label"
        else
            yaxis_label  = L"$\log_{10} \, $auto_label [%$(unit_label)]"
        end
        y_exp_factor = 0
    else
        y_log        = nothing
        y_unit       = y_plot_params.unit
        yaxis_label  = y_plot_params.axis_label
        y_exp_factor = y_plot_params.exp_factor
    end

    for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        if isone(n_sims)
            base_filename="$(y_quantity)_vs_$(x_quantity)"
        else
            base_filename="$(sim_name)_$(y_quantity)_vs_$(x_quantity)"
        end

        plotSnapshot(
            [simulation_path],
            request,
            [heatmap!];
            pf_kwargs=[(;)],
            # `plotSnapshot` configuration
            output_path,
            base_filename,
            output_format=".png",
            show_progress=true,
            # Data manipulation options
            slice,
            filter_function,
            da_functions=[daScatterWeightedDensity],
            da_args=[(x_quantity, y_quantity, z_quantity, z_unit)],
            da_kwargs=[
                (;
                    x_range,
                    y_range,
                    x_log,
                    y_log,
                    total,
                    n_bins,
                    filter_function=da_ff,
                ),
            ],
            post_processing=getNothing,
            pp_args=(),
            pp_kwargs=(;),
            transform_box=true,
            translation,
            rotation,
            smooth=0,
            x_unit,
            y_unit,
            x_exp_factor,
            y_exp_factor,
            x_trim=(-Inf, Inf),
            y_trim=(-Inf, Inf),
            x_edges=false,
            y_edges=false,
            x_func=identity,
            y_func=identity,
            # Axes options
            xaxis_label,
            yaxis_label,
            xaxis_var_name=x_plot_params.var_name,
            yaxis_var_name=y_plot_params.var_name,
            xaxis_scale_func=identity,
            yaxis_scale_func=identity,
            # Plotting and animation options
            save_figures=true,
            backup_results=false,
            theme,
            sim_labels=nothing,
            title,
            colorbar,
            # Animation options
            animation=false,
            animation_filename="animation.mp4",
            framerate=10,
        )

    end

    return nothing

end

"""
    atomicMolecularTransition(
        simulation_paths::Vector{String},
        slice::IndexType,
        ranges::Vector{<:Tuple{<:Real,<:Real}};
        <keyword arguments>
    )::Nothing

Plot the atomic gas to molecular gas transition for a set of metallicity ranges.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `ranges::Vector{<:Tuple{<:Real,<:Real}}`: Metallicity (as in the fractional mass of metals) ranges.
  - `plot_type::Symbol=:heatmap`: Type of plot. The options are:

      + `:heatmap` -> Heatmap. One figure per range will be produced.
      + `:scatter` -> Scatter plot. A single figure with every range will be produced.
  - `output_path::String="./"`: Path to the output folder.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function atomicMolecularTransition(
    simulation_paths::Vector{String},
    slice::IndexType,
    ranges::Vector{<:Tuple{<:Real,<:Real}};
    plot_type::Symbol=:heatmap,
    output_path::String="./",
    theme::Attributes=Theme(),
)::Nothing

    # Set some plotting parameters
    x_quantity = :atomic_number_density
    y_quantity = :molecular_neutral_fraction

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    filter_function, translation, rotation, request = selectFilter(
        :subhalo,
        mergeRequests(
            x_plot_params.request,
            y_plot_params.request,
            Dict(:gas => ["GZ  ", "GMET"], :stars => ["GZ2 ", "GME2"]),
        ),
    )

    for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)
        filename = "$(sim_name)_$(y_quantity)_vs_$(x_quantity)"

        if plot_type == :heatmap

            for range in ranges

                plotSnapshot(
                    fill(simulation_path, length(ranges)),
                    request,
                    [heatmap!];
                    pf_kwargs=[(;)],
                    # `plotSnapshot` configuration
                    output_path,
                    base_filename="$(filename)_$(range[1])_Z_$(range[2])",
                    output_format=".png",
                    show_progress=true,
                    # Data manipulation options
                    slice,
                    filter_function,
                    da_functions=[daScatterDensity],
                    da_args=[(x_quantity, y_quantity)],
                    da_kwargs=[
                        (;
                            x_log=x_plot_params.unit,
                            y_log=y_plot_params.unit,
                            filter_function=dd -> filterByQuantity(
                                dd,
                                :gas_metallicity,
                                :gas,
                                range[1],
                                range[2],
                            )
                        )
                    ],
                    post_processing=getNothing,
                    pp_args=(),
                    pp_kwargs=(;),
                    transform_box=true,
                    translation,
                    rotation,
                    smooth=0,
                    x_unit=Unitful.NoUnits,
                    y_unit=Unitful.NoUnits,
                    x_exp_factor=x_plot_params.exp_factor,
                    y_exp_factor=y_plot_params.exp_factor,
                    x_trim=(-Inf, Inf),
                    y_trim=(-Inf, Inf),
                    x_edges=false,
                    y_edges=false,
                    x_func=identity,
                    y_func=identity,
                    # Axes options
                    xaxis_label=L"$\log_{10} \, $auto_label $[\mathrm{cm}^{-3}]$",
                    yaxis_label=L"$\log_{10} \, $auto_label",
                    xaxis_var_name=x_plot_params.var_name,
                    yaxis_var_name=y_plot_params.var_name,
                    xaxis_scale_func=identity,
                    yaxis_scale_func=identity,
                    # Plotting and animation options
                    save_figures=true,
                    backup_results=false,
                    theme,
                    sim_labels=nothing,
                    title=L"%$(range[1]) \, < \, Z \, < \, %$(range[2])",
                    colorbar=false,
                    # Animation options
                    animation=false,
                    animation_filename="animation.mp4",
                    framerate=10,
                )

            end

        elseif plot_type == :scatter

            plotSnapshot(
                fill(simulation_path, length(ranges)),
                request,
                [scatter!];
                pf_kwargs=[(; markersize=4)],
                # `plotSnapshot` configuration
                output_path,
                base_filename=filename,
                output_format=".png",
                show_progress=true,
                # Data manipulation options
                slice,
                filter_function,
                da_functions=[daScatterGalaxy],
                da_args=[(x_quantity, y_quantity)],
                da_kwargs = [
                    (;
                        filter_function=dd -> filterByQuantity(
                            dd,
                            :gas_metallicity,
                            :gas,
                            range[1],
                            range[2],
                        ),
                    ) for range in ranges
                ],
                post_processing=getNothing,
                pp_args=(),
                pp_kwargs=(;),
                transform_box=true,
                translation,
                rotation,
                smooth=0,
                x_unit=x_plot_params.unit,
                y_unit=y_plot_params.unit,
                x_exp_factor=x_plot_params.exp_factor,
                y_exp_factor=y_plot_params.exp_factor,
                x_trim=(-Inf, Inf),
                y_trim=(-Inf, Inf),
                x_edges=false,
                y_edges=false,
                x_func=identity,
                y_func=identity,
                # Axes options
                xaxis_label=x_plot_params.axis_label,
                yaxis_label=y_plot_params.axis_label,
                xaxis_var_name=x_plot_params.var_name,
                yaxis_var_name=y_plot_params.var_name,
                xaxis_scale_func=log10,
                yaxis_scale_func=log10,
                # Plotting and animation options
                save_figures=true,
                backup_results=false,
                theme=merge(
                    theme,
                    Theme(Legend=(nbanks=1, halign=:left, valign=:top, padding=(0, 0, 0, 15)),),
                ),
                sim_labels= ["$(range[1]) < Z < $(range[2])" for range in ranges],
                title="",
                colorbar=false,
                # Animation options
                animation=false,
                animation_filename="animation.mp4",
                framerate=10,
            )

        else

            throw(ArgumentError("atomicMolecularTransition: `plot_type` can only be :heatmap or \
            :scatter, but I got :$(plot_type)"))

        end

    end

    return nothing

end

"""
    gasBarPlot(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol,
        edges::Vector{<:Number};
        <keyword arguments>
    )::Nothing

Plot a bar plot of the gas fractions, where the bins are a given gas `quantity`.

Only for gas cells that have entered out routine.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: Target quantity. The possibilities are:

      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`            -> Stellar gas mass (according to out SF model).
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:ionized_neutral_fraction`    -> Fraction of ionized gas to neutral gas.
      + `:stellar_gas_fraction`        -> Stellar gas fraction (according to out SF model).
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_td`                      -> Total gas depletion time.
      + `:molecular_td`                -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time.
      + `:br_molecular_td`             -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_td`                   -> Atomic hydrogen (``\\mathrm{HI}``) depletion time.
      + `:ionized_td`                  -> Ionized hydrogen (``\\mathrm{HII}``) depletion time.
      + `:neutral_td`                  -> Neutral hydrogen (``\\mathrm{HI + H_2}``) depletion time.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:stellar_eff`                 -> The star formation efficiency per free-fall time for the gas that has turn into stars.
      + `:gas_eff`                     -> The star formation efficiency per free-fall time for the gas.
      + `:molecular_eff`               -> The star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas.
      + `:br_molecular_eff`            -> The star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_eff`                  -> The star formation efficiency per free-fall time for the atomic hydrogen (``\\mathrm{HI}``) gas.
      + `:ionized_eff`                 -> The star formation efficiency per free-fall time for the ionized hydrogen (``\\mathrm{HII}``) gas.
      + `:neutral_eff`                 -> The star formation efficiency per free-fall time for the neutral hydrogen (``\\mathrm{HI + H_2}``) gas.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
      + `:ode_gas_it`                  -> Integration time.
      + `:ode_gas_accu_it`             -> Accumulated integration time.
      + `:ode_gas_tau_s`               -> Star formation time scale, ``\\tau_\\mathrm{S}``.
      + `:ode_gas_eta_d`               -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``.
      + `:ode_gas_eta_i`               -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``.
      + `:ode_gas_r`                   -> Mass recycling parameter, ``R``.
      + `:ode_gas_cold_mf`             -> Cold gas mass fraction.
      + `:ode_stellar_it`              -> Integration time, for the gas that form the stars.
      + `:ode_stellar_accu_it`         -> Accumulated integration time, for the gas that form the stars.
      + `:ode_stellar_tau_s`           -> Star formation time scale, ``\\tau_\\mathrm{S}``, for the gas that form the stars.
      + `:ode_stellar_eta_d`           -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``, for the gas that form the stars.
      + `:ode_stellar_eta_i`           -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``, for the gas that form the stars.
      + `:ode_stellar_r`               -> Mass recycling parameter, ``R``, for the gas that form the stars.
      + `:ode_stellar_cold_mf`         -> Cold gas mass fraction, for the gas that form the stars.
      + `:ode_stellar_gas_rho`         -> Gas mass density, for the gas that form the stars.
      + `:ode_stellar_gas_Z`           -> Gas metallicity, for the gas that form the stars (solar units).
      + `:ode_stellar_gas_mass`        -> Cell mass, for the gas that form the stars.
      + `:ode_stellar_gas_sfr`         -> SFR associated to the gas particles/cells within the code, for the gas that form the stars.
      + `:ode_stellar_gas_P`           -> Gas pressure, for the gas that form the stars.
  - `edges::Vector{<:Number}`: A sorted list of bin edges for `quantity`.
  - `axis_label::Union{AbstractString,Nothing}=nothing`: Label for the axis. It can contain the string `auto_label`, which will be replaced by the default label: `var_name` / 10^`exp_factor` `unit`. If set to `nothing` a label will be assigned automaticaly.
  - `exp_ticks::Bool=false`: If the axis ticks will be the ``\\log_{10}`` of `edges`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function gasBarPlot(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol,
    edges::Vector{<:Number};
    axis_label::Union{AbstractString,Nothing}=nothing,
    exp_ticks::Bool=false,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(quantity)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(
            plot_params.request,
            plotParams(:molecular_mass).request,
            plotParams(:atomic_mass).request,
            plotParams(:ionized_mass).request,
        ),
    )

    # Compute the number of bins for the gas quantity
    n_bins = length(edges) - 1

    # Number of bars per bin
    n_bars = 3

    # Compute the dodge argument for `barplot!`
    dodge = repeat(1:n_bars, outer=n_bins)

    # Set the color list
    colors = Makie.wong_colors()[[3,4,1,2]]

    # Compute the axis ticks
    if exp_ticks
        tick_nums = log10.(ustrip.(plot_params.unit, edges))
    else
        tick_nums = ustrip.(plot_params.unit, edges)
    end

    ticks = [string(round((tick_nums[i] + tick_nums[i + 1]) / 2, sigdigits=2)) for i in 1:n_bins]

    n_sims = length(simulation_paths)

    for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        if isone(n_sims)
            base_filename = "fractions_vs_$(quantity)_barplot"
        else
            base_filename = "$(sim_name)_fractions_vs_$(quantity)_barplot"
        end

        plotSnapshot(
            [simulation_path],
            request,
            [barplot!];
            pf_kwargs=[(; dodge, color=colors[dodge])],
            # `plotSnapshot` configuration
            output_path,
            base_filename,
            output_format=".png",
            show_progress=true,
            # Data manipulation options
            slice,
            filter_function,
            da_functions=[daGasFractions],
            da_args=[(quantity, edges)],
            da_kwargs=[(; filter_function=filterByELSFR)],
            post_processing=ppBarPlotLabels,
            pp_args=(include_stars,),
            pp_kwargs=(; colors),
            transform_box=true,
            translation,
            rotation,
            smooth=0,
            x_unit=Unitful.NoUnits,
            y_unit=Unitful.NoUnits,
            x_exp_factor=0,
            y_exp_factor=0,
            x_trim=(-Inf, Inf),
            y_trim=(-Inf, Inf),
            x_edges=false,
            y_edges=false,
            x_func=identity,
            y_func=identity,
            # Axes options
            xaxis_label=L"\mathrm{Fraction} \,\, [%]",
            yaxis_label=isnothing(axis_label) ? plot_params.axis_label : axis_label,
            xaxis_var_name="",
            yaxis_var_name=plot_params.var_name,
            xaxis_scale_func=identity,
            yaxis_scale_func=identity,
            # Plotting and animation options
            save_figures=true,
            backup_results=false,
            theme=merge(
                theme,
                Theme(
                    size=(850, 850),
                    Legend=(nbanks=1,),
                    Axis=(
                        limits=(nothing, 105, nothing, nothing),
                        xticks=([0, 50, 100], [L"0.0", L"50", L"100"]),
                        yticks=(1:n_bins, ticks),
                    ),
                    BarPlot=(
                        flip_labels_at=10,
                        label_formatter=barPlotLabelFormater,
                        label_size=include_stars ? 25 : 35,
                    ),
                ),
            ),
            sim_labels=nothing,
            title="",
            colorbar=false,
            # Animation options
            animation=false,
            animation_filename="animation.mp4",
            framerate=10,
        )

    end

    return nothing

end

"""
    timeSeries(
        simulation_paths::Vector{String},
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a time series.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`              -> Stellar mass.
      + `:gas_mass`                  -> Gas mass.
      + `:hydrogen_mass`             -> Hydrogen mass.
      + `:dm_mass`                   -> Dark matter mass.
      + `:bh_mass`                   -> Black hole mass.
      + `:molecular_mass`            -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`               -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`              -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`              -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`          -> Stellar gas mass (according to out SF model).
      + `:stellar_number`            -> Number of stellar particles.
      + `:gas_number`                -> Number of gas cells.
      + `:dm_number`                 -> Number of dark matter particles.
      + `:bh_number`                 -> Number of black hole particles.
      + `:molecular_fraction`        -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`     -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`           -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`          -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`          -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`-> Fraction of molecular hydrogen in the neutral gas.
      + `:ionized_neutral_fraction`  -> Fraction of ionized gas to neutral gas.
      + `:gas_mass_density`          -> Mean gas mass density.
      + `:stellar_gas_fraction`      -> Stellar gas fraction (according to out SF model).
      + `:stellar_area_density`      -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`          -> Gas mass surface density, for a radius of `DISK_R`.
      + `:molecular_area_density`    -> Molecular mass surface density, for a radius of `DISK_R`.
      + `:br_molecular_area_density` -> Molecular mass surface density, for a radius of `DISK_R`, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_area_density`       -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`      -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`      -> Neutral mass surface density, for a radius of `DISK_R`.
      + `:sfr_area_density`          -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:gas_td`                    -> The mean total gas depletion time.
      + `:molecular_td`              -> The mean molecular hydrogen (``\\mathrm{H_2}``) depletion time.
      + `:br_molecular_td`           -> The mean molecular hydrogen (``\\mathrm{H_2}``) depletion time, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_td`                 -> The mean atomic hydrogen (``\\mathrm{HI}``) depletion time.
      + `:ionized_td`                -> The mean ionized hydrogen (``\\mathrm{HII}``) depletion time.
      + `:neutral_td`                -> The mean neutral hydrogen (``\\mathrm{HI + H_2}``) depletion time.
      + `:gas_metallicity`           -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`       -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`           -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`       -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`       -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`           -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`            -> Norm of the dark matter specific angular momentum.
      + `:sfr`                       -> The star formation rate.
      + `:ssfr`                      -> The specific star formation rate.
      + `:observational_sfr`         -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`        -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:stellar_eff`               -> The mean star formation efficiency per free-fall time for the gas that has turn into stars.
      + `:gas_eff`                   -> The mean star formation efficiency per free-fall time for the gas.
      + `:molecular_eff`             -> The mean star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas.
      + `:br_molecular_eff`          -> The mean star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_eff`                -> The mean star formation efficiency per free-fall time for the atomic hydrogen (``\\mathrm{HI}``) gas.
      + `:ionized_eff`               -> The mean star formation efficiency per free-fall time for the ionized hydrogen (``\\mathrm{HII}``) gas.
      + `:neutral_eff`               -> The mean star formation efficiency per free-fall time for the neutral hydrogen (``\\mathrm{HI + H_2}``) gas.
      + `:scale_factor`              -> Scale factor.
      + `:redshift`                  -> Redshift.
      + `:physical_time`             -> Physical time since the Big Bang.
      + `:lookback_time`             -> Physical time left to reach the last snapshot.
      + `:ode_gas_it`                -> Integration time.
      + `:ode_gas_accu_it`           -> Accumulated integration time.
      + `:ode_gas_tau_s`             -> Star formation time scale, ``\\tau_\\mathrm{S}``.
      + `:ode_gas_eta_d`             -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``.
      + `:ode_gas_eta_i`             -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``.
      + `:ode_gas_r`                 -> Mass recycling parameter, ``R``.
      + `:ode_gas_cold_mf`           -> Cold gas mass fraction.
      + `:ode_stellar_it`            -> Integration time, for the gas that form the stars.
      + `:ode_stellar_accu_it`       -> Accumulated integration time, for the gas that form the stars.
      + `:ode_stellar_tau_s`         -> Star formation time scale, ``\\tau_\\mathrm{S}``, for the gas that form the stars.
      + `:ode_stellar_eta_d`         -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``, for the gas that form the stars.
      + `:ode_stellar_eta_i`         -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``, for the gas that form the stars.
      + `:ode_stellar_r`             -> Mass recycling parameter, ``R``, for the gas that form the stars.
      + `:ode_stellar_cold_mf`       -> Cold gas mass fraction, for the gas that form the stars.
      + `:ode_stellar_gas_rho`       -> Gas mass density, for the gas that form the stars.
      + `:ode_stellar_gas_Z`         -> Gas metallicity, for the gas that form the stars (solar units).
      + `:ode_stellar_gas_mass`      -> Cell mass, for the gas that form the stars.
      + `:ode_stellar_gas_sfr`       -> SFR associated to the gas particles/cells within the code, for the gas that form the stars.
      + `:ode_stellar_gas_P`         -> Gas pressure, for the gas that form the stars.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`              -> Stellar mass.
      + `:gas_mass`                  -> Gas mass.
      + `:hydrogen_mass`             -> Hydrogen mass.
      + `:dm_mass`                   -> Dark matter mass.
      + `:bh_mass`                   -> Black hole mass.
      + `:molecular_mass`            -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`         -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`               -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`              -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`              -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`          -> Stellar gas mass (according to out SF model).
      + `:stellar_number`            -> Number of stellar particles.
      + `:gas_number`                -> Number of gas cells.
      + `:dm_number`                 -> Number of dark matter particles.
      + `:bh_number`                 -> Number of black hole particles.
      + `:molecular_fraction`        -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`     -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`           -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`          -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`          -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`-> Fraction of molecular hydrogen in the neutral gas.
      + `:ionized_neutral_fraction`  -> Fraction of ionized gas to neutral gas.
      + `:gas_mass_density`          -> Mean gas mass density.
      + `:stellar_gas_fraction`      -> Stellar gas fraction (according to out SF model).
      + `:stellar_area_density`      -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`          -> Gas mass surface density, for a radius of `DISK_R`.
      + `:molecular_area_density`    -> Molecular mass surface density, for a radius of `DISK_R`.
      + `:br_molecular_area_density` -> Molecular mass surface density, for a radius of `DISK_R`, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_area_density`       -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`      -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`      -> Neutral mass surface density, for a radius of `DISK_R`.
      + `:sfr_area_density`          -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:gas_td`                    -> The mean total gas depletion time.
      + `:molecular_td`              -> The mean molecular hydrogen (``\\mathrm{H_2}``) depletion time.
      + `:br_molecular_td`           -> The mean molecular hydrogen (``\\mathrm{H_2}``) depletion time, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_td`                 -> The mean atomic hydrogen (``\\mathrm{HI}``) depletion time.
      + `:ionized_td`                -> The mean ionized hydrogen (``\\mathrm{HII}``) depletion time.
      + `:neutral_td`                -> The mean neutral hydrogen (``\\mathrm{HI + H_2}``) depletion time.
      + `:gas_metallicity`           -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`       -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`           -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`       -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_specific_am`       -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`           -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`            -> Norm of the dark matter specific angular momentum.
      + `:sfr`                       -> The star formation rate.
      + `:ssfr`                      -> The specific star formation rate.
      + `:observational_sfr`         -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`        -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:stellar_eff`               -> The mean star formation efficiency per free-fall time for the gas that has turn into stars.
      + `:gas_eff`                   -> The mean star formation efficiency per free-fall time for the gas.
      + `:molecular_eff`             -> The mean star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas.
      + `:br_molecular_eff`          -> The mean star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_eff`                -> The mean star formation efficiency per free-fall time for the atomic hydrogen (``\\mathrm{HI}``) gas.
      + `:ionized_eff`               -> The mean star formation efficiency per free-fall time for the ionized hydrogen (``\\mathrm{HII}``) gas.
      + `:neutral_eff`               -> The mean star formation efficiency per free-fall time for the neutral hydrogen (``\\mathrm{HI + H_2}``) gas.
      + `:scale_factor`              -> Scale factor.
      + `:redshift`                  -> Redshift.
      + `:physical_time`             -> Physical time since the Big Bang.
      + `:lookback_time`             -> Physical time left to reach the last snapshot.
      + `:ode_gas_it`                -> Integration time.
      + `:ode_gas_accu_it`           -> Accumulated integration time.
      + `:ode_gas_tau_s`             -> Star formation time scale, ``\\tau_\\mathrm{S}``.
      + `:ode_gas_eta_d`             -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``.
      + `:ode_gas_eta_i`             -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``.
      + `:ode_gas_r`                 -> Mass recycling parameter, ``R``.
      + `:ode_gas_cold_mf`           -> Cold gas mass fraction.
      + `:ode_stellar_it`            -> Integration time, for the gas that form the stars.
      + `:ode_stellar_accu_it`       -> Accumulated integration time, for the gas that form the stars.
      + `:ode_stellar_tau_s`         -> Star formation time scale, ``\\tau_\\mathrm{S}``, for the gas that form the stars.
      + `:ode_stellar_eta_d`         -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``, for the gas that form the stars.
      + `:ode_stellar_eta_i`         -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``, for the gas that form the stars.
      + `:ode_stellar_r`             -> Mass recycling parameter, ``R``, for the gas that form the stars.
      + `:ode_stellar_cold_mf`       -> Cold gas mass fraction, for the gas that form the stars.
      + `:ode_stellar_gas_rho`       -> Gas mass density, for the gas that form the stars.
      + `:ode_stellar_gas_Z`         -> Gas metallicity, for the gas that form the stars (solar units).
      + `:ode_stellar_gas_mass`      -> Cell mass, for the gas that form the stars.
      + `:ode_stellar_gas_sfr`       -> SFR associated to the gas particles/cells within the code, for the gas that form the stars.
      + `:ode_stellar_gas_P`         -> Gas pressure, for the gas that form the stars.
  - `y_log::Bool=true`: If the y axis is will have a log10 scale. Only works if `fraction` = false.
  - `cumulative::Bool=false`: If the `y_quantity` will be accumulated or not.
  - `fraction::Bool=false`: If the `y_quantity` will be represented as a fraction of the last value. If `cumulative` = true, this will apply to the accumulated values.
  - `slice::IndexType=(:)`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `extra_filter::Function=filterNothing`: Filter function that will be applied after the one given by `filter_mode`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `extra_filter` filter function.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `backup_results::Bool=false`: If the values to be plotted will be backup in a [JLD2](https://github.com/JuliaIO/JLD2.jl) file.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function timeSeries(
    simulation_paths::Vector{String},
    x_quantity::Symbol,
    y_quantity::Symbol;
    y_log::Bool=true,
    cumulative::Bool=false,
    fraction::Bool=false,
    slice::IndexType=(:),
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    extra_filter::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    backup_results::Bool=false,
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    y_var_name = y_plot_params.var_name

    if cumulative
        y_var_name = "Accumulated $(y_var_name)"
    end

    if fraction
        y_var_name = "Fractional $(y_var_name)"
        filename = "$(y_quantity)_vs_$(x_quantity)_fractional"
    else
        filename = "$(y_quantity)_vs_$(x_quantity)"
    end

    if fraction || !y_log
        yaxis_scale_func = identity
    else
        yaxis_scale_func = log10
    end

    plotTimeSeries(
        simulation_paths,
        [lines!];
        pf_kwargs=[(;)],
        # `plotTimeSeries` configuration
        output_path,
        filename,
        output_format=".png",
        show_progress=true,
        # Data manipulation options
        slice,
        da_functions=[daEvolution],
        da_args=[(x_quantity, y_quantity)],
        da_kwargs=[
            (;
                filter_mode,
                extra_filter,
                ff_request,
                smooth=0,
                cumulative,
                fraction,
                scaling=identity,
            )
        ],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        x_unit=x_plot_params.unit,
        y_unit=fraction ? Unitful.NoUnits : y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=fraction ? 0 : y_plot_params.exp_factor,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func,
        # Plotting options
        save_figures=!backup_results,
        backup_results,
        theme,
        sim_labels,
        title="",
    )

    return nothing

end

"""
    gasEvolution(
        simulation_paths::Vector{String};
        <keyword arguments>
    )::Nothing

Plot a time series of the gas components. Either their masses or their fractions.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `fractions::Bool=true`: If the fractions (default), or the masses, will be plotted.
  - `slice::IndexType=(:)`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `extra_filter::Function=filterNothing`: Filter function that will be applied after the one given by `filter_mode`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `extra_filter` filter function.
  - `filename::Union{String,Nothing}=nothing`: Name for the output file. If left as `nothing`, the filename will be chosen automaticaly.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function gasEvolution(
    simulation_paths::Vector{String};
    fractions::Bool=false,
    slice::IndexType=(:),
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    extra_filter::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    filename::Union{String,Nothing}=nothing,
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(:physical_time)

    if fractions
        quantities = [:ionized_fraction, :atomic_fraction, :molecular_fraction]
        sim_labels = ["Ionized fraction", "Atomic fraction", "Molecular fraction"]
        y_plot_params = plotParams(:generic_fraction)
    else
        quantities = [:stellar_mass, :hydrogen_mass, :ionized_mass, :atomic_mass, :molecular_mass]
        sim_labels = [
            "Stellar mass",
            "Hydrogen mass",
            "Ionized mass",
            "Atomic mass",
            "Molecular mass",
        ]
        y_plot_params = plotParams(:generic_mass)
    end

    for simulation_path in simulation_paths

        if isnothing(filename)
            if fractions
                filename = "gas_fractions_vs_physical_time_$(basename(simulation_path))"
            else
                filename = "gas_masses_vs_physical_time_$(basename(simulation_path))"
            end
        end

        plotTimeSeries(
            fill(simulation_path, length(quantities)),
            [lines!];
            pf_kwargs=[(;)],
            # `plotTimeSeries` configuration
            output_path,
            filename,
            output_format=".png",
            show_progress=true,
            # Data manipulation options
            slice,
            da_functions=[daEvolution],
            da_args=[(:physical_time, quantity) for quantity in quantities],
            da_kwargs=[(; filter_mode, extra_filter, ff_request, smooth=0, scaling=identity)],
            post_processing=getNothing,
            pp_args=(),
            pp_kwargs=(;),
            x_unit=x_plot_params.unit,
            y_unit=y_plot_params.unit,
            x_exp_factor=x_plot_params.exp_factor,
            y_exp_factor=y_plot_params.exp_factor,
            x_trim=(-Inf, Inf),
            y_trim=(-Inf, Inf),
            x_edges=false,
            y_edges=false,
            x_func=identity,
            y_func=identity,
            # Axes options
            xaxis_label=x_plot_params.axis_label,
            yaxis_label=y_plot_params.axis_label,
            xaxis_var_name=x_plot_params.var_name,
            yaxis_var_name=y_plot_params.var_name,
            xaxis_scale_func=identity,
            yaxis_scale_func=log10,
            # Plotting options
            save_figures=true,
            backup_results=false,
            theme,
            sim_labels,
            title="",
        )

    end

    return nothing

end

"""
    virialAccretionEvolution(
        simulation_paths::Vector{String};
        <keyword arguments>
    )::Nothing

Plot a time series of the accreted mass into the virial radius.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType=(:)`: Slice of the simulations, i.e. which snapshots will be plotted. It can be vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.
  - `tracers::Bool=false`: If tracers will be use to compute the mass accretion.
  - `smooth::Int=0`: The time series will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `output_path::String="./"`: Path to the output folder.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function virialAccretionEvolution(
    simulation_paths::Vector{String};
    slice::IndexType=(:),
    halo_idx::Int=1,
    tracers::Bool=false,
    smooth::Int=0,
    output_path::String="./",
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(:physical_time)
    y_plot_params = plotParams(:mass_accretion)

    if tracers
        filename="virial-mass-accretion_with_tracers"
    else
        filename="virial-mass-change_evolution"
    end

    plotTimeSeries(
        simulation_paths,
        [lines!];
        pf_kwargs=[(;)],
        # `plotTimeSeries` configuration
        output_path,
        filename,
        output_format=".png",
        show_progress=true,
        # Data manipulation options
        slice,
        da_functions=[daVirialAccretion],
        da_args=[()],
        da_kwargs=[(; filter_mode=:halo, halo_idx, tracers, smooth)],
        post_processing=ppHorizontalFlags!,
        pp_args=([0.0],),
        pp_kwargs=(; colors=[:gray65], line_styles=[nothing]),
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=y_plot_params.exp_factor,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=tracers ? y_plot_params.var_name : "Net mass change",
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        # Plotting options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
    )

    return nothing

end

"""
    discAccretionEvolution(
        simulation_paths::Vector{String};
        <keyword arguments>
    )::Nothing

Plot a time series of the accreted mass into the disc.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType=(:)`: Slice of the simulations, i.e. which snapshots will be plotted. It can be vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `max_r::Unitful.Length=DISK_R`: Radius of the cylinder.
  - `max_z::Unitful.Length=5.0u"kpc"`: Half height of the cylinder.
  - `smooth::Int=0`: The time series will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `output_path::String="./"`: Path to the output folder.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=nothing`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function discAccretionEvolution(
    simulation_paths::Vector{String};
    slice::IndexType=(:),
    max_r::Unitful.Length=DISK_R,
    max_z::Unitful.Length=5.0u"kpc",
    smooth::Int=0,
    output_path::String="./",
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(:physical_time)
    y_plot_params = plotParams(:mass_accretion)

    plotTimeSeries(
        simulation_paths,
        [lines!];
        pf_kwargs=[(;)],
        # `plotTimeSeries` configuration
        output_path,
        filename="disc_mass_accretion_with_tracers",
        output_format=".png",
        show_progress=true,
        # Data manipulation options
        slice,
        da_functions=[daDiscAccretion],
        da_args=[()],
        da_kwargs=[(; filter_mode=:halo, max_r, max_z, smooth)],
        post_processing=ppHorizontalFlags!,
        pp_args=([0.0],),
        pp_kwargs=(; colors=[:gray65], line_styles=[nothing]),
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=y_plot_params.exp_factor,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        # Plotting options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
    )

    return nothing

end

"""
    rotationCurve(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the galaxy rotation curve of a set of simulations.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `radius::Unitful.Length=DISK_R`: Maximum radial distance for the rotation curve.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function rotationCurve(
    simulation_paths::Vector{String},
    slice::IndexType;
    radius::Unitful.Length=DISK_R,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(:stellar_radial_distance)
    y_plot_params = plotParams(:stellar_vcirc)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(x_plot_params.request, y_plot_params.request, Dict(:stars => ["GAGE"])),
    )

    plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="rotation_curve",
        output_format=".png",
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daRotationCurve],
        da_args=[(radius,)],
        da_kwargs=[(;)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth=round(Int64, 5 * ustrip(u"kpc", radius)),
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=y_plot_params.exp_factor,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    densityProfile(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a density profile.

!!! note

    This method plots one quantity for several simulations in one figure.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantities::Vector{Symbol}`: Quantities for the y axis. The options are:

      + `:stellar_mass`       -> Stellar mass.
      + `:gas_mass`           -> Gas mass.
      + `:hydrogen_mass`      -> Hydrogen mass.
      + `:dm_mass`            -> Dark matter mass.
      + `:bh_mass`            -> Black hole mass.
      + `:molecular_mass`     -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`  -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`        -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`       -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`       -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`   -> Stellar gas mass (according to out SF model).
      + `:sfr`                -> The star formation rate.
      + `:ssfr`               -> The specific star formation rate.
      + `:observational_sfr`  -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr` -> The specific star formation rate of the last `AGE_RESOLUTION`.
  - `cumulative::Bool=false`: If the profile will be accumulated or not.
  - `yscale::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `radius::Unitful.Length=DISK_R`: Radius of the profile.
  - `n_bins::Int=100`: Number of bins.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function densityProfile(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol;
    cumulative::Bool=false,
    yscale::Function=identity,
    radius::Unitful.Length=DISK_R,
    n_bins::Int=100,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(quantity)
    filter_function, translation, rotation, request = selectFilter(filter_mode, plot_params.request)

    if quantity == :stellar_mass

        yaxis_var_name = L"\Sigma_\star"

    elseif quantity == :dm_mass

        yaxis_var_name = L"\Sigma_\mathrm{DM}"

    elseif quantity == :bh_mass

        yaxis_var_name = L"\Sigma_\mathrm{BH}"

    elseif quantity == :gas_mass

        yaxis_var_name = L"\Sigma_\mathrm{gas}"

    elseif quantity == :hydrogen_mass

        yaxis_var_name = L"\Sigma_\mathrm{H}"

    elseif quantity == :molecular_mass

        yaxis_var_name = L"\Sigma_\mathrm{H_2}"

    elseif quantity == :br_molecular_mass

        yaxis_var_name = L"\Sigma_\mathrm{H_2}^\mathrm{BR}"

    elseif quantity == :atomic_mass

        yaxis_var_name = L"\Sigma_\mathrm{HI}"

    elseif quantity == :ionized_mass

        yaxis_var_name = L"\Sigma_\mathrm{HII}"

    elseif quantity == :neutral_mass

        yaxis_var_name = L"\Sigma_\mathrm{HI + H_2}"

    elseif quantity ∈ [:sfr, :observational_sfr]

        yaxis_var_name = L"\Sigma_\mathrm{SFR}"

    elseif quantity ∈ [:ssfr, :observational_ssfr]

        yaxis_var_name = L"\Sigma_\mathrm{sSFR}"

    else

        throw(ArgumentError("densityProfile: I don't know the quantity :$(quantity)"))

    end

    grid = CircularGrid(radius, n_bins)

    # Draw the figures with CairoMakie
    plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="$(quantity)_density_profile",
        output_format=".pdf",
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daProfile],
        da_args=[(quantity, grid)],
        da_kwargs=[(; flat=true, total=true, cumulative, density=true)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=u"kpc",
        y_unit=plot_params.unit * u"kpc^-2",
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label="auto_label",
        yaxis_label=plot_params.axis_label,
        xaxis_var_name=L"r",
        yaxis_var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=yscale,
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    densityProfile(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantities::Vector{Symbol};
        <keyword arguments>
    )::Nothing

Plot a density profile.

!!! note

    This method plots several quantities for one simulations in one figure.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantities::Vector{Symbol}`: Quantities for the y axis. The options are:

      + `:stellar_mass`       -> Stellar mass.
      + `:gas_mass`           -> Gas mass.
      + `:hydrogen_mass`      -> Hydrogen mass.
      + `:dm_mass`            -> Dark matter mass.
      + `:bh_mass`            -> Black hole mass.
      + `:molecular_mass`     -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`  -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`        -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`       -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`       -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`   -> Stellar gas mass (according to out SF model).
      + `:sfr`                -> The star formation rate.
      + `:ssfr`               -> The specific star formation rate.
      + `:observational_sfr`  -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr` -> The specific star formation rate of the last `AGE_RESOLUTION`.
  - `cumulative::Bool=false`: If the profile will be accumulated or not.
  - `yscale::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `radius::Unitful.Length=DISK_R`: Radius of the profile.
  - `n_bins::Int=100`: Number of bins.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=string.(quantities)`: Labels for the plot legend, one per quantity. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function densityProfile(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantities::Vector{Symbol};
    cumulative::Bool=false,
    yscale::Function=identity,
    radius::Unitful.Length=DISK_R,
    n_bins::Int=100,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{<:AbstractString},Nothing}=string.(quantities),
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(:generic_area_density)
    filter_function, translation, rotation, request = selectFilter(filter_mode, plot_params.request)

    grid = CircularGrid(radius, n_bins)

    for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        # Draw the figures with CairoMakie
        plotSnapshot(
            fill(simulation_path, length(quantities)),
            request,
            [lines!];
            pf_kwargs=[(;)],
            # `plotSnapshot` configuration
            output_path,
            base_filename="$(sim_name)_density_profiles",
            output_format=".png",
            show_progress=true,
            # Data manipulation options
            slice,
            filter_function,
            da_functions=[daProfile],
            da_args=[(quantity, grid) for quantity in quantities],
            da_kwargs=[(; flat=true, total=true, cumulative, density=true)],
            post_processing=getNothing,
            pp_args=(),
            pp_kwargs=(;),
            transform_box=true,
            translation,
            rotation,
            smooth=0,
            x_unit=u"kpc",
            y_unit=plot_params.unit,
            x_exp_factor=0,
            y_exp_factor=0,
            x_trim=(-Inf, Inf),
            y_trim=(-Inf, Inf),
            x_edges=false,
            y_edges=false,
            x_func=identity,
            y_func=identity,
            # Axes options
            xaxis_label="auto_label",
            yaxis_label=plot_params.axis_label,
            xaxis_var_name=L"r",
            yaxis_var_name=plot_params.var_name,
            xaxis_scale_func=identity,
            yaxis_scale_func=yscale,
            # Plotting and animation options
            save_figures=true,
            backup_results=false,
            theme,
            sim_labels,
            title="",
            colorbar=false,
            # Animation options
            animation=false,
            animation_filename="animation.mp4",
            framerate=10,
        )

    end

    return nothing

end

"""
    massProfile(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantities::Vector{Symbol};
        <keyword arguments>
    )::Nothing

Plot a mass profile.

!!! note

    This method plots several quantities for one simulations in one figure.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantities::Vector{Symbol}`: Quantities for the y axis. The options are:

      + `:stellar_mass`      -> Stellar mass.
      + `:gas_mass`          -> Gas mass.
      + `:hydrogen_mass`     -> Hydrogen mass.
      + `:dm_mass`           -> Dark matter mass.
      + `:bh_mass`           -> Black hole mass.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`      -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`  -> Stellar gas mass (according to out SF model).
  - `cumulative::Bool=false`: If the profile will be accumulated or not.
  - `yscale::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `radius::Unitful.Length=DISK_R`: Radius of the profile.
  - `n_bins::Int=100`: Number of bins.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=string.(quantities)`: Labels for the plot legend, one per quantity. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function massProfile(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantities::Vector{Symbol};
    cumulative::Bool=false,
    yscale::Function=identity,
    radius::Unitful.Length=DISK_R,
    n_bins::Int=100,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{<:AbstractString},Nothing}=string.(quantities),
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(:generic_mass)
    filter_function, translation, rotation, request = selectFilter(filter_mode, plot_params.request)

    grid = CircularGrid(radius, n_bins)

    n_sims = length(simulation_paths)

    for simulation_path in simulation_paths

        # Get the simulation name as a string
        sim_name = basename(simulation_path)

        if isone(n_sims)
            if cumulative
                base_filename = "mass_profiles_cumulative"
            else
                base_filename = "mass_profiles"
            end
        else
            if cumulative
                base_filename = "$(sim_name)_mass_profiles_cumulative"
            else
                base_filename = "$(sim_name)_mass_profiles"
            end
        end

        # Draw the figures with CairoMakie
        plotSnapshot(
            fill(simulation_path, length(quantities)),
            request,
            [lines!];
            pf_kwargs=[(;)],
            # `plotSnapshot` configuration
            output_path,
            base_filename,
            output_format=".png",
            show_progress=true,
            # Data manipulation options
            slice,
            filter_function,
            da_functions=[daProfile],
            da_args=[(quantity, grid) for quantity in quantities],
            da_kwargs=[(; flat=true, total=true, cumulative, density=false)],
            post_processing=getNothing,
            pp_args=(),
            pp_kwargs=(;),
            transform_box=true,
            translation,
            rotation,
            smooth=0,
            x_unit=u"kpc",
            y_unit=plot_params.unit,
            x_exp_factor=0,
            y_exp_factor=0,
            x_trim=(-Inf, Inf),
            y_trim=(-Inf, Inf),
            x_edges=false,
            y_edges=false,
            x_func=identity,
            y_func=identity,
            # Axes options
            xaxis_label="auto_label",
            yaxis_label=plot_params.axis_label,
            xaxis_var_name=L"r",
            yaxis_var_name=plot_params.var_name,
            xaxis_scale_func=identity,
            yaxis_scale_func=yscale,
            # Plotting and animation options
            save_figures=true,
            backup_results=false,
            theme,
            sim_labels,
            title="",
            colorbar=false,
            # Animation options
            animation=false,
            animation_filename="animation.mp4",
            framerate=10,
        )
        end

    return nothing

end

"""
    velocityProfile(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a velocity profile.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `component::Symbol`: Which component will be calculated. The options are:

      + `:stellar_vradial`     -> Stellar radial speed (``v_r``).
      + `:stellar_vtangential` -> Stellar tangential speed (``v_\\theta``).
      + `:stellar_vzstar`      -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
  - `yscale::Function=identity`: Scaling function for the y axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function velocityProfile(
    simulation_paths::Vector{String},
    slice::IndexType,
    component::Symbol;
    yscale::Function=identity,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(component)
    filter_function, translation, rotation, request = selectFilter(filter_mode, plot_params.request)

    grid = CircularGrid(DISK_R, 25)

    # Draw the figures with CairoMakie
    plotSnapshot(
        simulation_paths,
        request,
        [scatterlines!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="$(component)_profile",
        output_format=".png",
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daProfile],
        da_args=[(component, grid)],
        da_kwargs=[(; flat=true, total=false, cumulative=false, density=false)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=u"kpc",
        y_unit=plot_params.unit,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label="auto_label",
        yaxis_label=plot_params.axis_label,
        xaxis_var_name=L"r",
        yaxis_var_name=plot_params.var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=yscale,
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    stellarHistory(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot the evolution of a given stellar `quantity` using the stellar ages at a given instant in time.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Snapshot at which the stellar ages will be read. If set to several snapshots, one plot per snapshot will be done. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: Quantity for the y axis. The options are:

      + `:sfr`                 -> The star formation rate.
      + `:ssfr`                -> The specific star formation rate.
      + `:stellar_mass`        -> Stellar mass.
      + `:stellar_metallicity` -> Mass fraction of all elements above He in the stars (solar units).
  - `y_log::Bool=true`: If the y axis is will have a log10 scale.
  - `n_bins::Int=20`: Number of bins (time intervals).
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `backup_results::Bool=false`: If the values to be plotted will be backup in a [JLD2](https://github.com/JuliaIO/JLD2.jl) file.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function stellarHistory(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol;
    y_log::Bool=true,
    n_bins::Int=20,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    backup_results::Bool=false,
    theme::Attributes=Theme(),
)::Nothing

    x_plot_params = plotParams(:physical_time)
    y_plot_params = plotParams(quantity)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        y_plot_params.request,
    )

    # Draw the figures with CairoMakie
    plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="$(quantity)_history",
        output_format=".png",
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daStellarHistory],
        da_args=[()],
        da_kwargs=[(; quantity, n_bins)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=x_plot_params.exp_factor,
        y_exp_factor=y_plot_params.exp_factor,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=y_log ? log10 : identity,
        # Plotting and animation options
        save_figures=!backup_results,
        backup_results,
        theme,
        sim_labels,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    lineHistogram(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol,
        type::Symbol,
        range::NTuple{2,<:Number};
        <keyword arguments>
    )::Nothing

Plot a histogram of `quantity`.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: The possibilities are:

      + `:stellar_mass`                -> Stellar mass.
      + `:gas_mass`                    -> Gas mass.
      + `:hydrogen_mass`               -> Hydrogen mass.
      + `:dm_mass`                     -> Dark matter mass.
      + `:bh_mass`                     -> Black hole mass.
      + `:molecular_mass`              -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass`           -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`                 -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:ionized_mass`                -> Ionized hydrogen (``\\mathrm{HII}``) mass.
      + `:neutral_mass`                -> Neutral hydrogen (``\\mathrm{HI + H_2}``) mass.
      + `:stellar_gas_mass`            -> Stellar gas mass (according to out SF model).
      + `:molecular_fraction`          -> Gas mass fraction of molecular hydrogen.
      + `:br_molecular_fraction`       -> Gas mass fraction of molecular hydrogen, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_fraction`             -> Gas mass fraction of atomic hydrogen.
      + `:ionized_fraction`            -> Gas mass fraction of ionized hydrogen.
      + `:neutral_fraction`            -> Gas mass fraction of neutral hydrogen.
      + `:molecular_neutral_fraction`  -> Fraction of molecular hydrogen in the neutral gas.
      + `:ionized_neutral_fraction`    -> Fraction of ionized gas to neutral gas.
      + `:stellar_gas_fraction`        -> Stellar gas fraction (according to out SF model).
      + `:mol_eq_quotient`             -> Equilibrium quotient for the molecular fraction equation of the SF model.
      + `:ion_eq_quotient`             -> Equilibrium quotient for the ionized fraction equation of the SF model.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:gas_td`                      -> Total gas depletion time.
      + `:molecular_td`                -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time.
      + `:br_molecular_td`             -> Molecular hydrogen (``\\mathrm{H_2}``) depletion time, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_td`                   -> Atomic hydrogen (``\\mathrm{HI}``) depletion time.
      + `:ionized_td`                  -> Ionized hydrogen (``\\mathrm{HII}``) depletion time.
      + `:neutral_td`                  -> Neutral hydrogen (``\\mathrm{HI + H_2}``) depletion time.
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:stellar_metallicity`         -> Mass fraction of all elements above He in the stars (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:X_stellar_abundance`         -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:stellar_radial_distance`     -> Distance of every stellar particle to the origin.
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:dm_radial_distance`          -> Distance of every dark matter particle to the origin.
      + `:stellar_xy_distance`         -> Projected distance of every stellar particle to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:dm_xy_distance`              -> Projected distance of every dark matter particle to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
      + `:stellar_circularity`         -> Stellar circularity.
      + `:stellar_vcirc`               -> Stellar circular velocity.
      + `:stellar_vradial`             -> Stellar radial speed.
      + `:stellar_vtangential`         -> Stellar tangential speed.
      + `:stellar_vzstar`              -> Stellar speed in the z direction, computed as ``v_z \\, \\mathrm{sign}(z)``.
      + `:stellar_age`                 -> Stellar age.
      + `:sfr`                         -> The star formation rate.
      + `:ssfr`                        -> The specific star formation rate.
      + `:observational_sfr`           -> The star formation rate of the last `AGE_RESOLUTION`.
      + `:observational_ssfr`          -> The specific star formation rate of the last `AGE_RESOLUTION`.
      + `:stellar_eff`                 -> The star formation efficiency per free-fall time for the gas that has turn into stars.
      + `:gas_eff`                     -> The star formation efficiency per free-fall time for the gas.
      + `:molecular_eff`               -> The star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas.
      + `:br_molecular_eff`            -> The star formation efficiency per free-fall time for the molecular hydrogen (``\\mathrm{H_2}``) gas, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_eff`                  -> The star formation efficiency per free-fall time for the atomic hydrogen (``\\mathrm{HI}``) gas.
      + `:ionized_eff`                 -> The star formation efficiency per free-fall time for the ionized hydrogen (``\\mathrm{HII}``) gas.
      + `:neutral_eff`                 -> The star formation efficiency per free-fall time for the neutral hydrogen (``\\mathrm{HI + H_2}``) gas.
      + `:temperature`                 -> Gas temperature, as ``\\log_{10}(T \\, / \\, \\mathrm{K})``.
      + `:pressure`                    -> Gas pressure.
      + `:ode_gas_it`                  -> Integration time.
      + `:ode_gas_accu_it`             -> Accumulated integration time.
      + `:ode_gas_tau_s`               -> Star formation time scale, ``\\tau_\\mathrm{S}``.
      + `:ode_gas_eta_d`               -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``.
      + `:ode_gas_eta_i`               -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``.
      + `:ode_gas_r`                   -> Mass recycling parameter, ``R``.
      + `:ode_gas_cold_mf`             -> Cold gas mass fraction.
      + `:ode_stellar_it`              -> Integration time, for the gas that form the stars.
      + `:ode_stellar_accu_it`         -> Accumulated integration time, for the gas that form the stars.
      + `:ode_stellar_tau_s`           -> Star formation time scale, ``\\tau_\\mathrm{S}``, for the gas that form the stars.
      + `:ode_stellar_eta_d`           -> Photodissociation efficiency, ``\\eta_\\mathrm{diss}``, for the gas that form the stars.
      + `:ode_stellar_eta_i`           -> Photoionization efficiency, ``\\ate_\\mathrm{ion}``, for the gas that form the stars.
      + `:ode_stellar_r`               -> Mass recycling parameter, ``R``, for the gas that form the stars.
      + `:ode_stellar_cold_mf`         -> Cold gas mass fraction, for the gas that form the stars.
      + `:ode_stellar_gas_rho`         -> Gas mass density, for the gas that form the stars.
      + `:ode_stellar_gas_Z`           -> Gas metallicity, for the gas that form the stars (solar units).
      + `:ode_stellar_gas_mass`        -> Cell mass, for the gas that form the stars.
      + `:ode_stellar_gas_sfr`         -> SFR associated to the gas particles/cells within the code, for the gas that form the stars.
      + `:ode_stellar_gas_P`           -> Gas pressure, for the gas that form the stars.
  - `type::Symbol`: Type of cell/particle.
  - ` range::NTuple{2,<:Number}`: Range of values for the histogram.
  - `n_bins::Int=100`: Number of bins.
  - `log::Bool=false`: If the bins will be logarithmic.
  - `norm::Int=0`: Number of count that will be use to normalize the histogram. If left as 0, the histogram will be normalize with the maximum bin count.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `extra_filter::Function=filterNothing`: Filter function that will be applied after the one given by `filter_mode`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `extra_filter` filter function.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function lineHistogram(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol,
    type::Symbol,
    range::NTuple{2,<:Number};
    n_bins::Int=100,
    log::Bool=false,
    norm::Int=0,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    extra_filter::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(quantity)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(plot_params.request, ff_request),
    )

    grid = LinearGrid(range..., n_bins; log)

    plotSnapshot(
        simulation_paths,
        request,
        [lines!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="$(quantity)_line_histogram",
        output_format=".png",
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daLineHistogram],
        da_args=[(quantity, grid, type)],
        da_kwargs=[(; filter_function=extra_filter, norm)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=plot_params.unit,
        y_unit=Unitful.NoUnits,
        x_exp_factor=plot_params.exp_factor,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=plot_params.axis_label,
        yaxis_label="auto_label",
        xaxis_var_name=plot_params.var_name,
        yaxis_var_name=L"\mathrm{Normalized \,\, counts}",
        xaxis_scale_func=log ? log10 : identity,
        yaxis_scale_func=identity,
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme=merge(
            theme,
            Theme(Legend=(nbanks=1, halign=:left, valign=:top, padding=(40, 0, 0, 0)),),
        ),
        sim_labels,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    compareFeldmann2020(
        simulation_paths::Vector{String},
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a time series plus the corresponding experimental results from Feldmann (2020).

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`      -> Stellar mass.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:sfr`               -> The star formation rate of the last `AGE_RESOLUTION`.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`      -> Stellar mass.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:sfr`               -> The star formation rate of the last `AGE_RESOLUTION`.
  - `slice::IndexType=(:)`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `scatter::Bool=false`: If the data will be presented as a line plot with error bands (default), or alternatively, a scatter plot.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)

R. Feldmann (2020). *The link between star formation and gas in nearby galaxies*. Communications Physics **3(226)**. [doi:10.1038/s42005-020-00493-0](https://doi.org/10.1038/s42005-020-00493-0)
"""
function compareFeldmann2020(
    simulation_paths::Vector{String},
    x_quantity::Symbol,
    y_quantity::Symbol;
    slice::IndexType=(:),
    scatter::Bool=false,
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    if x_quantity == :sfr
        x_quantity = :observational_sfr
    end

    if y_quantity == :sfr
        y_quantity = :observational_sfr
    end

    x_plot_params = plotParams(x_quantity)
    y_plot_params = plotParams(y_quantity)

    if x_quantity == :sfr
        x_quantity = :observational_sfr
    end

    plotTimeSeries(
        simulation_paths,
        [scatter!];
        pf_kwargs=[(;)],
        # `plotTimeSeries` configuration
        output_path,
        filename="$(y_quantity)_vs_$(x_quantity)_with_Feldmann2020",
        output_format=".png",
        show_progress=true,
        # Data manipulation options
        slice,
        da_functions=[daEvolution],
        da_args=[(x_quantity, y_quantity)],
        da_kwargs=[(; filter_mode, smooth=0, scaling=identity)],
        post_processing=ppFeldmann2020!,
        pp_args=(x_quantity, y_quantity),
        pp_kwargs=(; scatter),
        x_unit=x_plot_params.unit,
        y_unit=y_plot_params.unit,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label=x_plot_params.axis_label,
        yaxis_label=y_plot_params.axis_label,
        xaxis_var_name=x_plot_params.var_name,
        yaxis_var_name=y_plot_params.var_name,
        xaxis_scale_func=log10,
        yaxis_scale_func=log10,
        # Plotting options
        save_figures=true,
        backup_results=false,
        theme=merge(theme, Theme(size=(850, 850), Legend=(nbanks=1,))),
        sim_labels,
        title="",
    )

    return nothing

end

"""
    compareMolla2015(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot a Milky Way profile plus the corresponding experimental results from Mollá et al. (2015).

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_area_density`      -> Stellar area mass density.
      + `:molecular_area_density`    -> Molecular mass surface density.
      + `:br_molecular_area_density` -> Molecular mass surface density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_area_density`       -> Atomic hydrogen area mass density.
      + `:sfr_area_density`          -> Star formation rate area density, for the last `AGE_RESOLUTION`.
      + `:O_stellar_abundance`       -> Stellar abundance of oxygen, as ``12 + \\log_{10}(\\mathrm{O \\, / \\, H})``.
      + `:N_stellar_abundance`       -> Stellar abundance of nitrogen, as ``12 + \\log_{10}(\\mathrm{N \\, / \\, H})``.
      + `:C_stellar_abundance`       -> Stellar abundance of carbon, as ``12 + \\log_{10}(\\mathrm{C \\, / \\, H})``.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function*. Monthly Notices of the Royal Astronomical Society **451(4)**, 3693–3708. [doi:10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)
"""
function compareMolla2015(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol;
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    plot_params = plotParams(quantity)
    request = addRequest(plot_params.request, Dict(:gas => ["VEL "], :stars => ["VEL "]))
    filter_function, translation, rotation, request = selectFilter(filter_mode, plot_params.request)

    # Select the correct grid acording to the available data from M. Mollá et al. (2015)
    if quantity == :stellar_area_density

        grid = CircularGrid(16.5u"kpc", 14; shift=2.5u"kpc")

    elseif quantity ∈ [:molecular_area_density, :br_molecular_area_density, :sfr_area_density]

        grid = CircularGrid(19.5u"kpc", 20; shift=-0.5u"kpc")

    elseif quantity == :atomic_area_density

        grid = CircularGrid(20.5u"kpc", 21; shift=-0.5u"kpc")

    elseif quantity == :O_stellar_abundance

        grid = CircularGrid(18.5u"kpc", 19; shift=-0.5u"kpc")

    elseif quantity == :N_stellar_abundance

        grid = CircularGrid(17.5u"kpc", 18; shift=-0.5u"kpc")

    elseif quantity == :C_stellar_abundance

        grid = CircularGrid(15.5u"kpc", 16; shift=-0.5u"kpc")

    end

    # Draw the figures with CairoMakie
    plotSnapshot(
        simulation_paths,
        request,
        [scatterlines!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="$(quantity)_profile_with_Molla2015",
        output_format=".png",
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daMolla2015],
        da_args=[(grid, quantity)],
        da_kwargs=[(;)],
        post_processing=ppMolla2015!,
        pp_args=(quantity,),
        pp_kwargs=(; color=Makie.wong_colors()[6], linestyle=nothing, error_bars=true),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=u"kpc",
        y_unit=plot_params.unit,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label="auto_label",
        yaxis_label=plot_params.axis_label,
        xaxis_var_name=L"r",
        yaxis_var_name=plot_params.var_name,
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end

"""
    kennicuttSchmidtLaw(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the Kennicutt-Schmidt law.

!!! note

    Only stars younger than [`AGE_RESOLUTION`](@ref) are consider. The star formation surface density is the stellar mass surface density divided by [`AGE_RESOLUTION`](@ref).

!!! note

    This function uses physical units regardless of the [`PHYSICAL_UNITS`](@ref) global setting.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`. All the simulations will be plotted together.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored. All the selected snapshots will be plotted together.
  - `quantity::Symbol=:molecular_mass`: Quantity for the x axis. The options are:

      + `:gas_mass`          -> Total gas mass surface density.
      + `:molecular_mass`    -> Molecular mass surface density. This one can be plotted with the results of Bigiel et al. (2008) and Sun et al. (2023).
      + `:br_molecular_mass` -> Molecular mass surface density, computed using the pressure relation in Blitz et al. (2006). This one can be plotted with the results of Bigiel et al. (2008) and Sun et al. (2023).
      + `:neutral_mass`      -> Neutral mass surface density. This one can be plotted with the results of Bigiel et al. (2008), and Kennicutt (1998).
  - `gas_type::Symbol=:cells`: If the gas surface density will be calculated assuming the gas is in :particles or in Voronoi :cells.
  - `reduce_grid::Symbol=:square`: Grid for the density projection. The options are:

      + `:square`    -> The gas and stellar distributions will be projected into a regular cubic grid first and then into a flat square one, to emulate the way the surface densities are measured in observations.
      + `:circular` -> The gas and stellar distributions will be projected into a regular cubic grid first, then into a flat square one, and finally into a flat circular grid, formed by a series of concentric rings. This emulates the traditonal way the Kennicutt-Schmidt law is measured in simulations.
  - `grid_size::Unitful.Length=BOX_L`: Physical side length of the cubic and square grids, and diameter of the circular grid (if `reduce_grid` = :circular). As a reference, Bigiel et al. (2008) uses measurements up to the optical radius r25 (where the B-band magnitude drops below 25 mag arcsec^−2). This limits which cells/particles will be consider.
  - `bin_size::Unitful.Length=BIGIEL_PX_SIZE`: Target bin size for the grids. If `reduce_grid` = :square, it is the physical side length of the pixels in the final square grid. If `reduce_grid` = :circular, it is the ring width for the final circular grid. In both cases of `reduce_grid`, the result will only be exact if `bin_size` divides `grid_size` exactly, otherwise `grid_size` will take priority and the final sizes will only approximate `bin_size`. For the cubic grids a default value of 200 pc is always used.
  - `plot_type::Symbol=:scatter`: If the plot will be a :scatter plot or a :heatmap. Heatmaps will not show legends, experimental measurements or several simulations at once.
  - `integrated::Bool=false`: If the integrated (one point per galaxy) or resolved (several point per galaxy) Kennicutt-Schmidt law will be plotted. `integrated` = true only works with `plot_type` = :scatter. The central value is the weighted median and the error bars are the median absolute deviations.
  - `sfr_density::Bool=true`: If the quantity for the y axis will be the SFR surface density or, if set to false, the stellar mass surface density.
  - `gas_weights::Union{Symbol,Nothing}=nothing`: If `plot_type` = :scatter, each point (a bin in the 2D grid) can be weighted by a gas quantity. If `integrated` = true, the median will be computed with these weights in mind. If `integrated` = false, each point will have a color given by the weight. The posible weights are:

      + `:gas_mass_density` -> Gas mass surface density of each bin. See the documentation for the function [`daDensity2DProjection`](@ref).
      + `:gas_sfr`          -> The total gas SFR of the column associated with each bin. See the documentation for the function [`daGasSFR2DProjection`](@ref).
      + `:gas_metallicity`  -> The total metallicity of the column associated with each bin. See the documentation for the function [`daMetallicity2DProjection`](@ref).
      + `:temperature`      -> The mean gas temperature of the column associated with each bin. See the documentation for the function [`daTemperature2DProjection`](@ref).
  - `measurements::Bool=true`: If the experimental measurements from Kennicutt (1998), Bigiel et al. (2008) or Bigiel et al. (2010) will be plotted alongside the simulation results.
  - `measurement_type::Union{String,Symbol}=:fits`: Type of measurement to plot, only valid if `measurement` = true. The option are:

      + `:fits`: Fits from Bigiel et al. (2008) and/or Kennicutt (1998) depending on the quantity in the x axis. The fits will be plotted as lines with uncertanty bands.
      + `"NGC XXX"`: Plot the resolved data of the given galaxy as a scatter plot. Uses the data from Sun et al. (2023). See the documentation of [`ppSun2023!`](@ref) for options.
      + `:main`: Plot the data of the main galaxy distribution in Sun et al. (2023), as a scatter plot.
      + `:all`: Plot the data of every galaxy in Sun et al. (2023), as a scatter plot.
  - `fit::Bool=false`: If the simulation data law will be fitted with a power law. The fit will be plotted as a line. This option is only valid if `integrated` = false and `plot_type` = :scatter, otherwise it will be ignored.
  - `x_range::Union{NTuple{2,<:Number},Nothing}=nothing`: x axis range for the heatmap grid. If set to `nothing`, the extrema of the x values will be used. Only relevant if `plot_type` = :heatmap.
  - `y_range::Union{NTuple{2,<:Number},Nothing}=nothing`: y axis range for the heatmap grid. If set to `nothing`, the extrema of the y values will be used. Only relevant if `plot_type` = :heatmap.
  - `n_bins::Int=100`: Number of bins per side of the heatmap grid. Only relevant if `plot_type` = :heatmap.
  - `colorbar::Bool=false`: If a colorbar will be added.
  - `output_file::String="./kennicutt_schmidt_law.png"`: Path to the output file.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Selects which cells/particles will be consider, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)

J. Sun et al. (2023). *Star Formation Laws and Efficiencies across 80 Nearby Galaxies*. The Astrophysical Journal Letters, **945(2)**, L19. [doi:10.3847/2041-8213/acbd9c](https://doi.org/10.3847/2041-8213/acbd9c)
"""
function kennicuttSchmidtLaw(
    simulation_paths::Vector{String},
    slice::IndexType;
    quantity::Symbol=:molecular_mass,
    gas_type::Symbol=:cells,
    reduce_grid::Symbol=:square,
    grid_size::Unitful.Length=BOX_L,
    bin_size::Unitful.Length=BIGIEL_PX_SIZE,
    plot_type::Symbol=:scatter,
    integrated::Bool=false,
    sfr_density::Bool=true,
    gas_weights::Union{Symbol,Nothing}=nothing,
    measurements::Bool=true,
    measurement_type::Union{String,Symbol}=:fits,
    fit::Bool=false,
    x_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    y_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    n_bins::Int=100,
    colorbar::Bool=false,
    output_file::String="./kennicutt_schmidt_law.png",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    # Compute the number of simulations
    ns = length(simulation_paths)

    ################################################################################################
    # Default values
    ################################################################################################

    # Default voxel side length
    voxel_size = 200.0u"pc"

    # Default units for the gas surface density
    Σg_m_unit = u"Msun"
    Σg_l_unit = u"kpc"

    # Default units for the stellar/sfr surface density
    Σs_m_unit = u"Msun"
    Σs_l_unit = u"kpc"
    Σs_t_unit = u"yr"

    ################################################################################################
    # Physical units
    ################################################################################################

    # Save the origial value of the global `PHYSICAL_UNITS`
    og_pu_value = PHYSICAL_UNITS

    if !og_pu_value && logging[]

        @warn("kennicuttSchmidtLaw: The global `PHYSICAL_UNITS` is set to false, \
        but Kennicutt-Schmidt law plots must be in physical units, so the global \
        setting will be ignored and default to true just for this function")

    end

    # Kennicutt-Schmidt law plots must be in physical units even for cosmological simulations
    global PHYSICAL_UNITS = true

    ################################################################################################
    # Check arguments
    ################################################################################################

    (
        quantity ∈ [:gas_mass, :molecular_mass, :br_molecular_mass, :neutral_mass] ||
        throw(ArgumentError("kennicuttSchmidtLaw: `quantity` can only be :gas_mass, \
        :molecular_mass, :br_molecular_mass or :neutral_mass, but I got :$(quantity)"))
    )

    (
        plot_type ∈ [:scatter, :heatmap] ||
        throw(ArgumentError("kennicuttSchmidtLaw: `plot_type` can only be :scatter or :heatmap, \
        but I got :$(plot_type)"))
    )

    (
        reduce_grid ∈ [:square, :circular] ||
        throw(ArgumentError("kennicuttSchmidtLaw: `reduce_grid` can only be :square or :circular, \
        but I got :$(reduce_grid)"))
    )

    if integrated

        if plot_type == :heatmap

            (
                !logging[] ||
                @warn("kennicuttSchmidtLaw: If `integrated` is set to true, `plot_type` = :heatmap \
                will be ignored and default to :scatter")
            )

            plot_type = :scatter

        end

        if fit

            (
                !logging[] ||
                @warn("kennicuttSchmidtLaw: `integrated` is set to :true, so \
                `fit` = true will be ignored and default to false")
            )

            fit = false

        end

    end

    if bin_size < voxel_size

        (
            !logging[] ||
            @warn("kennicuttSchmidtLaw: `reduce_grid` is set to :square and `bin_size` is set to \
            a value lower than $(voxel_size). This is not allowed. `bin_size` will be ignored and \
            default to $(voxel_size)")
        )

        bin_size = voxel_size

    end

    if bin_size > grid_size / 2.0

        (
            !logging[] ||
            @warn("kennicuttSchmidtLaw: `bin_size` is set to a value larger than \
            `grid_size` / 2 = $(grid_size / 2.0). This makes no sense. `bin_size` \
            will be ignored and default to $(BIGIEL_PX_SIZE)")
        )

        bin_size = BIGIEL_PX_SIZE

    end

    if !isnothing(sim_labels)

        # Compute the number of labels
        nl = length(sim_labels)

        (
            ns == nl ||
            throw(ArgumentError("kennicuttSchmidtLaw: `sim_labels` must have as many elements as \
            `simulation_paths`, but I got length(simulation_paths) = $(ns) \
            != length(sim_labels) = $(nl)"))
        )

    end

    if measurements

        if !sfr_density

            (
                !logging[] ||
                @warn("kennicuttSchmidtLaw: `sfr_density` is set to false, so \
                `measurements` = true will be ignored and default to false. The experimental \
                measurements are only for the SFR surface density")
            )

            measurements = false

        end

        if quantity == :gas_mass && logging[]

            @warn("kennicuttSchmidtLaw: The measurements (fits or otherwise) are only available \
            for molecular and neutral gas. For `quantity` = :gas_mass the neutral gas measurements \
            will be used, even though this is technically not correct")

        end

        if isa(measurement_type, String) && integrated && logging[]

            @warn("kennicuttSchmidtLaw: `integrated` is set to true but you have set \
            `measurement_type` to plot the resolved measurements of galaxy $(measurement_type). \
            Are you sure you want this?")

        end

        if measurement_type ∈ [:all, :main] && integrated && logging[]

            @warn("kennicuttSchmidtLaw: `integrated` is set to true but you have set \
            `measurement_type` to plot the resolved measurements of several galaxies in \
            Sun et al. (2023). Are you sure you want this?")

        end

        if measurement_type ∈ [:all, :main] && quantity ∈ [:gas_mass, :neutral_mass] && logging[]

            @warn("kennicuttSchmidtLaw: The measurements from Sun et al. (2023) are only for \
            molecular gas, but the selected quantity is :$(quantity). Are you sure you want this?")

        end

    end

    if plot_type == :heatmap

        if !isnothing(gas_weights)

            (
                !logging[] ||
                @warn("kennicuttSchmidtLaw: `plot_type` is set to :heatmap, so `gas_weights` = \
                :$(gas_weights) will be ignored and default to nothing")
            )

            gas_weights = nothing

        end

        if fit

            (
                !logging[] ||
                @warn("kennicuttSchmidtLaw: `plot_type` is set to :heatmap, so \
                `fit` = true will be ignored and default to false")
            )

            fit = false

        end

        if measurements

            (
                !logging[] ||
                @warn("kennicuttSchmidtLaw: `plot_type` is set to :heatmap, so \
                `measurements` = true will be ignored and default to false")
            )

            measurements = false

        end

        if ns > 1

            (
                !logging[] || @warn("kennicuttSchmidtLaw: `plot_type` is set to :heatmap, so only \
                one simulation at a time can be plotted, but I got length(simulation_paths) = \
                $(ns) > 1. `plot_type` = :heatmap will be ignored and default to :scatter")
            )

            plot_type = :scatter

        end

        if reduce_grid == :circular && logging[]

            @warn("kennicuttSchmidtLaw: `plot_type` is set to :heatmap and `reduce_grid` to \
            :circular. Are you sure you want this?")

        end

    end

    if colorbar && ((plot_type == :scatter && isnothing(gas_weights)) || integrated)

        (
            !logging[] ||
            @warn("kennicuttSchmidtLaw: `colorbar` is set to true, but there is no color range in \
            the plot (either `plot_type` = :scatter and `gas_weights` = nothing or `integrated` = \
            true). `colorbar` = true will be ignored and default to false")
        )

        colorbar = false

    end

    ################################################################################################
    # Compute grids
    ################################################################################################

    # Compute the number of bins for the high resolution grids
    hr_n_bins = round(Int, uconvert(Unitful.NoUnits, grid_size / voxel_size))

    if reduce_grid == :square

        # Compute the number of bins for the low resolution grids
        lr_n_bins = round(Int, uconvert(Unitful.NoUnits, grid_size / bin_size))

        # Compute the interger factor between the high resolution grids (`hr_n_bins`px)
        # and the low resolution grids (`lr_n_bins`px)
        reduce_factor = hr_n_bins ÷ lr_n_bins

        stellar_grid = CubicGrid(grid_size, reduce_factor * lr_n_bins)
        gas_grid     = CubicGrid(grid_size, reduce_factor * lr_n_bins)

    else

        stellar_grid = CubicGrid(grid_size, hr_n_bins)
        gas_grid     = CubicGrid(grid_size, hr_n_bins)

        # Compute the ring width for the circular grid
        reduce_factor = round(Int, uconvert(Unitful.NoUnits, (grid_size / 2.0) / bin_size))

    end

    ################################################################################################
    # Compute the density maps and save them as JLD2 files
    ################################################################################################

    # Set a folder for the JLD2 files
    temp_folder = joinpath(dirname(output_file), "_temp_jld2")

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        plotParams(:stellar_mass).request,
    )

    ##########################
    # Compute the stellar map
    ##########################

    plotSnapshot(
        simulation_paths,
        request,
        [heatmap!];
        output_path=temp_folder,
        base_filename=string(:stellar_mass),
        slice,
        filter_function,
        da_functions=[daDensity2DProjection],
        da_args=[(stellar_grid, :stellar_mass, :particles)],
        da_kwargs=[
            (;
                reduce_factor,
                reduce_grid,
                filter_function=dd->filterByStellarAge(dd),
                m_unit=Σs_m_unit,
                l_unit=Σs_l_unit,
            ),
        ],
        transform_box=true,
        translation,
        rotation,
        x_unit=u"kpc",
        y_unit=u"kpc",
        save_figures=false,
        backup_results=true,
    )

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        plotParams(quantity).request,
    )

    #############################
    # Compute the `quantity` map
    #############################

    plotSnapshot(
        simulation_paths,
        request,
        [heatmap!];
        output_path=temp_folder,
        base_filename=string(quantity),
        slice,
        filter_function,
        da_functions=[daDensity2DProjection],
        da_args=[(gas_grid, quantity, gas_type)],
        da_kwargs=[(; reduce_factor, reduce_grid, m_unit=Σg_m_unit, l_unit=Σg_l_unit,)],
        transform_box=true,
        translation,
        rotation,
        x_unit=u"kpc",
        y_unit=u"kpc",
        save_figures=false,
        backup_results=true,
    )

    ##########################
    # Compute the weights map
    ##########################

    if !isnothing(gas_weights)

        if gas_weights == :gas_mass_density

            da_function = daDensity2DProjection
            da_args     = [(gas_grid, :gas_mass, gas_type)]
            c_label     = L"\log_{10} \, \Sigma_\mathrm{gas} \,\, [\mathrm{M_\odot \, kpc^{-2}}]"

        elseif gas_weights == :gas_sfr

            da_function = daGasSFR2DProjection
            da_args     = [(gas_grid, gas_type)]
            c_label     = L"\log_{10} \, \mathrm{SFR_{gas} \,\, [M_\odot \, yr^{-1}]}"

        elseif gas_weights == :gas_metallicity

            da_function = daMetallicity2DProjection
            da_args     = [(gas_grid, :gas, gas_type)]
            c_label     = L"$\log_{10}$ %$(plotParams(:gas_metallicity).var_name)"

        elseif gas_weights == :temperature

            da_function = daTemperature2DProjection
            da_args     = [(gas_grid, gas_type)]
            c_label     = plotParams(:temperature).axis_label

        else

            throw(ArgumentError("kennicuttSchmidtLaw: `gas_weights` can only be \
            :gas_mass_density, :gas_sfr, :gas_metallicity or :temperature, but I got \
            :$(gas_weights)"))

        end

        filter_function, translation, rotation, request = selectFilter(
            filter_mode,
            mergeRequests(plotParams(gas_weights).request, plotParams(:gas_mass_density).request),
        )

        plotSnapshot(
            simulation_paths,
            request,
            [heatmap!];
            output_path=temp_folder,
            base_filename="gas_weights",
            slice,
            filter_function,
            da_functions=[da_function],
            da_args,
            da_kwargs=[(; reduce_factor, reduce_grid)],
            transform_box=true,
            translation,
            rotation,
            x_unit=u"kpc",
            y_unit=u"kpc",
            save_figures=false,
            backup_results=true,
        )

    end

    ################################################################################################
    # Set the axis labels
    ################################################################################################

    if quantity == :gas_mass

        x_label = getLabel(plotParams(:gas_area_density).var_name, 0, Σg_m_unit * Σg_l_unit^-2)

    elseif quantity == :molecular_mass

        x_label = getLabel(
            plotParams(:molecular_area_density).var_name, 0, Σg_m_unit * Σg_l_unit^-2,
        )

    elseif quantity == :br_molecular_mass

        x_label = getLabel(
            plotParams(:br_molecular_area_density).var_name, 0, Σg_m_unit * Σg_l_unit^-2,
        )

    elseif quantity == :neutral_mass

        x_label = getLabel(plotParams(:neutral_area_density).var_name, 0, Σg_m_unit * Σg_l_unit^-2)

    end

    if sfr_density

        y_label = getLabel(
            plotParams(:sfr_area_density).var_name, 0, Σs_m_unit * Σs_t_unit^-1 * Σs_l_unit^-2,
        )

    else

        y_label = getLabel(plotParams(:stellar_area_density).var_name, 0, Σs_m_unit * Σs_l_unit^-2)

    end

    ################################################################################################
    # Read and plot the data in the JLD2 files
    ################################################################################################

    if sfr_density
        # Factor to go from stellar surface density to SFR surface density
        # log10(Σsfr) = log10(Σ*) - log10Δt
        log10Δt = log10(ustrip(Σs_t_unit, AGE_RESOLUTION))
    end

    # Set the plot theme
    if integrated || reduce_grid == :circular
        markersize = 15
    else
        markersize = 6
    end

    current_theme = merge(
        theme,
        Theme(
            Legend=(
                nbanks=1,
                rowgap=-15,
                halign=:left,
                valign=:top,
                padding=(15, 0, 0, 0),
            ),
            Text=(fontsize=30,),
            Scatter=(; markersize),
        ),
        DEFAULT_THEME,
        theme_latexfonts(),
    )

    with_theme(current_theme) do

        f = Figure()

        ax = CairoMakie.Axis(
            f[1, 1];
            xlabel=L"$\log_{10}$ %$(x_label)",
            ylabel=L"$\log_{10}$ %$(y_label)",
        )

        colors = [:gray15, current_theme[:palette][:color][][2:ns]...]

        for (sim_idx, simulation) in pairs(simulation_paths)

            simulation_table = DataFrame(makeSimulationTable(simulation)[slice, :])
            sim_name         = "simulation_$(lpad(string(sim_idx), 3, "0"))"
            snapshot_numbers = simulation_table[!, :numbers]

            # Allocate memory for the heatmap. For heatmaps we need to accumulate the values
            # for every snapshot before plotting
            if plot_type == :heatmap
                x_heatmap = Float64[]
                y_heatmap = Float64[]
            end

            for snapshot_number in snapshot_numbers

                ############################################
                # Read the JLD2 files and sanitize the data
                ############################################

                ##############
                # Gas density
                ##############

                x_address = "$(quantity)_$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                x_file    = jldopen(joinpath(temp_folder, "$(string(quantity)).jld2"), "r")

                if reduce_grid == :square
                    x_data = vec(x_file[x_address][3])
                else
                    x_data = x_file[x_address][3]
                end

                x_idxs = map(x -> isnan(x) || iszero(x), x_data)

                ##############
                # SFR density
                ##############

                y_address = "stellar_mass_$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                y_file    = jldopen(joinpath(temp_folder, "stellar_mass.jld2"), "r")

                if reduce_grid == :square
                    y_data = vec(y_file[y_address][3])
                else
                    y_data = y_file[y_address][3]
                end

                y_idxs = map(x -> isnan(x) || iszero(x), y_data)

                delete_idxs = x_idxs ∪ y_idxs

                ##########
                # Weights
                ##########

                if !isnothing(gas_weights)

                    z_address = "gas_weights_$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                    z_file    = jldopen(joinpath(temp_folder, "gas_weights.jld2"), "r")

                    if reduce_grid == :square
                        z_data = vec(z_file[z_address][3])
                    else
                        z_data = z_file[z_address][3]
                    end

                    z_idxs = map(x -> isnan(x) || iszero(x), z_data)

                    delete_idxs = delete_idxs ∪ z_idxs

                    deleteat!(z_data, delete_idxs)

                end

                deleteat!(x_data, delete_idxs)
                deleteat!(y_data, delete_idxs)

                if sfr_density
                    y_data .-= log10Δt
                end

                # For the integrated Kennicutt-Schmidt law, compute the median and median absolute
                # deviation of the gas and stellar densities
                if integrated

                    lin_x = exp10.(x_data)
                    lin_y = exp10.(y_data)

                    if !isnothing(gas_weights)

                        w_z = weights(exp10.(z_data))

                        μx = median(lin_x, w_z)
                        μy = median(lin_y, w_z)

                    else

                        μx = median(lin_x)
                        μy = median(lin_y)

                    end

                    σx = mad(lin_x; center=μx, normalize=false)
                    σy = mad(lin_y; center=μy, normalize=false)

                    x_data = [log10(μx ± σx)]
                    y_data = [log10(μy ± σy)]

                end

                if plot_type == :heatmap

                    append!(x_heatmap, x_data)
                    append!(y_heatmap, y_data)

                end

                close(x_file)
                close(y_file)
                isnothing(gas_weights) || close(z_file)

                #################################
                # Plot the Kennicutt-Schmidt law
                #################################

                if  plot_type == :scatter

                    if integrated

                        scatter!(
                            ax,
                            Measurements.value.(x_data),
                            Measurements.value.(y_data);
                            color=colors[sim_idx],
                        )

                        errorbars!(
                            ax,
                            Measurements.value.(x_data),
                            Measurements.value.(y_data),
                            Measurements.uncertainty.(y_data);
                            color=colors[sim_idx],
                            direction=:y,
                        )

                        errorbars!(
                            ax,
                            Measurements.value.(x_data),
                            Measurements.value.(y_data),
                            Measurements.uncertainty.(x_data);
                            color=colors[sim_idx],
                            direction=:x,
                        )

                    else

                        if isnothing(gas_weights)

                            scatter!(
                                ax,
                                x_data,
                                y_data;
                                color=(colors[sim_idx], 0.8),
                            )

                            if fit
                                ppFitLine!(
                                    f;
                                    top_position=(0.7, 0.99),
                                    color=Makie.wong_colors()[1],
                                )
                            end

                        else

                            scatter!(
                                ax,
                                x_data,
                                y_data;
                                color=z_data,
                                colormap=:nipy_spectral,
                            )

                            if fit
                                ppFitLine!(
                                    f;
                                    top_position=(0.7, 0.99),
                                    wts=exp10.(z_data),
                                    color=Makie.wong_colors()[1],
                                )
                            end

                        end

                    end

                end

            end

            if plot_type == :heatmap

                # If there is no specified range, use the extrema of the x values
                if isnothing(x_range)
                    xrange = extrema(x_heatmap)
                else
                    xrange = x_range
                end

                # If there is no specified range, use the extrema of the y values
                if isnothing(y_range)
                    yrange = extrema(y_heatmap)
                else
                    yrange = y_range
                end

                # Compute the bin half width for each axis
                x_bin_h_width = 0.5 * (xrange[2] - xrange[1]) / n_bins
                y_bin_h_width = 0.5 * (yrange[2] - yrange[1]) / n_bins

                # Compute the center value of each bin for each axis
                x_axis = collect(
                    range(
                        xrange[1] + x_bin_h_width;
                        length=n_bins,
                        step=2 * x_bin_h_width,
                    ),
                )
                y_axis = collect(
                    range(
                        yrange[1] + y_bin_h_width;
                        length=n_bins,
                        step=2 * y_bin_h_width,
                    ),
                )

                # Compute the 2D histogram (number of pixels in each bin)
                values = histogram2D(
                    permutedims(hcat(x_heatmap, y_heatmap), (2, 1)),
                    collect(range(xrange[1], xrange[2]; length=n_bins + 1)),
                    collect(range(yrange[1], yrange[2]; length=n_bins + 1));
                )

                # The transpose and reverse operation are used to conform to the way heatmap!
                # expect the matrix to be structured
                z_axis = reverse!(transpose(values), dims=2)

                heatmap!(ax, x_axis, y_axis, z_axis; colormap=:nipy_spectral)

                if logging[]

                    clean_values = filter(!isnan, z_axis)

                    if isempty(clean_values)

                        min_max_v = (NaN, NaN)
                        mean_v    = NaN
                        meadian_v = NaN
                        mode_v    = NaN

                    else

                        min_max_v = extrema(clean_values)
                        mean_v    = mean(clean_values)
                        meadian_v = median(clean_values)
                        mode_v    = mode(clean_values)

                    end

                    # Print the count range
                    @info(
                        "\nCount range \
                        \n  Simulation: $(basename(simulation)) \
                        \n  Quantity:   $(quantity) \
                        \n  Min - Max:  $(min_max_v) \
                        \n  Mean:       $(mean_v) \
                        \n  Median:     $(meadian_v) \
                        \n  Mode:       $(mode_v)"
                    )

                end

            end

        end

        ############################################################################################
        # Plot the experimental fits and the legend
        ############################################################################################

        if !isnothing(sim_labels) && plot_type == :scatter

            if !isnothing(gas_weights) && !integrated

                markers = [
                    MarkerElement(; color=(colors[1], 0.8), marker=:circle, markersize=20) for
                    _ in eachindex(sim_labels)
                ]

            else

                markers = [
                    MarkerElement(; color, marker=:circle, markersize=20) for color in colors
                ]

            end

        end

        if measurements

            if measurement_type == :fits

                if quantity ∈ [:molecular_mass, :br_molecular_mass]

                    pp_legend = ppBigiel2008!(
                        f,
                        true;
                        x_unit=Σg_m_unit * Σg_l_unit^-2,
                        y_unit=Σs_m_unit * Σs_t_unit^-1 * Σs_l_unit^-2,
                        colors=[Makie.wong_colors()[1], Makie.wong_colors()[2]],
                    )

                else

                    legend_bigiel = ppBigiel2008!(
                        f,
                        false;
                        x_unit=Σg_m_unit * Σg_l_unit^-2,
                        y_unit=Σs_m_unit * Σs_t_unit^-1 * Σs_l_unit^-2,
                        colors=[Makie.wong_colors()[1], Makie.wong_colors()[2]],
                    )

                    legend_kennicut = ppKennicutt1998!(
                        f;
                        x_unit=Σg_m_unit * Σg_l_unit^-2,
                        y_unit=Σs_m_unit * Σs_t_unit^-1 * Σs_l_unit^-2,
                        colors=[Makie.wong_colors()[3], Makie.wong_colors()[4]],
                    )

                    pp_legend = (
                        vcat(legend_bigiel[1], legend_kennicut[1]),
                        vcat(legend_bigiel[2], legend_kennicut[2]),
                    )

                end

            else

                pp_legend = ppSun2023!(
                    f;
                    galaxy=measurement_type,
                    x_unit=Σg_m_unit * Σg_l_unit^-2,
                    y_unit=Σs_m_unit * Σs_t_unit^-1 * Σs_l_unit^-2,
                )

            end

            if !isnothing(sim_labels) && plot_type == :scatter

                Makie.Legend(f[1, 1], vcat(markers, pp_legend[1]), vcat(sim_labels, pp_legend[2]))

            end

        else

            if !isnothing(sim_labels) && plot_type == :scatter

                Makie.Legend(f[1, 1], markers, sim_labels)

            end

        end

        if !isnothing(sim_labels) && plot_type == :heatmap

            ppAnnotation!(f, sim_labels[1]; color=:white, fontsize=30)

        end

        ############################################################################################
        # Print the colorbar
        ############################################################################################

        if colorbar

            if plot_type == :heatmap
                Colorbar(f[1, 2], f.content[1].scene.plots[1]; label=L"\mathrm{Counts}")
            else
                Colorbar(f[1, 2], f.content[1].scene.plots[1]; label=c_label)
            end

            rowsize!(f.layout, 1, Makie.Fixed(pixelarea(ax.scene)[].widths[2]))

        end

        ############################################################################################
        # Save the plot
        ############################################################################################

        Makie.save(output_file, f)

    end

    rm(temp_folder; recursive=true)

    # Restore the original value of `PHYSICAL_UNITS`
    global PHYSICAL_UNITS = og_pu_value

    return nothing

end

"""
    fitVSFLaw(
        simulation_path::String,
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the resolved volumetric star formation (VSF) law with an optional linear fit.

!!! note

    Only stars younger than [`AGE_RESOLUTION`](@ref) are consider. The star formation surface density is just the stellar mass surface density divided by [`AGE_RESOLUTION`](@ref).

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: Quantity for the x axis. The options are:

      + `:gas_mass`          -> Gas density.
      + `:hydrogen_mass`     -> Hydrogen density.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) density.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) density.
      + `:ionized_mass`      -> Ionized hydrogen (``\\mathrm{HII}``) density.
      + `:neutral_mass`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) density.
  - `type::Symbol=:cells`: If the gas surface density will be calculated assuming the gas is in `:particles` or in Voronoi `:cells`.
  - `fit::Bool=true`: If a fit of the plotted values will be added on top of the scatter plot.
  - `box_size::Unitful.Length=BOX_L`: Physical side length for the grids
  - `x_range::NTuple{2,<:Real}=(-Inf, Inf)`: Only the data withing this range (for the x coordinates) will be fitted.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_label::Union{String,Nothing}=basename(simulation_path)`: Label for the scatter plot. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function fitVSFLaw(
    simulation_path::String,
    slice::IndexType;
    quantity::Symbol=:molecular_mass,
    type::Symbol=:cells,
    fit::Bool=true,
    box_size::Unitful.Length=BOX_L,
    x_range::NTuple{2,<:Real}=(-Inf, Inf),
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_label::Union{String,Nothing}=basename(simulation_path),
    theme::Attributes=Theme(),
)::Nothing

    grid = CubicGrid(box_size, 400)

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(plotParams(quantity).request, plotParams(:stellar_mass).request),
    )

    # Choose the correct x label
    if quantity == :gas_mass

        x_label = getLabel(
            plotParams(:gas_mass_density).var_name,
            0,
            u"Msun * pc^-3";
            latex=true,
        )

    elseif quantity == :molecular_mass

        x_label = getLabel(
            L"\rho_\mathrm{H_2}",
            0,
            u"Msun * pc^-3";
            latex=true,
        )

    elseif quantity == :br_molecular_mass

        x_label = getLabel(
            L"\rho_\mathrm{H_2^{BR}}",
            0,
            u"Msun * pc^-3";
            latex=true,
        )

    elseif quantity == :neutral_mass

        x_label = getLabel(
            L"\rho_\mathrm{H_I + H_2}",
            0,
            u"Msun * pc^-3";
            latex=true,
        )

    end

    # Set the y label
    y_label = getLabel(
        L"\rho_\mathrm{SFR}",
        0,
        u"Msun * yr^-1 * kpc^-3";
        latex=true,
    )

    plotSnapshot(
        [simulation_path],
        request,
        [scatter!];
        pf_kwargs=[(; color=Makie.wong_colors()[1], markersize=6, marker=:circle)],
        output_path,
        base_filename="vsf_law",
        output_format=".png",
        slice,
        filter_function,
        da_functions=[daVSFLaw],
        da_args=[(grid, quantity)],
        da_kwargs=[(; type, stellar_ff=dd->filterByStellarAge(dd))],
        post_processing=fit ? ppFitLine! : getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        x_trim=x_range,
        xaxis_label=L"$\log_{10}$ %$(x_label)",
        yaxis_label=L"$\log_{10}$ %$(y_label)",
        xaxis_var_name="",
        yaxis_var_name="",
        theme,
        sim_labels=[sim_label],
    )

    return nothing

end

"""
    massMetallicityRelation(
        simulation_paths::Vector{String},
        slice::IndexType;
        <keyword arguments>
    )::Nothing

Plot the resolved mass-metallicity relation. This method plots the M-Z relation at a fix moment in time.

!!! note

    Only stars younger than [`AGE_RESOLUTION`](@ref) and gas cells/particles within a sphere of radius `DISK_R` are consider.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `element::Symbol=:all`: Which metallicity to use. The options are:

      + `:all` -> Metallicity considering all elements, as ``Z / Z_\\odot``.
      + `:X`   -> Element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
  - `mass::Bool=true`: If the x axis will be the stellar mass density or the SFR density.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection, averaging the value of neighboring pixels. It has to divide the size of `grid` exactly.
  - `output_path::String="./resolvedKSLawZScatter"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).
"""
function massMetallicityRelation(
    simulation_paths::Vector{String},
    slice::IndexType;
    element::Symbol=:all,
    mass::Bool=true,
    reduce_factor::Int=1,
    output_path::String="./massMetallicityRelation",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    grid = CubicGrid(BOX_L, 400)

    (
        element ∈ [:all, keys(ELEMENT_INDEX)...] ||
        throw(ArgumentError("massMetallicityRelation: `quantity` can only be :all or any \
        of the keys of `ELEMENT_INDEX` (see `./src/constants/globals.jl`), but I got :$(quantity)"))
    )

    # Set a temporal folder for the JLD2 files
    temp_folder = joinpath(output_path, "_temp_jld2")

    # Write the JLD2 files with the density maps
    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        plotParams(:stellar_mass).request,
    )

    plotSnapshot(
        simulation_paths,
        request,
        [heatmap!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path=temp_folder,
        base_filename=string(:stellar_mass),
        output_format=".png",
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daDensity2DProjection],
        da_args=[(grid, :stellar_mass, :particles)],
        da_kwargs=[(; reduce_factor, filter_function=dd->filterByStellarAge(dd))],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=u"kpc",
        y_unit=u"kpc",
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label="auto_label",
        yaxis_label="auto_label",
        xaxis_var_name="x",
        yaxis_var_name="y",
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        # Plotting options
        save_figures=false,
        backup_results=true,
        theme=Theme(),
        sim_labels=nothing,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=15,
    )

    if element == :all
        metal_request = plotParams(:gas_metallicity).request
    else
        metal_request = plotParams(Symbol(element, "_gas_abundance")).request
    end

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(plotParams(:gas_mass_density).request, metal_request),
    )

    # Write the JLD2 file with the metallicity density
    plotSnapshot(
        simulation_paths,
        request,
        [heatmap!];
        pf_kwargs=[(;)],
        # `plotSnapshot` configuration
        output_path=temp_folder,
        base_filename="gas_metallicity",
        output_format=".png",
        show_progress=true,
        # Data manipulation options
        slice,
        filter_function,
        da_functions=[daMetallicity2DProjection],
        da_args=[(grid, :gas, :cells)],
        da_kwargs=[(; element)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth=0,
        x_unit=u"kpc",
        y_unit=u"kpc",
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim=(-Inf, Inf),
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label="auto_label",
        yaxis_label="auto_label",
        xaxis_var_name="x",
        yaxis_var_name="y",
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        # Plotting options
        save_figures=false,
        backup_results=true,
        theme=Theme(),
        sim_labels=nothing,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=15,
    )

    # Set the x label
    if mass
        x_label = getLabel(
            plotParams(:stellar_area_density).var_name,
            0,
            u"Msun * kpc^-2";
            latex=true,
        )
    else
        x_label = getLabel(
            plotParams(:sfr_area_density).var_name,
            0,
            u"Msun * yr^-1 * kpc^-2";
            latex=true,
        )
    end

    # Set the y label
    if element == :all
        ylabel = L"$\log_{10}$ %$(plotParams(:gas_metallicity).var_name)"
    else
        ylabel = plotParams(Symbol(element, "_gas_abundance")).axis_label
    end

    with_theme(theme, merge(DEFAULT_THEME, theme_latexfonts())) do

        for (sim_idx, simulation) in pairs(simulation_paths)

            simulation_table = DataFrame(makeSimulationTable(simulation)[slice, :])
            sim_name         = "simulation_$(lpad(string(sim_idx), 3, "0"))"
            times            = ustrip.(u"Gyr", simulation_table[!, :physical_times])
            snapshot_numbers = simulation_table[!, :numbers]

            for (time, snapshot_number) in zip(times, snapshot_numbers)

                f = Figure()

                ax = CairoMakie.Axis(
                    f[1, 1];
                    xlabel=L"$\log_{10}$ %$(x_label)",
                    ylabel,
                    aspect=AxisAspect(1),
                )

                x_address = "stellar_mass_$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"
                y_address = "gas_metallicity_$(SNAP_BASENAME)_$(snapshot_number)/$(sim_name)"

                jldopen(joinpath(temp_folder, "stellar_mass.jld2"), "r") do x_file

                    jldopen(joinpath(temp_folder, "gas_metallicity.jld2"), "r") do y_file

                        # Read the JLD2 files
                        x_data = vec(x_file[x_address][3])
                        y_data = vec(y_file[y_address][3])

                        # Delete 0s and NaNs in the data vectors
                        x_idxs = map(x -> isnan(x) || iszero(x), x_data)
                        y_idxs = map(x -> isnan(x) || iszero(x), y_data)

                        deleteat!(x_data, x_idxs ∪ y_idxs)
                        deleteat!(y_data, x_idxs ∪ y_idxs)

                        if mass
                            x_axis = x_data
                        else
                            x_axis = x_data .- log10(ustrip(u"yr", AGE_RESOLUTION))
                        end

                        scatter!(
                            ax,
                            x_axis,
                            y_data;
                            markersize=6,
                            color=Makie.wong_colors()[6],
                        )

                    end

                end

                ppAnnotation!(
                    f,
                    L"t = %$(rpad(round(time, sigdigits=3), 4, '0')) \, \mathrm{Gyr}",
                    position=(0.04, 0.98),
                    fontsize=25,
                )

                if !isnothing(sim_labels)
                    Makie.Legend(
                        f[1, 1],
                        [
                            MarkerElement(;
                                color=Makie.wong_colors()[6],
                                marker=:circle,
                                markersize=20,
                            ),
                        ],
                        sim_labels[sim_idx:sim_idx],
                        nbanks=1,
                        labelsize=22,
                        rowgap=-20,
                        halign=:left,
                        valign=:top,
                        padding=(13, 0, 0, 45),
                        patchlabelgap=2,
                    )
                end

                path = mkpath(joinpath(output_path, basename(simulation)))

                Makie.save(joinpath(path, "$(snapshot_number).png"), f)

            end

        end

    end

    rm(temp_folder; recursive=true)

    return nothing

end

"""
    gasVelocityCubes(
        simulation_paths::Vector{String},
        slice::ReducedIndexType;
        <keyword arguments>
    )::Nothing

Create a HDF5 file with the position, gas mass, velocity, and velocity dispersion of each voxel in a rectangular 3D grid.

The metadata for each snapshot in the HDF5 file includes the physical time in Gyr, the scale factor, and the redshift of that snapshot.

By default, the grid is centered at coordinates (0, 0, 0), has 300x300x300 voxels, and has a side length of [`BOX_L`](@ref). There are as many rows as there are voxels (27000000 by default).

The quantities in the HDF5 file for each voxel are:

Column 01: x coordinate [`l_unit`]
Column 02: y coordinate [`l_unit`]
Column 03: z coordinate [`l_unit`]
Column 04: Molecular mass [`m_unit`]
Column 05: Atomic mass [`m_unit`]
Column 06: Ionized mass [`m_unit`]
Column 07: Velocity in the x direction [`v_unit`]
Column 08: Velocity in the y direction [`v_unit`]
Column 09: Velocity in the z direction [`v_unit`]
Column 10: Velocity dispersion in the x direction [`v_unit`]
Column 11: Velocity dispersion in the y direction [`v_unit`]
Column 12: Velocity dispersion in the z direction [`v_unit`]

For gas represented by Voronoi cells (e.g. Arepo):

The mass is the mass of molecular, atomic or ionized gas intersecting the voxel, so it only considers the cell that is closest to the voxel. The velocity is given by the weighted mean of the velocities of the `n_neighbors` nearest cells. And the velocity dispersion, by the weighted standard deviation.

Notice that for Voronoi cells, the mass will be sample at a sub-cell resolution (as long as voxel size < cell size), while the velocities are sample at a locally lower-than-cell resolution (as long as `n_neighbors` > 1). The weights are given by the distance (in kpc) to each neighbor, using a Gaussian kernel.

For gas represented by particles (e.g. SPH codes):

The mass is the accumulated mass of the particles within each voxel. The velocity is the mean of the velocities of those particles, and the velocity dispersion is the standard deviation.

If there are no particles, the mass is 0, and the velocity and velocity dispersion are NaN. If there is only one particle, the mass and velocity are the ones from that particle, and the velocity dispersion is NaN.

By default (`filter_mode` = :subhalo) we use the following reference system:

  - The origin is in the position of the particle/cell at the potencial minimum of the main subhalo.
  - The x, y, and z axis form a right-handed cartesian reference system (x × y = z), where the z axis has the orientation of the stellar angular momentum, and the x and y axis are roughly in the direction of the corresponding principal axis.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::ReducedIndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13) or an `StepRange` (e.g. 5:2:13). Starts at 1.
  - `type::Symbol=:cells`: If the gas density will be calculated assuming the gas is in `:particles` or in Voronoi `:cells`.
  - `n_neighbors::Int=32`: Number of neighbors for the mean and standard deviation of the velocity. Setting this value to 1 maximizes the resolution for the velocity, and sets the standard deviation (columns 8, 9, and 10) to NaN. This is only relevant for simulations where gas is represented by Voronoi cells (`type` = :cells).
  - `grid::CubicGrid=CubicGrid(BOX_L, 300)`: Cubic grid.
  - `row_major_order::Bool=true`: Store the results in row-major order (C and Python) instead of column-major order (Julia, Fortran, and MATLAB). See [Row- and column-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order).
  - `m_unit::Unitful.Units=u"Msun"`: Mass unit
  - `l_unit::Unitful.Units=u"kpc"`: Length unit.
  - `v_unit::Unitful.Units=u"km * s^-1"`: Velocity unit.
  - `output_file::String="./gas_velocity_cube.hdf5"`: Path to the output HDF5 file. This file will be created, and the full path to it too, if it doesn't exist.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:subhalo`: Which cells/particles will be consider, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `show_progress::Bool=true`: If a progress bar will be shown.
"""
function gasVelocityCubes(
    simulation_paths::Vector{String},
    slice::ReducedIndexType;
    type::Symbol=:cells,
    n_neighbors::Int=32,
    grid::CubicGrid=CubicGrid(BOX_L, 300),
    row_major_order::Bool=true,
    m_unit::Unitful.Units=u"Msun",
    l_unit::Unitful.Units=u"kpc",
    v_unit::Unitful.Units=u"km * s^-1",
    output_file::String="./gas_velocity_cube.hdf5",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:subhalo,
    show_progress::Bool=true,
)::Nothing

    # Set the number of columns and rows
    n_rows = grid.n_bins^3
    n_cols = 12

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(
            plotParams(:molecular_mass).request,
            plotParams(:atomic_mass).request,
            plotParams(:ionized_mass).request,
            Dict(:gas => ["POS ", "VEL ", "RHO ", "MASS"]),
        ),
    )

    # For gas cells, reshape the grid to conform to the way `knn` expect the matrix to be structured
    if type == :cells

        physical_grid = Matrix{Float64}(undef, 3, n_rows)

        for i in eachindex(grid.grid)
            physical_grid[1, i] = ustrip(l_unit, grid.grid[i][1])
            physical_grid[2, i] = ustrip(l_unit, grid.grid[i][2])
            physical_grid[3, i] = ustrip(l_unit, grid.grid[i][3])
        end

    end

    # Create the output folder
    mkpath(dirname(output_file))

    # Create the output HDF5 file
    hdf5_file = h5open(output_file, "w")

    for simulation_path in simulation_paths

        simulation_name = basename(simulation_path)

        prog_bar = Progress(
            length(slice),
            dt=0.5,
            desc="Writing the velocity cube for simulation $(simulation_name)... ",
            color=:blue,
            barglyphs=BarGlyphs("|#  |"),
            enabled=show_progress,
        )

        # Create an HDF5 group for each simulation
        hdf5_group = create_group(hdf5_file, simulation_name)

        for snap_n in slice

            data_dict = makeDataDict(simulation_path, snap_n, request)

            snapshot_number = lpad(string(data_dict[:snap_data].global_index), 3, "0")

            # Filter the data
            filterData!(data_dict; filter_function)

            # Translate the data
            translateData!(data_dict, translation)

            # Rotate the data
            rotateData!(data_dict, rotation)

            # Load the gas quantities
            gd = data_dict[:gas]

            # Load the cell/particle positions
            positions = gd["POS "]

            # Load the cell/particle velocities
            velocities = ustrip.(v_unit, gd["VEL "])

            # Compute the mass of molecular, atomic, and ionized gas in each cell
            mol_masses = scatterQty(data_dict, :molecular_mass)
            ato_masses = scatterQty(data_dict, :atomic_mass)
            ion_masses = scatterQty(data_dict, :ionized_mass)

            if any(isempty, [mol_masses, ato_masses, ion_masses, velocities, positions])
                throw(ArgumentError("gasVelocityCubes: Some data is missing (there appears to be \
                no gas in the snapshot), so I cannot construct the velocity cube"))
            end

            # Alocate memory for:
            # Column 01: x coordinate [l_unit]
            # Column 02: y coordinate [l_unit]
            # Column 03: z coordinate [l_unit]
            # Column 04: Molecular mass [m_unit]
            # Column 05: Atomic mass [m_unit]
            # Column 06: Ionized mass [m_unit]
            # Column 07: Velocity in the x direction [v_unit]
            # Column 08: Velocity in the y direction [v_unit]
            # Column 09: Velocity in the z direction [v_unit]
            # Column 10: Velocity dispersion in the x direction [v_unit]
            # Column 11: Velocity dispersion in the y direction [v_unit]
            # Column 12: Velocity dispersion in the z direction [v_unit]
            data_matrix = Matrix{Float64}(undef, n_rows, n_cols)

            if type == :cells

                # Compute the volume of each cell
                cell_volumes = gd["MASS"] ./ gd["RHO "]

                # Compute the gas densities
                mol_densities = ustrip.(m_unit * l_unit^-3, mol_masses ./ cell_volumes)
                ato_densities = ustrip.(m_unit * l_unit^-3, ato_masses ./ cell_volumes)
                ion_densities = ustrip.(m_unit * l_unit^-3, ion_masses ./ cell_volumes)

                # Load the volume of the voxels
                voxel_volume = ustrip(l_unit^3, grid.bin_volume)

                # Compute the tree for a nearest neighbor search
                kdtree = KDTree(ustrip.(l_unit, positions))

                # Find the `n_neighbors` nearest cells to each voxel
                idxs, dists = knn(kdtree, physical_grid, n_neighbors, true)

                for i in eachindex(grid.grid)

                    # Physical coordinates of the voxel [l_unit]
                    data_matrix[i, 1:3] .= ustrip.(l_unit, grid.grid[i])

                    # Molecular gas mass [m_unit]
                    data_matrix[i, 4] = mol_densities[idxs[i][1]] * voxel_volume
                    # Atomic gas mass [m_unit]
                    data_matrix[i, 5] = ato_densities[idxs[i][1]] * voxel_volume
                    # Ionized gas mass [m_unit]
                    data_matrix[i, 6] = ion_densities[idxs[i][1]] * voxel_volume

                    if isone(n_neighbors)

                        # Neighbor velocity in the x direction [v_unit]
                        data_matrix[i, 7] = velocities[1, idxs[i]]
                        # Neighbor velocity in the y direction [v_unit]
                        data_matrix[i, 8] = velocities[2, idxs[i]]
                        # Neighbor velocity in the z direction [v_unit]
                        data_matrix[i, 9] = velocities[3, idxs[i]]

                        # For the case of only one neighbor, set the standard deviations to NaN
                        data_matrix[i, 10]  = NaN
                        data_matrix[i, 11]  = NaN
                        data_matrix[i, 12] = NaN

                    else

                        # Compute the analytic weights using a Gaussian kernel
                        neighbor_weights = aweights(evaluateNormal(dists[i]))

                        # Neighbor velocities in the x direction [v_unit]
                        vxs = velocities[1, idxs[i]]
                        # Neighbor velocities in the y direction [v_unit]
                        vys = velocities[2, idxs[i]]
                        # Neighbor velocities in the z direction [v_unit]
                        vzs = velocities[3, idxs[i]]

                        # Mean and standard deviation of the neighbor velocities in the x direction [v_unit]
                        data_matrix[i, 7], data_matrix[i, 10] = mean_and_std(vxs, neighbor_weights)
                        # Mean and standard deviation of the neighbor velocities in the y direction [v_unit]
                        data_matrix[i, 8], data_matrix[i, 11] = mean_and_std(vys, neighbor_weights)
                        # Mean and standard deviation of the neighbor velocities in the z direction [v_unit]
                        data_matrix[i, 9], data_matrix[i, 12] = mean_and_std(vzs, neighbor_weights)

                    end

                end

            elseif type == :particles

                # Find which particles are within each voxel
                idxs = listHistogram3D(positions, grid)

                for i in eachindex(grid.grid)

                    # Physical coordinates of the voxel [l_unit]
                    data_matrix[i, 1:3] .= ustrip.(l_unit, grid.grid[i])

                    # Molecular gas mass [m_unit]
                    data_matrix[i, 4] = ustrip(m_unit, sum(mol_masses[idxs[i]]; init=0.0*m_unit))
                    # Atomic gas mass [m_unit]
                    data_matrix[i, 5] = ustrip(m_unit, sum(ato_masses[idxs[i]]; init=0.0*m_unit))
                    # Ionized gas mass [m_unit]
                    data_matrix[i, 6] = ustrip(m_unit, sum(ion_masses[idxs[i]]; init=0.0*m_unit))

                    if isempty(idxs[i])

                        # If the voxel has no particles set the velocity to NaN
                        data_matrix[i, 7] = NaN
                        data_matrix[i, 8] = NaN
                        data_matrix[i, 9] = NaN

                        # If the voxel has no particles set the velocity dispersion to NaN
                        data_matrix[i, 10] = NaN
                        data_matrix[i, 11] = NaN
                        data_matrix[i, 12] = NaN

                    elseif isone(length(idxs[i]))

                        # Velocity in the x direction [v_unit]
                        data_matrix[i, 7] = velocities[1, idxs[i][1]]
                        # Velocity in the y direction [v_unit]
                        data_matrix[i, 8] = velocities[2, idxs[i][1]]
                        # Velocity in the z direction [v_unit]
                        data_matrix[i, 9] = velocities[3, idxs[i][1]]

                        # If the voxel has a single particle set the velocity dispersion to NaN
                        data_matrix[i, 10] = NaN
                        data_matrix[i, 11] = NaN
                        data_matrix[i, 12] = NaN

                    else

                        # Velocities in the x direction of the particles within the voxel [v_unit]
                        vxs = velocities[1, idxs[i]]
                        # Velocities in the y direction of the particles within the voxel [v_unit]
                        vys = velocities[2, idxs[i]]
                        # Velocities in the z direction of the particles within the voxel [v_unit]
                        vzs = velocities[3, idxs[i]]

                        # Mean and standard deviation of the velocities in the x direction [v_unit]
                        data_matrix[i, 7], data_matrix[i, 10] = mean_and_std(vxs)
                        # Mean and standard deviation of the velocities in the y direction [v_unit]
                        data_matrix[i, 8], data_matrix[i, 11] = mean_and_std(vys)
                        # Mean and standard deviation of the velocities in the z direction [v_unit]
                        data_matrix[i, 9], data_matrix[i, 12] = mean_and_std(vzs)

                    end

                end

            else

                throw(ArgumentError("gasVelocityCubes: The argument `type` must be :cells or \
                :particles, but I got :$(type)"))

            end

            if row_major_order
                # Go from column-major order (Julia, MATLAB, and Fortran) to
                # row-major order (Python and C), for interoperability
                hdf5_group["snap_$(snapshot_number)", shuffle=(), deflate=5] = permutedims(
                    data_matrix,
                    reverse(1:ndims(data_matrix)),
                )
            else
                # Stay in column-major order (Julia, MATLAB, and Fortran)
                hdf5_group["snap_$(snapshot_number)", shuffle=(), deflate=5] = data_matrix
            end

            # Read the time, scale factor, and redshift
            pt = ustrip.(u"Gyr", data_dict[:snap_data].physical_time)
            sf = data_dict[:snap_data].scale_factor
            rs = data_dict[:snap_data].redshift

            # Write the time metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Time [Gyr]"]   = pt
            attrs(hdf5_group["snap_$(snapshot_number)"])["Scale factor"] = sf
            attrs(hdf5_group["snap_$(snapshot_number)"])["Redshift"]     = rs

            # Write the unit metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Mass unit"]     = string(m_unit)
            attrs(hdf5_group["snap_$(snapshot_number)"])["Length unit"]   = string(l_unit)
            attrs(hdf5_group["snap_$(snapshot_number)"])["Velocity unit"] = string(v_unit)

            # Write the grid metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Grid size [length unit]"] = ustrip.(
                l_unit,
                grid.physical_size,
            )
            attrs(hdf5_group["snap_$(snapshot_number)"])["Grid size [# voxels]"] = grid.n_bins

            # Write the column metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Columns"] = [
                "x",    # Column 01: x coordinate [l_unit]
                "y",    # Column 02: y coordinate [l_unit]
                "z",    # Column 03: z coordinate [l_unit]
                "MH2",  # Column 04: Molecular mass [m_unit]
                "MHI",  # Column 05: Atomic mass [m_unit]
                "MHII", # Column 06: Ionized mass [m_unit]
                "Vx",   # Column 07: Velocity in the x direction [v_unit]
                "Vy",   # Column 08: Velocity in the y direction [v_unit]
                "Vz",   # Column 09: Velocity in the z direction [v_unit]
                "Sx",   # Column 10: Velocity dispersion in the x direction [v_unit]
                "Sy",   # Column 11: Velocity dispersion in the y direction [v_unit]
                "Sz",   # Column 12: Velocity dispersion in the z direction [v_unit]
            ]

            next!(prog_bar)

        end

    end

    close(hdf5_file)

    return nothing

end

"""
    stellarVelocityCubes(
        simulation_paths::Vector{String},
        slice::ReducedIndexType;
        <keyword arguments>
    )::Nothing

Create a HDF5 file with the position, stellar mass, velocity, and velocity dispersion of each voxel in a rectangular 3D grid.

The metadata for each snapshot in the HDF5 file includes the physical time in Gyr, the scale factor, and the redshift of that snapshot.

By default, the grid is centered at coordinates (0, 0, 0), has 100x100x100 voxels, and has a side length of [`BOX_L`](@ref). There are as many rows as there are voxels (1000000 by default).

The quantities in the HDF5 file for each voxel are:

Column 01: x coordinate [`l_unit`]
Column 02: y coordinate [`l_unit`]
Column 03: z coordinate [`l_unit`]
Column 04: Stellar mass [`m_unit`]
Column 05: Velocity in the x direction [`v_unit`]
Column 06: Velocity in the y direction [`v_unit`]
Column 07: Velocity in the z direction [`v_unit`]
Column 08: Velocity dispersion in the x direction [`v_unit`]
Column 09: Velocity dispersion in the y direction [`v_unit`]
Column 10: Velocity dispersion in the z direction [`v_unit`]

The mass is the accumulated mass of the particles within each voxel. The velocity is the mean of the velocities of those particles, and the velocity dispersion is the standard deviation.

If there are no particles, the mass is 0, and the velocity and velocity dispersion are NaN. If there is only one particle, the mass and velocity are the ones from that particle, and the velocity dispersion is NaN.

By default (`filter_mode` = :subhalo) we use the following reference system:

  - The origin is in the position of the particle/cell at the potencial minimum of the main subhalo.
  - The x, y, and z axis form a right-handed cartesian reference system (x × y = z), where the z axis has the orientation of the stellar angular momentum, and the x and y axis are roughly in the direction of the corresponding principal axis.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::ReducedIndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13) or an `StepRange` (e.g. 5:2:13). Starts at 1.
  - `grid::CubicGrid=CubicGrid(BOX_L, 100)`: Cubic grid.
  - `row_major_order::Bool=true`: Store the results in row-major order (C and Python) instead of column-major order (Julia, Fortran, and MATLAB). See [Row- and column-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order).
  - `m_unit::Unitful.Units=u"Msun"`: Mass unit
  - `l_unit::Unitful.Units=u"kpc"`: Length unit.
  - `v_unit::Unitful.Units=u"km * s^-1"`: Velocity unit.
  - `output_file::String="./stellar_velocity_cube.hdf5"`: Path to the output HDF5 file. This file will be created, and the full path to it too, if it doesn't exist.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:subhalo`: Which cells/particles will be consider, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `show_progress::Bool=true`: If a progress bar will be shown.
"""
function stellarVelocityCubes(
    simulation_paths::Vector{String},
    slice::ReducedIndexType;
    grid::CubicGrid=CubicGrid(BOX_L, 100),
    row_major_order::Bool=true,
    m_unit::Unitful.Units=u"Msun",
    l_unit::Unitful.Units=u"kpc",
    v_unit::Unitful.Units=u"km * s^-1",
    output_file::String="./stellar_velocity_cube.hdf5",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:subhalo,
    show_progress::Bool=true,
)::Nothing

    # Set the number of columns and rows
    n_rows = grid.n_bins^3
    n_cols = 10

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(plotParams(:stellar_mass).request, Dict(:stars => ["POS ", "VEL "])),
    )

    # Create the output folder
    mkpath(dirname(output_file))

    # Create the output HDF5 file
    hdf5_file = h5open(output_file, "w")

    for simulation_path in simulation_paths

        simulation_name = basename(simulation_path)

        prog_bar = Progress(
            length(slice),
            dt=0.5,
            desc="Writing the velocity cube for simulation $(simulation_name)... ",
            color=:blue,
            barglyphs=BarGlyphs("|#  |"),
            enabled=show_progress,
        )

        # Create an HDF5 group for each simulation
        hdf5_group = create_group(hdf5_file, simulation_name)

        for snap_n in slice

            data_dict = makeDataDict(simulation_path, snap_n, request)

            snapshot_number = lpad(string(data_dict[:snap_data].global_index), 3, "0")

            # Filter the data
            filterData!(data_dict; filter_function)

            # Translate the data
            translateData!(data_dict, translation)

            # Rotate the data
            rotateData!(data_dict, rotation)

            # Load the particle positions
            positions = data_dict[:stars]["POS "]

            # Load the stellar velocities
            velocities = ustrip.(v_unit, data_dict[:stars]["VEL "])

            # Compute the stellar mass in each cell
            masses = scatterQty(data_dict, :stellar_mass)

            if any(isempty, [masses, velocities, positions])
                throw(ArgumentError("stellarVelocityCubes: Some data is missing (there appears to \
                be no stars in the snapshot), so I cannot construct the velocity cube"))
            end

            # Alocate memory for:
            # Column 01: x coordinate [l_unit]
            # Column 02: y coordinate [l_unit]
            # Column 03: z coordinate [l_unit]
            # Column 04: Stellar mass [m_unit]
            # Column 05: Velocity in the x direction [v_unit]
            # Column 06: Velocity in the y direction [v_unit]
            # Column 07: Velocity in the z direction [v_unit]
            # Column 08: Velocity dispersion in the x direction [v_unit]
            # Column 09: Velocity dispersion in the y direction [v_unit]
            # Column 10: Velocity dispersion in the z direction [v_unit]
            data_matrix = Matrix{Float64}(undef, n_rows, n_cols)

            # Find which particles are within each voxel
            idxs = listHistogram3D(positions, grid)

            for i in eachindex(grid.grid)

                # Physical coordinates of the voxel [l_unit]
                data_matrix[i, 1:3] .= ustrip.(l_unit, grid.grid[i])

                # Stellar mass [m_unit]
                data_matrix[i, 4] = ustrip(m_unit, sum(masses[idxs[i]]; init=0.0*m_unit))

                if isempty(idxs[i])

                    # If the voxel has no particles set the velocity to NaN
                    data_matrix[i, 5] = NaN
                    data_matrix[i, 6] = NaN
                    data_matrix[i, 7] = NaN

                    # If the voxel has no particles set the velocity dispersion to NaN
                    data_matrix[i, 8]  = NaN
                    data_matrix[i, 9]  = NaN
                    data_matrix[i, 10] = NaN

                elseif isone(length(idxs[i]))

                    # Velocity in the x direction [v_unit]
                    data_matrix[i, 5] = velocities[1, idxs[i][1]]
                    # Velocity in the y direction [v_unit]
                    data_matrix[i, 6] = velocities[2, idxs[i][1]]
                    # Velocity in the z direction [v_unit]
                    data_matrix[i, 7] = velocities[3, idxs[i][1]]

                    # If the voxel has a single particle set the velocity dispersion to NaN
                    data_matrix[i, 8]  = NaN
                    data_matrix[i, 9]  = NaN
                    data_matrix[i, 10] = NaN

                else

                    # Velocities in the x direction of the particles within the voxel [v_unit]
                    vxs = velocities[1, idxs[i]]
                    # Velocities in the y direction of the particles within the voxel [v_unit]
                    vys = velocities[2, idxs[i]]
                    # Velocities in the z direction of the particles within the voxel [v_unit]
                    vzs = velocities[3, idxs[i]]

                    # Mean and standard deviation of the velocities in the x direction [v_unit]
                    data_matrix[i, 5], data_matrix[i, 8] = mean_and_std(vxs)
                    # Mean and standard deviation of the velocities in the y direction [v_unit]
                    data_matrix[i, 6], data_matrix[i, 9] = mean_and_std(vys)
                    # Mean and standard deviation of the velocities in the z direction [v_unit]
                    data_matrix[i, 7], data_matrix[i, 10] = mean_and_std(vzs)

                end

            end

            if row_major_order
                # Go from column-major order (Julia, MATLAB, and Fortran) to
                # row-major order (Python and C), for interoperability
                hdf5_group["snap_$(snapshot_number)", shuffle=(), deflate=5] = permutedims(
                    data_matrix,
                    reverse(1:ndims(data_matrix)),
                )
            else
                # Stay in column-major order (Julia, MATLAB, and Fortran)
                hdf5_group["snap_$(snapshot_number)", shuffle=(), deflate=5] = data_matrix
            end

            # Read the time, scale factor, and redshift
            pt = ustrip.(u"Gyr", data_dict[:snap_data].physical_time)
            sf = data_dict[:snap_data].scale_factor
            rs = data_dict[:snap_data].redshift

            # Write the time metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Time [Gyr]"]   = pt
            attrs(hdf5_group["snap_$(snapshot_number)"])["Scale factor"] = sf
            attrs(hdf5_group["snap_$(snapshot_number)"])["Redshift"]     = rs

            # Write the unit metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Mass unit"]     = string(m_unit)
            attrs(hdf5_group["snap_$(snapshot_number)"])["Length unit"]   = string(l_unit)
            attrs(hdf5_group["snap_$(snapshot_number)"])["Velocity unit"] = string(v_unit)

            # Write the grid metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Grid size [length unit]"] = ustrip.(
                l_unit,
                grid.physical_size,
            )
            attrs(hdf5_group["snap_$(snapshot_number)"])["Grid size [# voxels]"] = grid.n_bins

            # Write the column metadata
            attrs(hdf5_group["snap_$(snapshot_number)"])["Columns"] = [
                "x",  # Column 01: x coordinate [l_unit]
                "y",  # Column 02: y coordinate [l_unit]
                "z",  # Column 03: z coordinate [l_unit]
                "M*", # Column 04: Stellar mass [m_unit]
                "Vx", # Column 05: Velocity in the x direction [v_unit]
                "Vy", # Column 06: Velocity in the y direction [v_unit]
                "Vz", # Column 07: Velocity in the z direction [v_unit]
                "Sx", # Column 08: Velocity dispersion in the x direction [v_unit]
                "Sy", # Column 09: Velocity dispersion in the y direction [v_unit]
                "Sz", # Column 10: Velocity dispersion in the z direction [v_unit]
            ]

            next!(prog_bar)

        end

    end

    close(hdf5_file)

    return nothing

end

"""
    clumpingFactor(
        simulation_paths::Vector{String},
        slice::IndexType,
        quantity::Symbol;
        <keyword arguments>
    )::Nothing

Plot the clumping factor of `quantity` for different volume scales.

# Arguments

  - `simulation_paths::Vector{String}`: Paths to the simulation directories, set in the code variable `OutputDir`.
  - `slice::IndexType`: Slice of the simulations, i.e. which snapshots will be plotted. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots). Starts at 1 and out of bounds indices are ignored.
  - `quantity::Symbol`: The number density of which quantity will be used. The options are:

      + `:gas`          -> Gas number density.
      + `:molecular`    -> Molecular hydrogen number density.
      + `:br_molecular` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic`       -> Atomic hydrogen number density.
      + `:ionized`      -> Ionized hydrogen number density.
      + `:neutral`      -> Neutral hydrogen number density.
  - `nn::Int=32`: Number of neighbors.
  - `smooth::Int=0`: The result will be average out using `smooth` bins for the volume. Set it to 0 if you want no smoothing.
  - `x_trim::NTuple{2,<:Real}=(-Inf, Inf)`: The data will be trim down so the x coordinates fit within `x_trim`.
  - `output_path::String="./"`: Path to the output folder.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Consider every cell/particle within the simulation box.
      + `:halo`            -> Consider only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Consider only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Consider only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Consider only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
      + A dictionary with three entries:

          + `:filter_function` -> The filter function.
          + `:translation`     -> Translation for the simulation box. The posibilites are:

              + `:zero`                       -> No translation is applied.
              + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
              + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc, after filtering) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
              + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
              + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
              + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
          + `:rotation`        -> Rotation for the simulation box. The posibilites are:

              + `:zero`                       -> No rotation is applied.
              + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
              + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
              + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
              + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `da_ff::Function=filterNothing`: A function with the signature:

    `da_ff(data_dict) -> indices`

    where

      + `data_dict::Dict`: A dictionary with the following shape:

        * `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
        * `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
        * `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
        * ...
      + `indices::Dict`: A dictionary with the following shape:

        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * `cell/particle type` -> idxs::IndexType
        * ...
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `da_ff` filter function.
  - `sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths)`: Labels for the plot legend, one per simulation. Set it to `nothing` if you don't want a legend.
  - `theme::Attributes=Theme()`: Plot theme that will take precedence over [`DEFAULT_THEME`](@ref).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function clumpingFactor(
    simulation_paths::Vector{String},
    slice::IndexType,
    quantity::Symbol;
    nn::Int=32,
    smooth::Int=100,
    x_trim::NTuple{2,<:Real}=(-Inf, Inf),
    output_path::String="./",
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    da_ff::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    sim_labels::Union{Vector{<:AbstractString},Nothing}=basename.(simulation_paths),
    theme::Attributes=Theme(),
)::Nothing

    (
        quantity ∈ [:gas, :molecular, :br_molecular, :atomic, :ionized, :neutral] ||
        throw(ArgumentError("clumpingFactor: `quantity` can only be :gas, :molecular, \
        :br_molecular, :atomic, :ionized or :neutral, but I got :$(quantity)"))
    )

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(plotParams(Symbol(quantity, "_number_density")).request, ff_request),
    )

    plotSnapshot(
        simulation_paths,
        request,
        [scatter!];
        pf_kwargs=[(; markersize=10)],
        # `plotSnapshot` configuration
        output_path,
        base_filename="$(quantity)_clumping_factor",
        output_format=".png",
        show_progress=true,
        # Data manipulation options
        slice=slice,
        filter_function,
        da_functions=[daClumpingFactor],
        da_args=[(quantity,)],
        da_kwargs=[(; nn, filter_function=da_ff)],
        post_processing=getNothing,
        pp_args=(),
        pp_kwargs=(;),
        transform_box=true,
        translation,
        rotation,
        smooth,
        x_unit=u"kpc^3",
        y_unit=Unitful.NoUnits,
        x_exp_factor=0,
        y_exp_factor=0,
        x_trim,
        y_trim=(-Inf, Inf),
        x_edges=false,
        y_edges=false,
        x_func=identity,
        y_func=identity,
        # Axes options
        xaxis_label="auto_label",
        yaxis_label="auto_label",
        xaxis_var_name=L"\bar{V}",
        yaxis_var_name=L"C_\rho",
        xaxis_scale_func=identity,
        yaxis_scale_func=identity,
        # Plotting and animation options
        save_figures=true,
        backup_results=false,
        theme,
        sim_labels,
        title="",
        colorbar=false,
        # Animation options
        animation=false,
        animation_filename="animation.mp4",
        framerate=10,
    )

    return nothing

end
