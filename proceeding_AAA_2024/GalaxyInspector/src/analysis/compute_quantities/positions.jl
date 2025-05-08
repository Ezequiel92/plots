####################################################################################################
# Compute characteristic positions
####################################################################################################

"""
    computeCenter(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Length}

Read the position of the particle/cell at the potencial minimum of a given halo or subhalo.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
  - `subfind_idx::NTuple{2,Int}`: Tuple with two elements:

      + Index of the target halo (FoF group). Starts at 1.
      + Index of the target subhalo (subfind), relative the target halo. Starts at 1. If it is set to 0, the potencial minimum of the halo with index `halo_idx` is returned.

# Returns

  - The position of the potencial minimum.
"""
function computeCenter(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Length}

    # If there are no subfind data, return the origin
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return zeros(typeof(1.0u"kpc"), 3)
    end

    halo_idx, subhalo_rel_idx = subfind_idx

    # Load the necessary data
    n_subhalos_in_halo = data_dict[:group]["G_Nsubs"]
    g_pos = data_dict[:group]["G_Pos"]
    s_pos = data_dict[:subhalo]["S_Pos"]

    # Check that the requested halo index is within bounds
    n_halos = data_dict[:gc_data].header.n_groups_total

    (
        !iszero(n_halos) && !any(isempty, [n_subhalos_in_halo, g_pos, s_pos]) ||
        return zeros(typeof(1.0u"kpc"), 3)
    )

    (
        0 < halo_idx <= n_halos ||
        throw(ArgumentError("computeCenter: There is only $(n_halos) FoF groups in \
        $(data_dict[:gc_data].path), so halo_idx = $(halo_idx) is out of bounds"))
    )

    # Select the halo potencial minimum if `subhalo_rel_idx` == 0
    isPositive(subhalo_rel_idx) || return g_pos[:, halo_idx]

    # Check that the requested subhalo index is within bounds
    n_subfinds = n_subhalos_in_halo[halo_idx]

    if iszero(n_subfinds)

        (
            !logging[] ||
            @info("computeCenter: There are 0 subhalos in the FoF group $(halo_idx) from \
            $(data_dict[:gc_data].path), so the center will be the halo potencial minimum")

        )


        return g_pos[:, halo_idx]

    end

    (
        subhalo_rel_idx <= n_subfinds ||
        throw(ArgumentError("computeCenter: There is only $(n_subfinds) subhalos for the FoF \
        group $(halo_idx) in $(data_dict[:gc_data].path), so `subhalo_rel_idx` = \
        $(subhalo_rel_idx) is out of bounds"))
    )

    # Compute the number of subhalos and particles up to the last halo before `halo_idx`
    if isone(halo_idx)
        n_subs_floor = 0
    else
        n_subs_floor = sum(n_subhalos_in_halo[1:(halo_idx - 1)]; init=0)
    end

    # Compute the subhalo absolute index
    subhalo_abs_idx = n_subs_floor + subhalo_rel_idx

    # Select the subhalo potencial minimum
    return s_pos[:, subhalo_abs_idx]

end

"""
    computeCenter(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Length}

Read the position of the potencial minimum of a given subhalo.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
  - `subhalo_abs_idx::Int`: Absolute index of the target subhalo (subfind). Starts at 1.

# Returns

  - The specified potencial minimum.
"""
function computeCenter(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Length}

    # If there are no subfind data, return the origin
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return zeros(typeof(1.0u"kpc"), 3)
    end

    s_pos = data_dict[:subhalo]["S_Pos"]

    # Check that the requested subhalo index is within bounds
    n_subgroups_total = data_dict[:gc_data].header.n_subgroups_total

    !iszero(n_subgroups_total) && !isempty(s_pos) || return zeros(typeof(1.0u"kpc"), 3)

    (
        0 < subhalo_abs_idx <= n_subgroups_total ||
        throw(ArgumentError("computeCenter: There is only $(n_subgroups_total) subhalos in \
        $(data_dict[:gc_data].path), so `subhalo_abs_idx` = $(subhalo_abs_idx) is out of bounds"))
    )

    # Select the subhalo potencial minimum
    return s_pos[:, subhalo_abs_idx]

end

"""
    computeCenter(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Length}

Compute a characteristic center for the system.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
  - `cm_type::Symbol`: It can be:

      + `:global_cm`   -> Center of mass of the whole system.
      + `:{component}` -> Center of mass of the given component (e.g. :stars, :gas, :halo, etc). It can be any of the keys of [`PARTICLE_INDEX`](@ref).
      + `:zero`        -> Origin.

# Returns

  - The specified center of mass.
"""
function computeCenter(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Length}

    if cm_type == :global_cm

        return computeGlobalCenterOfMass(data_dict)

    elseif cm_type ∈ keys(PARTICLE_INDEX)

        return computeCenterOfMass(data_dict[cm_type]["POS "], data_dict[cm_type]["MASS"])

    elseif cm_type == :zero

        return zeros(typeof(1.0u"kpc"), 3)

    end

    throw(ArgumentError("computeCenter: `cm_type` can only be :global_cm, :zero or one of the keys \
    of `PARTICLE_INDEX` but I got :$(center_type)"))

end

"""
    computeDistance(
        positions::Matrix{<:Number};
        <keyword arguments>
    )::Vector{<:Number}

Compute the distance of a group of points to `center`.

# Arguments

  - `positions::Matrix{<:Number}`: Positions of the points. Each column is a point and each row a dimension.
  - `center::Union{Vector{<:Number},Nothing}=nothing`: Origin used to compute the distances. If set to `nothing`, 0 is used.

# Returns

  - The distance of every point to `center`.
"""
function computeDistance(
    positions::Matrix{<:Number};
    center::Union{Vector{<:Number},Nothing}=nothing,
)::Vector{<:Number}

    if isempty(positions)
        return eltype(positions)[]
    end

    if center === nothing
        return [norm(col) for col in eachcol(positions)]
    end

    (
        length(center) == size(positions, 1) ||
        throw(ArgumentError("computeDistance: `center` must have as many elements as `positions` \
        has rows, but I got length(`center`) = $(length(center)) and size(`positions`, 1) = \
        $(size(positions, 1))"))
    )

    return [norm(col .- center) for col in eachcol(positions)]

end

"""
    computeCenterOfMass(
        positions::Matrix{<:Unitful.Length},
        mass::Vector{<:Unitful.Mass},
    )::Vector{<:Unitful.Length}

Compute the center of mass of a group of cells/particles.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.

# Returns

  - The center of mass.
"""
function computeCenterOfMass(
    positions::Matrix{<:Unitful.Length},
    masses::Vector{<:Unitful.Mass},
)::Vector{<:Unitful.Length}

    # Check for missing data
    !any(isempty, [positions, masses]) || return zeros(typeof(1.0u"kpc"), 3)

    # Compute the total mass
    M = sum(masses)

    # Compute the center of mass
    center_of_mass = [sum(row .* masses) / M for row in eachrow(positions)]

    return center_of_mass

end

"""
    computeGlobalCenterOfMass(data_dict::Dict)::Vector{<:Unitful.Length}

Compute the center of mass of the whole system.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...

# Returns

  - The center of mass.
"""
function computeGlobalCenterOfMass(data_dict::Dict)::Vector{<:Unitful.Length}

    components = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), components)

    # Concatenate the position and masses of all the cells and particles in the system
    positions = hcat([data_dict[component]["POS "] for component in components]...)
    masses    = vcat([data_dict[component]["MASS"] for component in components]...)

    # Check for missing data
    !any(isempty, [positions, masses]) || return zeros(typeof(1.0u"kpc"), 3)

    (
        !logging[] ||
        @info("computeGlobalCenterOfMass: The center of mass will be computed using $(components)")
    )

    return computeCenterOfMass(positions, masses)

end

"""
    findHaloSubhalo(
        data_dict::Dict,
        star_idxs::Vector{Int},
        real_stars_idxs::Vector{Bool},
    )::NTuple{2,Vector{Int}}

Find in which halo and subhalo of `data_dict` each star in `star_idxs` was born.

For stars with no halo or subhalo, an index of -1 is given. The subhalo index is relative to the corresponding halo.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
  - `star_idxs::Vector{Int}`: Indices of the target stars in `data_dict`.
  - `real_stars_idxs::Vector{Bool}`: Boolean list of stellar particles. True for a real star and false for a wind particle.

# Returns

  - A tuple with two elements:

      + A vector with the birth halo (index starting at 1) of each star (in the order of `star_idxs`).
      + A vector with the birth subhalo (index starting at 1) of each star (in the order of `star_idxs`).
"""
function findHaloSubhalo(
    data_dict::Dict,
    star_idxs::Vector{Int},
    real_stars_idxs::Vector{Bool},
)::NTuple{2,Vector{Int}}

    ################################################################################################
    # Read the subfind metadata
    ################################################################################################

    # Read the total number of halos
    n_halos = data_dict[:gc_data].header.n_groups_total

    # Read the number of subhalos in each halo
    n_subhalos_in_halo = data_dict[:group]["G_Nsubs"]

    # Read the number of stars in each halo
    n_stars_in_halo = data_dict[:group]["G_LenType"][PARTICLE_INDEX[:stars] + 1, :]

    # Read the number of stars in each subhalo
    n_stars_in_subhalo = data_dict[:subhalo]["S_LenType"][PARTICLE_INDEX[:stars] + 1, :]

    ################################################################################################
    # Allocate memory
    ################################################################################################

    # If each halo was the birth place of a star form `star_idxs`
    # So as to compute the relevant indices of each halo only once
    born_in_this_halo = fill(false, n_halos)

    # Index of the last real star particle belonging to each subhalo
    # Each element of this vector corresponds to a halo
    last_idxs_in_subhalo_list = Vector{Vector{Int}}(undef, n_halos)

    for i in 1:n_halos

        # Number of subhalos in halo `i`
        n_subfinds = n_subhalos_in_halo[i]

        # Index of the last real star particle belonging to each subhalo
        # Each element of this vector corresponds to a subhalo of halo `i`
        last_idxs_in_subhalo_list[i] = Vector{Int}(undef, n_subfinds)

    end

    # Output vectors
    halo_idxs = fill(-1, length(star_idxs))
    subhalo_idxs = fill(-1, length(star_idxs))

    ################################################################################################
    # Compute the index of the last real star particle belonging to each halo
    ################################################################################################

    # Compute the index of the last star/wind particle in each of the halos
    last_idxs_in_halo = cumsum(n_stars_in_halo)

    for (i, idx) in enumerate(last_idxs_in_halo)

        # Compute the number of wind particles up to the particle with index `idx`
        n_wind = count(x -> !(x), real_stars_idxs[1:idx])

        # Shift `last_idxs_in_halo` to ignore wind particles
        last_idxs_in_halo[i] = idx - n_wind

    end

    ############################################################################################
    # Compute in which halo and subhalo each star was born
    ############################################################################################
    for (i, star_idx) in enumerate(star_idxs)

        ############################################################################################
        # Compute in which halo each star was born
        ############################################################################################

        # Find the halo where the target star was born
        halo_idx = searchsortedfirst(last_idxs_in_halo, star_idx)

        # If the star does not belong to any halo, leave the index as -1
        halo_idx <= length(last_idxs_in_halo) || continue

        halo_idxs[i] = halo_idx

        ############################################################################################
        # Compute in which subhalo each star was born
        ############################################################################################

        # Index of the last real star particle belonging to each subhalo
        # Each element of this vector corresponds to a subhalo of halo `halo_idx`
        last_idxs_in_subhalo = last_idxs_in_subhalo_list[halo_idx]

        # If it is the first time checking a star born in the halo `halo_idx`,
        # compute the index of the last real star particle in each subhalo of halo `halo_idx`
        if !born_in_this_halo[halo_idx]

            if isone(halo_idx)
                # Absolute index of the first subhalo
                first_subhalo_abs_idx = 1

                # Absolute index of the last subhalo
                last_subhalo_abs_idx = n_subhalos_in_halo[1]

                # Number of stars up to, but not including, the halo `halo_idx`
                n_star_floor = 0
            else
                # Absolute index of the first subhalo
                first_subhalo_abs_idx = sum(n_subhalos_in_halo[1:(halo_idx - 1)]) + 1

                # Absolute index of the last subhalo
                last_subhalo_abs_idx = sum(n_subhalos_in_halo[1:halo_idx])

                # Number of stars up to, but not including, the halo `halo_idx`
                n_star_floor = sum(n_stars_in_halo[1:(halo_idx - 1)])
            end

            # Compute the index of the last star/wind particle in each of the subhalos of the halo `halo_idx`
            cumsum!(last_idxs_in_subhalo, n_stars_in_subhalo[first_subhalo_abs_idx:last_subhalo_abs_idx])
            map!(x -> x + n_star_floor, last_idxs_in_subhalo, last_idxs_in_subhalo)

            # Compute the index of the last real star particle in each of the subhalos of the halo `halo_idx`
            for (i, idx) in enumerate(last_idxs_in_subhalo)

                # Compute the number of wind particles up to the particle with index `idx`
                n_wind = count(x -> !(x), real_stars_idxs[(n_star_floor + 1):idx])

                # Shift `last_idxs_in_subhalo` to ignore wind particles
                last_idxs_in_subhalo[i] = idx - n_wind

            end

            # Set to true to compute the relevant indices of this halo only once
            born_in_this_halo[halo_idx] = true

        end

        # Find the (relative index of the) subhalo where the target star was born
        subhalo_idx = searchsortedfirst(last_idxs_in_subhalo, star_idx)

        # If the star does not belong to any subhalo, leave the index as -1
        subhalo_idx <= length(last_idxs_in_subhalo) || continue

        subhalo_idxs[i] = subhalo_idx

    end

    return halo_idxs, subhalo_idxs

end

"""
    locateStellarBirthPlace(data_dict::Dict)::NTuple{2,Vector{Int}}

Find in which halo and subhalo each star in `data_dict` was born.

For stars with no halo or subhalo, an index of -1 is given. The subhalo index is relative to the corresponding halo.

# Arguments

  - `data_dict::Dict`: A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...

# Returns

  - A tuple with two elements:

      + A vector with the birth halo (index starting at 1) of each star (in the order of `data_dict`).
      + A vector with the birth subhalo (index starting at 1) of each star (in the order of `data_dict`).
"""
function locateStellarBirthPlace(data_dict::Dict)::NTuple{2,Vector{Int}}

    ################################################################################################
    # Read the data in `data_dict`
    ################################################################################################

    # Read the birth time of each star
    birth_ticks = data_dict[:stars]["GAGE"]

    if data_dict[:sim_data].cosmological
        # Go from scale factor to physical time
        birth_times = computeTime(birth_ticks, data_dict[:snap_data].header)
    else
        birth_times = birth_ticks
    end

    # Read the time stamp of each snapshot
    times = data_dict[:sim_data].table[!, :physical_times]

    (
        length(times) >= 2 ||
        throw(ArgumentError("locateStellarBirthPlace: I found less that two snapshots in \
        $(data_dict[:sim_data].path). But I need more to locate the birth place of the stars"))
    )

    # Read the ID of each star
    ids = data_dict[:stars]["ID  "]

    # Compute the number of snapshots
    n_snaps = length(times)

    ################################################################################################
    # Compute the indices and IDs of the stars born between each of the snapshots
    ################################################################################################

    # Allocate memory
    present_star_idxs = [Int[] for _ in 1:n_snaps]

    for (star_idx, birth_time) in enumerate(birth_times)

        snap_idx = searchsortedfirst(times, birth_time)

        if snap_idx > n_snaps
            push!(present_star_idxs[n_snaps], star_idx)
        else
            push!(present_star_idxs[snap_idx], star_idx)
        end

    end

    # Read the IDs of the stars born between each of the snapshots
    star_ids = [ids[idxs] for idxs in present_star_idxs]

    ################################################################################################
    # Read each snapshot and find the original halo and subhalos of the stars born there
    ################################################################################################

    # Make a dataframe with the following columns:
    #   - 1. DataFrame index
    #   - 2. Number in the file name
    #   - 3. Scale factor
    #   - 4. Redshift
    #   - 5. Physical time
    #   - 6. Lookback time
    #   - 7. Snapshot path
    #   - 8. Group catalog path
    simulation_table = makeSimulationTable(data_dict[:sim_data].path)

    # Allocate memory
    birth_halo    = fill(-1, length(birth_times))
    birth_subhalo = fill(-1, length(birth_times))

    request = Dict(
        :stars => ["ID  "],
        :group => ["G_Nsubs", "G_LenType"],
        :subhalo => ["S_LenType"],
    )

    for (global_idx, snapshot_row) in pairs(eachrow(simulation_table))

        # Select the IDs of the stars born in this snapshot
        ids = star_ids[global_idx]

        # Select the present index of the stars born in this snapshot
        present_idxs = present_star_idxs[global_idx]

        # Skip snapshots with no stars born in them
        !isempty(ids) || continue

        # Get the snapshot file path
        snapshot_path = snapshot_row[:snapshot_paths]

        # Get the group catalog file path
        groupcat_path = snapshot_row[:groupcat_paths]

        # Skip missing snapshots
        !ismissing(snapshot_path) || continue

        # Store the metadata of the current snapshot and simulation
        metadata = Dict(
            :sim_data => data_dict[:sim_data],
            :snap_data => Snapshot(
                snapshot_path,
                global_idx,
                global_idx,
                snapshot_row[:physical_times],
                snapshot_row[:lookback_times],
                snapshot_row[:scale_factors],
                snapshot_row[:redshifts],
                readSnapHeader(snapshot_path),
            ),
            :gc_data => GroupCatalog(
                groupcat_path,
                readGroupCatHeader(groupcat_path),
            ),
        )

        # Read the data in the snapshot
        past_data_dict = merge(
            metadata,
            readSnapshot(snapshot_path, request),
            readGroupCatalog(groupcat_path, snapshot_path, request),
        )

        # Get the birth index of the stars born in this snapshot
        past_idxs = parentIDToIndex(past_data_dict, ids)[:stars]

        (
            length(ids) == length(past_idxs) ||
            throw(DimensionMismatch("locateStellarBirthPlace:  There are IDs in `ids` that are not \
            present in the birth snapshot or are from other cell/particle type. \
            This should be impossible!"))
        )

        # Find the halo and subhalo where each star was born, for the stars born in this snapshot
        halo_idxs, subhalo_idxs = findHaloSubhalo(
            past_data_dict,
            past_idxs,
            findRealStars(past_data_dict[:snap_data].path),
        )

        # Store the halo and subhalos indices in the current position of each star
        birth_halo[present_idxs] .= halo_idxs
        birth_subhalo[present_idxs] .= subhalo_idxs

    end

    (
        allequal(length, [birth_halo, birth_subhalo, birth_ticks]) ||
        throw(ArgumentError("locateStellarBirthPlace: The results do not have as many elements as
        there are stars. Something went wrong!"))
    )

    return birth_halo, birth_subhalo

end
