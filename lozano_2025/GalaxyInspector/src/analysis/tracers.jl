####################################################################################################
# Tracer functions
####################################################################################################

"""
    parentToTracerID(data_dict::Dict, target_ids::Vector{UInt})::Vector{UInt}

Find the IDs of the tracers corresponding to a given list of parent IDs.

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
  - `target_ids::Vector{UInt}`: List of parent IDs.

# Returns

  - A vector with the IDs of the tracers.
"""
function parentToTracerID(data_dict::Dict, target_ids::Vector{UInt})::Vector{UInt}

    # Read the full list of parent and tracer IDs
    parent_ids = data_dict[:tracer]["PAID"]
    tracer_ids = data_dict[:tracer]["TRID"]

    # Find the indices of the target IDs in the parent ID list, ignoring invalid targets
    idxs = filter!(!isnothing, indexin(target_ids, parent_ids))

    return tracer_ids[idxs]

end

"""
    tracerToParentID(data_dict::Dict, target_ids::Vector{UInt})::Vector{UInt}

Find the IDs of the parents corresponding to a given list of tracer IDs.

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
  - `target_ids::Vector{UInt}`: List of tracer IDs.

# Returns

  - A vector with the IDs of the parents.
"""
function tracerToParentID(data_dict::Dict, target_ids::Vector{UInt})::Vector{UInt}

    # Read the full list of parent and tracer IDs
    parent_ids = data_dict[:tracer]["PAID"]
    tracer_ids = data_dict[:tracer]["TRID"]

    # Find the indices of the target IDs in the tracer ID list, ignoring invalid targets
    idxs = filter!(!isnothing, indexin(target_ids, tracer_ids))

    return parent_ids[idxs]

end

"""
    parentIDToIndex(
        data_dict::Dict,
        target_ids::Vector{UInt},
    )::Dict{Symbol,Vector{Int}}

Find the indices of the cell/particles with IDs given by `target_ids`.

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
  - `target_ids::Vector{UInt}`: List of tracer IDs.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::Int
      + `cell/particle type` -> idxs::Int
      + `cell/particle type` -> idxs::Int
      + ...
"""
function parentIDToIndex(
    data_dict::Dict,
    target_ids::Vector{UInt},
)::Dict{Symbol,Vector{Int}}

    # Check which cell/particle types are present
    components = filter!(
        ts -> !isempty(data_dict[ts]["ID  "]),
        snapshotTypes(data_dict) ∩ [:stars, :gas, :black_hole],
    )

    # Allocate memory
    index_dict = Dict{Symbol,Vector{Int}}()

    for ts in components

        # Read the IDs of the cell/particles
        parent_ids = data_dict[ts]["ID  "]

        # Find the indices of the target IDs in the cell/particle ID list, ignoring invalid targets
        index_dict[ts] = filter!(!isnothing, indexin(target_ids, parent_ids))

    end

    return index_dict

end

"""
    findTracers(data_dict::Dict; <keyword arguments>)::Vector{UInt}

Find the tracers whose parents are allowed by `filter_function`.

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
  - `filter_function::Function=filterNothing`: A function with the signature:

    `filter_function(data_dict) -> indices`

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

# Returns

  - A vector with the IDs of the tracers.
"""
function findTracers(data_dict::Dict; filter_function::Function=filterNothing)::Vector{UInt}

    # Check which cell/particle types are present
    components = filter!(
        ts -> !isempty(data_dict[ts]["ID  "]),
        snapshotTypes(data_dict) ∩ [:stars, :gas, :black_hole],
    )

    # Find the indices of the cells and particles that are allowed by `filter_function`
    filter_idxs = filter_function(data_dict)

    # Read the IDs of the cells and particles that are allowed by `filter_function`
    parent_ids = vcat([data_dict[ts]["ID  "][filter_idxs[ts]] for ts in components]...)

    !isempty(parent_ids) || return Vector{UInt}[]

    # Find the IDs of the tracers whose parent are allowed by `filter_function`,
    # ignoring cells and particles with no tracers
    return parentToTracerID(data_dict, parent_ids)

end

"""
    tracersWithinR200(data_dict::Dict; <keyword arguments>)::Vector{UInt}

Find the tracers that are within the virial radius (``R_{200}``).

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
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.

# Returns

  - A vector with the IDs of the tracers.
"""
function tracersWithinR200(data_dict::Dict; halo_idx::Int=1)::Vector{UInt}

    if isempty(data_dict[:group]["G_R_Crit200"])

        filter_function = filterNothing

    else

        # Read the virial radius
        r200 = data_dict[:group]["G_R_Crit200"][halo_idx]

        # Construct a filter function that only allows cells and particles within the virial radius
        filter_function = dd -> filterWithinSphere(dd, (0.0u"kpc", r200), :zero)

    end

    return findTracers(data_dict; filter_function)

end

"""
    tracersWithinDisc(data_dict::Dict; <keyword arguments>)::Vector{UInt}

Find the tracers that are within a given cylinder.

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
  - `max_r::Unitful.Length=DISK_R`: Radius of the cylinder.
  - `max_z::Unitful.Length=5.0u"kpc"`: Half height of the cylinder.

# Returns

  - A vector with the IDs of the tracers.
"""
function tracersWithinDisc(
    data_dict::Dict;
    max_r::Unitful.Length=DISK_R,
    max_z::Unitful.Length=5.0u"kpc",
)::Vector{UInt}

    # Construct a filter function that only allows cells and particles within a given cylinder
    filter_function = dd -> filterWithinCylinder(dd, max_r, max_z, :zero)

    return findTracers(data_dict; filter_function)

end

"""
    tracersToMass(
        data_dict::Dict,
        target_ids::Vector{UInt},
    )::Dict{Symbol,Vector{<:Unitful.Mass}}

Find the masses of the parent cell/particles of the tracers with IDs `target_ids`.

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
  - `target_ids::Vector{UInt}`: List of tracer IDs.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> masses::Vector{<:Unitful.Mass}
      + `cell/particle type` -> masses::Vector{<:Unitful.Mass}
      + `cell/particle type` -> masses::Vector{<:Unitful.Mass}
      + ...
"""
function tracersToMass(
    data_dict::Dict,
    target_ids::Vector{UInt},
)::Dict{Symbol,Vector{<:Unitful.Mass}}

    # Find the parent ID corresponding to each tracer
    parent_ids = tracerToParentID(data_dict, target_ids)

    # Find the index and cell/particle type of each parent
    index_dict = parentIDToIndex(data_dict, parent_ids)

    return Dict(ts => data_dict[ts]["MASS"][idx] for (ts, idx) in index_dict)

end
