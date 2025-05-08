####################################################################################################
# Filters
####################################################################################################

"""
    filterData!(data_dict::Dict; <keyword arguments>)::Nothing

Filter `data_dict` using the indices provided by `filter_function`.

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
"""
function filterData!(data_dict::Dict; filter_function::Function=filterNothing)::Nothing

    # Compute the filter dictionary
    indices = filter_function(data_dict)

    for component in snapshotTypes(data_dict)

        idxs = indices[component]

        for (block, values) in data_dict[component]
            if !isempty(values)
                data_dict[component][block] = collect(selectdim(values, ndims(values), idxs))
            end
        end

    end

    return nothing

end

"""
    filterData(data_dict::Dict; <keyword arguments>)::Dict

Return a filtered copy of `data_dict` using the indices provided by `filter_function`.

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

  - The filtered data.
"""
function filterData(data_dict::Dict; filter_function::Function=filterNothing)::Dict

    dd_copy = deepcopy(data_dict)

    # Compute the filter dictionary
    indices = filter_function(dd_copy)

    for component in snapshotTypes(dd_copy)

        idxs = indices[component]

        for (block, values) in dd_copy[component]
            if !isempty(values)
                dd_copy[component][block] = collect(selectdim(values, ndims(values), idxs))
            end
        end

    end

    return dd_copy

end

"""
    selectFilter(
        filter_mode::Symbol,
        request::Dict{Symbol,Vector{String}},
    )::Tuple{Function,Union{Symbol,NTuple{2,Int}},Symbol,Dict{Symbol,Vector{String}}}

Select a filter function, and the corresponding translation and rotation for the simulation box, from a list of premade ones.

Creates a request dictionary, using `request` as a base, adding what is necessary for the filter function and corresponding transformations.

# Arguments

  - `filter_mode::Symbol`: Which cells/particles will be plotted, the options are:

      + `:all`             -> Plot every cell/particle within the simulation box.
      + `:halo`            -> Plot only the cells/particles that belong to the main halo.
      + `:subhalo`         -> Plot only the cells/particles that belong to the main subhalo.
      + `:sphere`          -> Plot only the cell/particle inside a sphere with radius `DISK_R` (see `./src/constants/globals.jl`).
      + `:stellar_subhalo` -> Plot only the cells/particles that belong to the main subhalo.
      + `:all_subhalo`     -> Plot every cell/particle centered around the main subhalo.
  - `request::Dict{Symbol,Vector{String}}`: Base request dictionary, nothing will be deleted from it.

# Returns

  - A Tuple with four elements:

      + The filter function.
      + Translation for the simulation box. The posibilites are:

          + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
          + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
          + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new origin.
          + `(halo_idx, 0)`               -> Selects the center of mass of the `halo_idx::Int` halo, as the new origin.
      + Rotation for the simulation box. The posibilites are:

          + `:global_am`          -> Sets the angular momentum of the whole system as the new z axis.
          + `:stellar_am`         -> Sets the stellar angular momentum as the new z axis.
          + `:stellar_pa`         -> Sets the stellar principal axis as the new coordinate system.
          + `:stellar_subhalo_pa` -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
      + New request dictionary.
"""
function selectFilter(
    filter_mode::Symbol,
    request::Dict{Symbol,Vector{String}},
)::Tuple{Function,Union{Symbol,NTuple{2,Int}},Symbol,Dict{Symbol,Vector{String}}}

    if filter_mode == :all

        # Plot every cell/particle
        filter_function = filterNothing
        translation = :global_cm
        rotation = :global_am

        new_request = mergeRequests(
            addRequest(
                request,
                Dict(component => ["POS ", "MASS", "VEL "] for component in keys(PARTICLE_INDEX)),
            ),
            Dict(:stars => ["POS ", "MASS", "VEL ", "GAGE"]),
        )

    elseif filter_mode == :halo

        # Plot only the cells/particles that belong to the main halo
        filter_function = dd -> filterBySubhalo(dd; halo_idx=1, subhalo_rel_idx=0)
        translation = (1, 0)
        rotation = :stellar_pa

        new_request = mergeRequests(
            addRequest(
                request,
                Dict(
                    component => ["POS ", "MASS", "VEL "] for component in keys(PARTICLE_INDEX)
                ),
            ),
            Dict(
                :group   => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
                :subhalo => ["S_LenType", "S_Pos", "S_Vel"],
                :stars   => ["POS ", "MASS", "VEL ", "GAGE"],
            ),
        )

    elseif filter_mode == :subhalo

        # Plot only the cells/particles that belong to the main subhalo
        filter_function = dd -> filterBySubhalo(dd; halo_idx=1, subhalo_rel_idx=1)
        translation = (1, 1)
        rotation = :stellar_pa

        new_request = mergeRequests(
            addRequest(
                request,
                Dict(
                    component => ["POS ", "MASS"] for component in keys(PARTICLE_INDEX)
                ),
            ),
            Dict(
                :group   => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
                :subhalo => ["S_LenType", "S_Pos", "S_Vel"],
                :stars   => ["POS ", "MASS", "VEL ", "GAGE"],
            ),
        )

    elseif filter_mode == :sphere

        # Plot only the cell/particle inside a sphere with radius `DISK_R`
        filter_function = dd -> filterWithinSphere(dd, (0.0u"kpc", DISK_R), :global_cm)
        translation = :global_cm
        rotation = :global_am

        new_request = addRequest(
            request,
            Dict(component => ["POS ", "MASS", "VEL "] for component in keys(PARTICLE_INDEX)),
        )

    elseif filter_mode == :stellar_subhalo

        # Plot only the cells/particles that belong to the main subhalo
        filter_function = dd -> filterBySubhalo(dd; halo_idx=1, subhalo_rel_idx=1)
        translation = :stars
        rotation = :stellar_pa

        new_request = mergeRequests(
            mergeRequests(request, Dict(:stars => ["POS ", "MASS", "VEL ", "GAGE"])),
            Dict(:group => ["G_Nsubs", "G_LenType"], :subhalo => ["S_LenType"]),
        )

    elseif filter_mode == :all_subhalo

        # Plot every cell/particle centered around the main subhalo
        filter_function = filterNothing
        translation = (1, 1)
        rotation = :stellar_subhalo_pa

        new_request = mergeRequests(
            addRequest(
                request,
                Dict(
                    component => ["POS ", "MASS"] for component in keys(PARTICLE_INDEX)
                ),
            ),
            Dict(
                :group   => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
                :subhalo => ["S_LenType", "S_Pos", "S_Vel"],
                :stars   => ["POS ", "MASS", "VEL ", "GAGE"],
            ),
        )

    else

        throw(ArgumentError("selectFilter: `filter_mode` can only be :all, :halo, :subhalo, \
        :stellar_subhalo, :all_subhalo, or :sphere, but I got :$(filter_mode)"))

    end

    return filter_function, translation, rotation, new_request

end

"""
    selectFilter(
        filter_mode::Dict{Symbol,Any},
        request::Dict{Symbol,Vector{String}},
    )::Tuple{
        Function,
        Union{Symbol,NTuple{2,Int},Int},
        Union{Symbol,NTuple{2,Int},Int},
        Dict{Symbol,Vector{String}},
    }

Select the filter function, translation, and rotation from `filter_mode`.

Creates a request dictionary, using `request` as a base, adding what is necessary for the filter function and corresponding transformations.

# Arguments

  - `filter_mode::Dict{Symbol,Any}`: A dictionary with three entries:

      + `:filter_function` -> The filter function.
      + `:translation`     -> Translation for the simulation box. The posibilites are:

          + `:zero`                       -> No translation is applied.
          + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
          + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
          + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
          + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
          + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
      + `:rotation`        -> Rotation for the simulation box. The posibilites are:

          + `:zero`                       -> No rotation is applied.
          + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
          + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
          + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
          + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
          + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new coordinate system.
          + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo as the new coordinate system.
          + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `request::Dict{Symbol,Vector{String}}`: Base request dictionary, nothing will be deleted from it.

# Returns

  - A Tuple with four elements:

      + The filter function.
      + Translation for the simulation box. The posibilites are:

          + `:global_cm`                  -> Selects the center of mass of the whole system as the new origin.
          + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
          + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
          + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
          + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
      + Rotation for the simulation box. The posibilites are:

          + `:zero`                       -> No rotation is applied.
          + `:global_am`                  -> Sets the angular momentum of the whole system as the new z axis.
          + `:stellar_am`                 -> Sets the stellar angular momentum as the new z axis.
          + `:stellar_pa`                 -> Sets the stellar principal axis as the new coordinate system.
          + `:stellar_subhalo_pa`         -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
          + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new coordinate system.
          + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo as the new coordinate system.
          + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
      + New request dictionary.
"""
function selectFilter(
    filter_mode::Dict{Symbol,Any},
    request::Dict{Symbol,Vector{String}},
)::Tuple{
    Function,
    Union{Symbol,NTuple{2,Int},Int},
    Union{Symbol,NTuple{2,Int},Int},
    Dict{Symbol,Vector{String}},
}

    new_request = mergeRequests(
        addRequest(
            request,
            Dict(component => ["POS ", "MASS", "VEL "] for component in keys(PARTICLE_INDEX)),
        ),
        Dict(
            :group   => ["G_Nsubs", "G_LenType", "G_Pos", "G_Vel"],
            :subhalo => ["S_LenType", "S_Pos", "S_Vel"],
            :stars   => ["POS ", "MASS", "VEL ", "GAGE"],
        ),
    )

    return (
        filter_mode[:filter_function],
        filter_mode[:translation],
        filter_mode[:rotation],
        new_request,
    )

end

"""
    intersectFilters(filters::Dict{Symbol,IndexType}...)::Dict{Symbol,IndexType}

Generate the filter resulting from intersecting `filters` (AND in boolean logic) .

# Arguments

  - `filters::Dict{Symbol,IndexType}`: Filter, as a dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function intersectFilters(filters::Dict{Symbol,IndexType}...)::Dict{Symbol,IndexType}

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    (
        allequal(keys.(filters)) ||
        throw(ArgumentError("intersectFilters: The filters must have the same list of components \
        (their keys)"))
    )

    for component in keys(filters[1])

        indices[component] = intersect(
            filters[1][component],
            getindex.(filters[2:end], component)...,
        )

    end

    return indices

end

intersectFilters(filters::Dict{Symbol,IndexType})::Dict{Symbol,IndexType} = filters

####################################################################################################
#
# A filter function must take a data dictionary, and return a filter dictionary
#
# These functions are for the second argument of `filterData` in `./src/analysis/filters.jl`
#
# Expected signature:
#
#   filter_function(data_dict) -> indices
#
# where:
#
#   - `data_dict::Dict`: A dictionary with the following shape:
#
#      + `:sim_data`          -> ::Simulation (see `Simulation` in `./src/constants/globals.jl`).
#      + `:snap_data`         -> ::Snapshot (see `Snapshot` in `./src/constants/globals.jl`).
#      + `:gc_data`           -> ::GroupCatalog (see `GroupCatalog` in `./src/constants/globals.jl`).
#      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#      + ...
#      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
#      + ...
#   - indices::Dict{Symbol,IndexType}: A dictionary with the following shape:
#
#      + `cell/particle type` -> idxs::IndexType
#      + `cell/particle type` -> idxs::IndexType
#      + `cell/particle type` -> idxs::IndexType
#      + ...
#
####################################################################################################

"""
Default filter function that does not filter any cells/particles.
"""
filterNothing(x...; y...)::Dict{Symbol,IndexType} = PASS_ALL

"""
    filterWithinSphere(
        data_dict::Dict,
        range::NTuple{2,<:Unitful.Length},
        origin...,
    )::Dict{Symbol,IndexType}

Filter out the cell/particles outside a given spherical shell.

# Arguments

  - `data::Dict`: A dictionary with the following shape:

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
  - `range::NTuple{2,<:Unitful.Length}`: Internal and external radius of the spherical shell.
  - `origin`: It can be any number and type of argument compatible with the second to last arguments of a [`computeCenter`](@ref) method.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterWithinSphere(
    data_dict::Dict,
    range::NTuple{2,<:Unitful.Length},
    origin...,
)::Dict{Symbol,IndexType}

    center = computeCenter(data_dict, origin...)

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    for component in snapshotTypes(data_dict)

        positions = data_dict[component]["POS "]

        if isempty(positions)
            indices[component] = (:)
        else
            distances = computeDistance(positions; center)
            indices[component] = map(x -> range[1] < x <= range[2], distances)
        end

    end

    return indices

end

"""
    filterWithinCylinder(
        data_dict::Dict,
        max_r::Unitful.Length,
        max_z::Unitful.Length,
        origin...,
    )::Dict{Symbol,IndexType}

Filter out the cell/particles outside a given cylinder.

# Arguments

  - `data::Dict`: A dictionary with the following shape:

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
  - `max_r::Unitful.Length`: Radius of the cylinder.
  - `max_z::Unitful.Length`: Half height of the cylinder.
  - `origin`: It can be any number and type of argument compatible with the second to last arguments of a [`computeCenter`](@ref) method.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterWithinCylinder(
    data_dict::Dict,
    max_r::Unitful.Length,
    max_z::Unitful.Length,
    origin...,
)::Dict{Symbol,IndexType}

    center  = computeCenter(data_dict, origin...)

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    for component in snapshotTypes(data_dict)

        positions = data_dict[component]["POS "]

        if isempty(positions)
            indices[component] = (:)
        else
            distances = computeDistance(positions; center)
            heights   = abs.(positions[3, :])
            indices[component] = map(r -> r <= max_r, distances) ∩ map(z -> z <= max_z, heights)
        end

    end

    return indices

end

"""
    filterBySubhalo(
        data_dict::Dict;
        <keyword arguments>
    )::Dict{Symbol,IndexType}

Filter out cells/particles that do not belong to a given halo and subhalo.

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
  - `halo_idx::Int`: Index of the target halo (FoF group). Starts at 1.
  - `subhalo_rel_idx::Int`: Index of the target subhalo (subfind), relative to the target halo. Starts at 1. If it is set to 0, all subhalos of the target halo are included.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterBySubhalo(
    data_dict::Dict;
    halo_idx::Int=1,
    subhalo_rel_idx::Int=1,
)::Dict{Symbol,IndexType}

    # If there are no subfind data, filter out every cell/particle
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return PASS_NONE
    end

    # Load the necessary data
    g_n_subs = data_dict[:group]["G_Nsubs"]
    g_len_type = data_dict[:group]["G_LenType"]
    s_len_type = data_dict[:subhalo]["S_LenType"]

    # If any of the data is missing return an empty filter dictionary
    n_groups_total = data_dict[:gc_data].header.n_groups_total

    !iszero(n_groups_total) && !any(isempty, [g_n_subs, g_len_type, s_len_type]) || return PASS_NONE

    # Check that the requested halo index is within bounds
    (
        0 < halo_idx <= n_groups_total ||
        throw(ArgumentError("filterBySubhalo: There is only $(n_groups_total) FoF groups in \
        $(data_dict[:gc_data].path), so `halo_idx` = $(halo_idx) is out of bounds"))
    )

    # Compute the number of subhalos and particles up to the last halo before `halo_idx`
    if isone(halo_idx)
        n_subs_floor = 0
        len_type_floor = zeros(Int32, size(s_len_type, 1))
    else
        n_subs_floor = sum(g_n_subs[1:(halo_idx - 1)]; init=0)
        len_type_floor = sum(g_len_type[:, 1:(halo_idx - 1)], dims=2; init=0)
    end

    # Check that the requested subhalo index is within bounds
    n_subfinds = g_n_subs[halo_idx]

    if iszero(n_subfinds)

        (
            !logging[] ||
            @info("filterBySubhalo: There are 0 subhalos in the FoF group $(halo_idx) from
            $(data_dict[:gc_data].path), so every particle will be filtered out")
        )

        return PASS_NONE

    end

    (
        subhalo_rel_idx <= n_subfinds ||
        throw(ArgumentError("filterBySubhalo: There is only $(n_subfinds) subhalos for the \
        FoF group $(halo_idx) from $(data_dict[:gc_data].path), so `subhalo_rel_idx` = \
        $(subhalo_rel_idx) is out of bounds"))
    )

    if subhalo_rel_idx <= 0

        # Consider all subhalos within the target halo
        first_idxs = len_type_floor .+ 1
        last_idxs  = first_idxs .+ g_len_type[:, halo_idx] .- 1

    else

        # Compute the subhalo absolute index
        subhalo_abs_idx = n_subs_floor + subhalo_rel_idx

        # Compute the number of particles in the current halo,
        # upto the last subhalo before `subhalo_rel_idx`
        if isone(subhalo_abs_idx)
            len_type_floor_in_halo = zeros(Int, size(s_len_type, 1))
        else
            len_type_floor_in_halo = sum(
                s_len_type[:, (n_subs_floor + 1):(subhalo_abs_idx - 1)], dims=2; init=0,
            )
        end

        # Compute the first and last index of the selected
        # cells/particles (for each cell/particle type)
        first_idxs = len_type_floor .+ len_type_floor_in_halo .+ 1
        last_idxs  = first_idxs .+ s_len_type[:, subhalo_abs_idx] .- 1

    end

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    # Find which cell/particle types are part of the keys of `data_dict`
    components_in_dd = snapshotTypes(data_dict)

    # Fill the filter dictionary
    for (i, (first_idx, last_idx)) in enumerate(zip(first_idxs, last_idxs))

        component = INDEX_PARTICLE[i - 1]

        # Only compute the indices for components in `data_dict`
        if component ∈ components_in_dd

            if first_idx == last_idx || iszero(last_idx)
                indices[component] = Int[]
            end

            if component == :stars

                # Find the indices of the stars, excluding wind particles
                real_stars_idxs = findRealStars(data_dict[:snap_data].path)

                n_wind_before = count(x -> !(x), real_stars_idxs[1:(first_idx - 1)])
                n_wind_between = count(x -> !(x), real_stars_idxs[first_idx:last_idx])

                stars_first_idx = first_idx - n_wind_before
                stars_last_idx = last_idx - n_wind_before - n_wind_between

                indices[component] = stars_first_idx:stars_last_idx

            else

                indices[component] = first_idx:last_idx

            end

        end

    end

    return indices

end

"""
    filterBySubhalo(data_dict::Dict, subhalo_abs_idx::Int)::Dict{Symbol,IndexType}

Filter out cells/particles that do not belong to a given subhalo.

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
  - `subhalo_abs_idx::Int`: Index of the target subhalo (subfind). Starts at 1.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterBySubhalo(data_dict::Dict, subhalo_abs_idx::Int)::Dict{Symbol,IndexType}

    # If there are no subfind data, filter out every cell/particle
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return PASS_NONE
    end

    # Load the necessary data
    s_len_type = data_dict[:subhalo]["S_LenType"]

    # If any of the data is missing return an empty filter dictionary
    n_subgroups_total = data_dict[:gc_data].header.n_subgroups_total

    !iszero(n_subgroups_total) && !isempty(s_len_type) || return PASS_NONE

    # Check that the requested subhalo index is within bounds
    (
        0 < subhalo_abs_idx <= n_subgroups_total ||
        throw(ArgumentError("filterBySubhalo: There is only $(n_subgroups_total) subhalos in \
        $(data_dict[:gc_data].path), so subhalo_abs_idx = $(subhalo_abs_idx) is out of bounds"))
    )

    # Compute the number of particles upto the last subhalo before `subhalo_abs_idx`
    if isone(subhalo_abs_idx)
        len_type_floor = zeros(Int, size(s_len_type, 1))
    else
        len_type_floor = sum(s_len_type[:, 1:(subhalo_abs_idx - 1)], dims=2; init=0)
    end

    # Compute the first and last index of the selected cells/particles (for each cell/particle type)
    first_idxs = len_type_floor .+ 1
    last_idxs  = len_type_floor .+ s_len_type[:, subhalo_abs_idx]

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    # Find which cell/particle types are part of the keys of `data_dict`
    components_in_dd = snapshotTypes(data_dict)

    # Fill the filter dictionary
    for (i, (first_idx, last_idx)) in enumerate(zip(first_idxs, last_idxs))

        component = INDEX_PARTICLE[i - 1]

        # Only compute the indices for components in `data_dict`
        if component ∈ components_in_dd

            if first_idx == last_idx || iszero(last_idx)
                indices[component] = Int[]
            end

            if component == :stars

                # Find the indices of the stars, excluding wind particles
                real_stars_idxs = findRealStars(data_dict[:snap_data].path)

                n_wind_before = count(x -> !(x), real_stars_idxs[1:(first_idx - 1)])
                n_wind_between = count(x -> !(x), real_stars_idxs[first_idx:last_idx])

                stars_first_idx = first_idx - n_wind_before
                stars_last_idx = last_idx - n_wind_before - n_wind_between

                indices[component] = stars_first_idx:stars_last_idx

            else

                indices[component] = first_idx:last_idx

            end

        end

    end

    return indices

end

"""
    filterByQuantity(
        data_dict::Dict,
        quantity::Symbol,
        component::Symbol,
        minimum::Number,
        maximum::Number,
    )::Dict{Symbol,IndexType}

Filter out particles/cells with `quantity` outside the range [`minimum:`, `maximum`].

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
  - `quantity::Symbol`: Target quantity. For the possibilities see the documentation of [`scatterQty`](@ref).
  - `component::Symbol`: Type of particle/cell. The possibilities are the keys of [`PARTICLE_INDEX`](@ref).
  - `minimum::Number`: Minimum value of `quantity`.
  - `maximum::Number`: Maximum value of `quantity`.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterByQuantity(
    data_dict::Dict,
    quantity::Symbol,
    component::Symbol,
    minimum::Number,
    maximum::Number,
)::Dict{Symbol,IndexType}

    # Compute the `quantity`
    values = scatterQty(data_dict, quantity)

    (
        minimum >= maximum &&
        throw(ArgumentError("filterByQuantity: `maximum` should be larger than `minimum`, \
        but I got `minimum` = $(minimum) >= `maximum` = $(maximum)"))
    )

    if isempty(values)

        (
            logging[] &&
            @warn("filterByQuantity: Quantity :$(quantity) is empty, no particles/cells will be \
            filtered out")
        )

        return PASS_ALL

    end

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    for type in snapshotTypes(data_dict)

        (
            length(data_dict[type]["MASS"]) != length(values) && logging[] &&
            @warn("filterByQuantity: The amount of particles/cells of type :$(component) is not \
            the same as the amount of values for quantity :$(quantity). Are you sure that \
            `component` = :$(component) is correct?")
        )

        if type == component
            indices[type] = map(x -> minimum <= x <= maximum, values)
        else
            indices[type] = (:)
        end

    end

    return indices

end

"""
    filterByEquilibrium(
        data_dict::Dict;
        <keyword arguments>
    )::Dict{Symbol,IndexType}

Filter out gas cells that have the molecular or ionized equation in or out of equilibrium, according to `equation` and `filtered_phase`.

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
  - `limit::Float64=1.0e-2`: Allowed deviation from equilibrium, as the fraction |RS - LS| / LS.
  - `equation::Symbol=:molecular`: Which equilibrium equation will be used. The options are :molecular and :ionized.
  - `filtered_phase::Symbol=:non_eq`: Which phase will be filtered out, the equilibrium phase (:eq) or the non-equilibrium phase (:non_eq).

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterByEquilibrium(
    data_dict::Dict;
    limit::Float64=1.0e-2,
    equation::Symbol=:molecular,
    filtered_phase::Symbol=:non_eq,
)::Dict{Symbol,IndexType}

    if equation == :molecular

        eq_quotient = GalaxyInspector.scatterQty(data_dict, :mol_eq_quotient)

    elseif equation == :ionized

        eq_quotient = GalaxyInspector.scatterQty(data_dict, :ion_eq_quotient)

    else

        throw(ArgumentError("filterByEquilibrium: `equation` can only be :molecular or :ionized, \
        but I got :$(equation)"))

    end

    (
        0.0 <= limit ||
        throw(ArgumentError("filterByEquilibrium: `limit` must be > 0, but I got \
        :$(limit_percent)"))
    )

    # Compute |RS - LS| / LS
    eq_values = abs.(exp10.(-eq_quotient) .- 1.0)

    if filtered_phase == :non_eq

        idxs = map(x -> x < limit, eq_values)

    elseif filtered_phase == :eq

        idxs = map(x -> limit < x, eq_values)

    else

        throw(ArgumentError("filterByEquilibrium: `filtered_phase` can only be :non_eq or :eq, \
        but I got :$(filtered_phase)"))

    end

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    for component in snapshotTypes(data_dict)

        if component == :gas
            indices[component] = idxs
        else
            indices[component] = (:)
        end

    end

    return indices

end

"""
    filterByELSFR(data_dict::Dict)::Dict{Symbol,IndexType}

Filter out gas cells that have not entered our star formation routine.

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

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterByELSFR(data_dict::Dict)::Dict{Symbol,IndexType}

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    for component in snapshotTypes(data_dict)

        if component == :gas
            if isempty(data_dict[:gas]["FRAC"])
                indices[component] = Int[]
            else
                indices[component] = map(!isnan, data_dict[:gas]["FRAC"][1, :])
            end
        else
            indices[component] = (:)
        end

    end

    return indices

end

"""
    filterOldStars(data_dict::Dict)::Dict{Symbol,IndexType}

Filter out stars that where born one or more snapshots ago.

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

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterOldStars(data_dict::Dict)::Dict{Symbol,IndexType}

    birth_ticks = data_dict[:stars]["GAGE"]

    # Get the global index (index in the context of the whole simulation) of the current snapshot
    present_idx = data_dict[:snap_data].global_index

    if present_idx == 1 || isempty(birth_ticks)

        new_stars_idxs = (:)

    else

        # Compute the stellar birth dates
        if data_dict[:sim_data].cosmological
            # Go from scale factor to physical time
            birth_times = computeTime(birth_ticks, data_dict[:snap_data].header)
        else
            birth_times = birth_ticks
        end

        # Get the physical times
        times = data_dict[:sim_data].table[:, 5]

        new_stars_idxs = map(t -> t > times[present_idx - 1], birth_times)

    end

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    for component in snapshotTypes(data_dict)

        if component == :stars
            indices[component] = new_stars_idxs
        else
            indices[component] = (:)
        end

    end

    return indices

end

"""
    filterByStellarAge(
        data_dict::Dict;
        <keyword arguments>
    )::Dict{Symbol,IndexType}

Filter out stars with an age outside the range [`min_age`, `max_age`].

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
  - `min_age::Unitful.Time=0.0u"Gyr"`: Minimum age.
  - `max_age::Unitful.Time=AGE_RESOLUTION`: Maximum age.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterByStellarAge(
    data_dict::Dict;
    min_age::Unitful.Time=0.0u"Gyr",
    max_age::Unitful.Time=AGE_RESOLUTION,
)::Dict{Symbol,IndexType}

    ages = computeStellarAge(data_dict)

    isempty(ages) && PASS_ALL

    new_stars_idxs = map(t -> min_age <= t <= max_age, ages)

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    for component in snapshotTypes(data_dict)

        if component == :stars
            indices[component] = new_stars_idxs
        else
            indices[component] = (:)
        end

    end

    return indices

end

"""
    filterByBirthPlace(
        data_dict::Dict,
        exclude::Symbol;
        <keyword arguments>
    )::Dict{Symbol,IndexType}

Filter out stars that where born either outside the given halo and subhalo (`exclude`= :exsitu), or inside (`exclude`= :insitu).

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
  - ´exclude::Symbol´: Which stars will be filtered out, either the ones born outside the given halo and subhalo (:exsitu), or inside (:insitu).
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.
  - `subhalo_rel_idx::Int=1`: Index of the target subhalo (subfind), relative to the target halo. Starts at 1. If it is set to 0, all subhalos of the target halo are consider insitu.

# Returns

  - A dictionary with the following shape:

      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + `cell/particle type` -> idxs::IndexType
      + ...
"""
function filterByBirthPlace(
    data_dict::Dict,
    exclude::Symbol;
    halo_idx::Int=1,
    subhalo_rel_idx::Int=1,
)::Dict{Symbol,IndexType}

    birth_halo, birth_subhalo = locateStellarBirthPlace(data_dict)

    stars_born_in_halo = map(isequal(halo_idx), birth_halo)

    if iszero(subhalo_rel_idx)
        stars_born_in_subhalo = (:)
    else
        stars_born_in_subhalo = map(isequal(subhalo_rel_idx), birth_subhalo)
    end

    insitu_stars = stars_born_in_halo ∩ stars_born_in_subhalo

    if exclude == :insitu
        stars_idxs = Vector{Bool}(.!(insitu_stars))
    elseif exclude == :exsitu
        stars_idxs = insitu_stars
    else
        throw(ArgumentError("filterByBirthPlace: `exclude` can only be :insitu or :exsitu, \
        but I got :$(exclude)"))
    end

    # Allocate memory
    indices = Dict{Symbol,IndexType}()

    for component in snapshotTypes(data_dict)

        if component == :stars
            indices[component] = stars_idxs
        else
            indices[component] = (:)
        end

    end

    return indices

end
