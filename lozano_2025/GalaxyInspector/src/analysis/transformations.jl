####################################################################################################
# Coordinate transformations
####################################################################################################

"""
    translatePoints(
        positions::Matrix{<:Number},
        new_origin::Vector{<:Number},
    )::Matrix{<:Number}

Translate a system of points, moving `new_origin` to the origin.

# Arguments

  - `positions::Matrix{<:Number}`: Points to be translated. Each column is a point and each row a dimension.
  - `new_origin::Vector{<:Number}`: Target origin.

# Returns

  - Matrix with the translated points.
"""
function translatePoints(
    positions::Matrix{<:Number},
    new_origin::Vector{<:Number},
)::Matrix{<:Number}

    !all(iszero, new_origin) || return positions

    return positions .- new_origin

end

"""
    translateData!(data_dict::Dict, translation::Union{Symbol,NTuple{2,Int},Int})::Nothing

Translate the positions of the cells/particles in `data_dict`.

!!! note

    The velocities will be boosted to the stellar center of mass of the system. If there are no stars, no transformation in applied to the velocities.

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
  - `translation::Union{Symbol,NTuple{2,Int},Int}=:zero`: Type of translation. The options are:

      + `:zero`                       -> No translation is applied.
      + `:global_cm`                  -> Sets the center of mass of the whole system as the new origin.
      + `:{component}`                -> Sets the center of mass of the given component (e.g. :stars, :gas, :halo, etc) as the new origin. It can be any of the keys of [`PARTICLE_INDEX`](@ref).
      + `(halo_idx, subhalo_rel_idx)` -> Sets the position of the potencial minimum for the `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new origin.
      + `(halo_idx, 0)`               -> Sets the center of mass of the `halo_idx::Int` halo as the new origin.
      + `subhalo_abs_idx`             -> Sets the center of mass of the `subhalo_abs_idx::Int` as the new origin.
"""
function translateData!(data_dict::Dict, translation::Union{Symbol,NTuple{2,Int},Int})::Nothing

    translation != :zero || return nothing

    new_origin = computeCenter(data_dict, translation)
    stellar_vcm = computeVcm(data_dict, translation)

    for component in snapshotTypes(data_dict)

        for (block, values) in data_dict[component]

            if !isempty(values)
                if block == "POS "
                    data_dict[component]["POS "] = translatePoints(values, new_origin)
                elseif block == "VEL "
                    data_dict[component]["VEL "] = translatePoints(values, stellar_vcm)
                end
            end

        end

    end

    return nothing

end

"""
    rotateData!(data_dict::Dict, axis_type::Symbol)::Nothing

Rotate the positions and velocities of the cells/particles in `data_dict`.

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
  - `rotation::Symbol`: Type of rotation. The options are:

      + `:zero`               -> No rotation is applied.
      + `:global_am`          -> Sets the angular momentum of the whole system as the new z axis.
      + `:stellar_am`         -> Sets the stellar angular momentum as the new z axis.
      + `:stellar_pa`         -> Sets the stellar principal axis as the new coordinate system.
      + `:stellar_subhalo_pa` -> Sets the principal axis of the stars in the main subhalo as the new coordinate system.
"""
function rotateData!(data_dict::Dict, rotation::Symbol)::Nothing

    # Compute the rotation matrix
    if rotation == :zero

        return nothing

    elseif rotation == :global_am

        rotation_matrix = computeGlobalAMRotationMatrix(data_dict)

    elseif rotation == :stellar_am

        !isempty(data_dict[:stars]["MASS"]) || return nothing

        rotation_matrix = computeAMRotationMatrix(
            data_dict[:stars]["POS "],
            data_dict[:stars]["VEL "],
            data_dict[:stars]["MASS"],
        )

    elseif rotation == :stellar_pa

        !isempty(data_dict[:stars]["MASS"]) || return nothing

        rotation_matrix = computePARotationMatrix(
            data_dict[:stars]["POS "],
            data_dict[:stars]["VEL "],
            data_dict[:stars]["MASS"],
        )

    elseif rotation == :stellar_subhalo_pa

        star_data = deepcopy(data_dict)

        filterData!(
            star_data,
            filter_function=dd -> filterBySubhalo(dd; halo_idx=1, subhalo_rel_idx=1),
        )

        !isempty(star_data[:stars]["MASS"]) || return nothing

        rotation_matrix = computePARotationMatrix(
            star_data[:stars]["POS "],
            star_data[:stars]["VEL "],
            star_data[:stars]["MASS"],
        )

    else

        throw(ArgumentError("rotateData!: I don't recognize the rotation :$(rotation)"))

    end

    for component in snapshotTypes(data_dict)

        for (block, values) in data_dict[component]

            if block ∈ ["POS ", "VEL "] && !isempty(values)
                data_dict[component][block] = rotation_matrix * values
            end

        end

    end

    return nothing

end

"""
    rotateData!(data_dict::Dict, axis_type::NTuple{2,Int})::Nothing

Rotate the positions and velocities of the cells/particles in `data_dict`.

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
  - `rotation::NTuple{2,Int}`: Type of rotation. The options are:

      + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo) as the new coordinate system.
      + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo as the new coordinate system.
"""
function rotateData!(data_dict::Dict, rotation::NTuple{2,Int})::Nothing

    # Compute the rotation matrix
    star_data = deepcopy(data_dict)

    filterData!(
        star_data,
        filter_function=dd -> filterBySubhalo(dd; halo_idx=rotation[1], subhalo_rel_idx=rotation[2]),
    )

    !isempty(star_data[:stars]["MASS"]) || return nothing

    rotation_matrix = computePARotationMatrix(
        star_data[:stars]["POS "],
        star_data[:stars]["VEL "],
        star_data[:stars]["MASS"],
    )

    for component in snapshotTypes(data_dict)

        for (block, values) in data_dict[component]

            if block ∈ ["POS ", "VEL "] && !isempty(values)
                data_dict[component][block] = rotation_matrix * values
            end

        end

    end

    return nothing

end

"""
    rotateData!(data_dict::Dict, axis_type::Int)::Nothing

Rotate the positions and velocities of the cells/particles in `data_dict`.

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
  - `rotation::Int`: Target subhalo absolute index, starting at 1. Sets the principal axis of the stars in the subhalo as the new coordinate system.
"""
function rotateData!(data_dict::Dict, rotation::Int)::Nothing

    # Compute the rotation matrix
    star_data = deepcopy(data_dict)

    filterData!(
        star_data,
        filter_function=dd -> filterBySubhalo(dd, rotation),
    )

    !isempty(star_data[:stars]["MASS"]) || return nothing

    rotation_matrix = computePARotationMatrix(
        star_data[:stars]["POS "],
        star_data[:stars]["VEL "],
        star_data[:stars]["MASS"],
    )

    for component in snapshotTypes(data_dict)

        for (block, values) in data_dict[component]

            if block ∈ ["POS ", "VEL "] && !isempty(values)
                data_dict[component][block] = rotation_matrix * values
            end

        end

    end

    return nothing

end

"""
    computeInertiaTensor(
        positions::Matrix{<:Unitful.Length},
        masses::Vector{<:Unitful.Mass},
    )::Matrix{Float64}

Compute the inertia tensor of a group of cells/particles.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Masses of the cells/particles.

# Returns

  - The inertia tensor.
"""
function computeInertiaTensor(
    positions::Matrix{<:Unitful.Length},
    masses::Vector{<:Unitful.Mass},
)::Matrix{Float64}

    # Check for missing data
    (
        !any(isempty, [positions, masses])  ||
        throw(ArgumentError("computeInertiaTensor: The inertia tensor is not defined for an \
        empty system"))
    )

    # Allocate memory
    J = zeros(Float64, (3, 3))

    # Compute the inertia tensor
    for (r, m) in zip(eachcol(positions), masses)

        x, y, z = r

        J[1, 1] += m * (y * y + z * z)
        J[2, 2] += m * (x * x + z * z)
        J[3, 3] += m * (x * x + y * y)

        J[1, 2] -= m * x * y
        J[1, 3] -= m * x * z
        J[2, 3] -= m * y * z

    end

    J[2, 1] = J[1, 2]
    J[3, 1] = J[1, 3]
    J[3, 2] = J[2, 3]

    return J

end

"""
    computeAMRotationMatrix(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass},
    )::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix that will turn the total angular momentum into the z axis; when view as an active (alibi) transformation.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.

# Returns

  - The rotation matrix.
"""
function computeAMRotationMatrix(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass},
)::Union{Matrix{Float64},UniformScaling{Bool}}

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return I

    # Compute the total angular momentum
    L = computeTotalAngularMomentum(positions, velocities, masses)

    # Rotation vector
    n = [L[2], -L[1], 0.0]

    # Angle of rotation
    θ = acos(L[3])

    return Matrix{Float64}(AngleAxis(θ, n...))

end

"""
    computePARotationMatrix(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass},
    )::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix that will turn the principal axis into the new coordinate system; when view as an passive (alias) transformation.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.

# Returns

  - The rotation matrix.
"""
function computePARotationMatrix(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass},
)::Union{Matrix{Float64},UniformScaling{Bool}}

    # Check for missing data
    !isempty(positions) || return I

    # Principal axis operator
    R = ustrip.(positions * positions')

    # Reverse the order of the eigenvectors, making the last column the eigenvector
    # with the largest eigenvalue, which should correspond to the new z axis
    pa = eigvecs(R)[:, end:-1:1]

    # Compute the total angular momentum
    L = computeTotalAngularMomentum(positions, velocities, masses; normal=true)

    # 3rd principal axis ≡ new z axis
    pa_z = pa[:, 3]

    # Rotate the principal axis as to align the third component with the angular momentum
    θ = acos(L ⋅ pa_z)
    n = cross(pa_z, L)
    aligned_pa = AngleAxis(θ, n...) * pa

    # The rotation matrix is made from the principal axis as rows
    rotation_matrix = aligned_pa'

    if det(rotation_matrix) < 0.0
        # If the determinant is < 0, that means that the chosen principal axis for the x and y
        # directions form a left-handed Cartesian reference system (x × y = -z). When applying
        # this as a rotation, the z axis will be flipped. So, in this case we swap the x and y
        # axis to get a right-handed Cartesian reference system (x × y = z) and generate the
        # correct rotation
        rotation_matrix[1, :], rotation_matrix[2, :] = rotation_matrix[2, :], rotation_matrix[1, :]
    end

    return Matrix{Float64}(rotation_matrix)

end

"""
    computeGlobalAMRotationMatrix(data_dict::Dict)::Union{Matrix{Float64},UniformScaling{Bool}}

Compute the rotation matrix that will turn the total angular momentum of the whole system, into the z axis; when view as an active (alibi) transformation.

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

  - The rotation matrix.
"""
function computeGlobalAMRotationMatrix(data_dict::Dict)::Union{Matrix{Float64},UniformScaling{Bool}}

    components = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), components)

    # Concatenate the positions, velocities, and masses of all the cells and particles in the system
    positions  = hcat([data_dict[component]["POS "] for component in components]...)
    velocities = hcat([data_dict[component]["VEL "] for component in components]...)
    masses     = vcat([data_dict[component]["MASS"] for component in components]...)

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return I

    (
        !logging[] ||
        @info("computeGlobalAMRotationMatrix: The rotation matrix will be computed using \
        $(components)")
    )

    return computeAMRotationMatrix(positions, velocities, masses)

end
