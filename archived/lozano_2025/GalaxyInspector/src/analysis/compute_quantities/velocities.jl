####################################################################################################
# Compute characteristic velocities and momentums
####################################################################################################

"""
    computeVcm(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Velocity}

Read the velocity of the center of mass of a given halo or subhalo.

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

  - The specified velocity.
"""
function computeVcm(data_dict::Dict, subfind_idx::NTuple{2,Int})::Vector{<:Unitful.Velocity}

    # If there are no subfind data, return the origin
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return zeros(typeof(1.0u"km*s^-1"), 3)
    end

    halo_idx, subhalo_rel_idx = subfind_idx

    # Load the necessary data
    n_subhalos_in_halo = data_dict[:group]["G_Nsubs"]
    g_vel = data_dict[:group]["G_Vel"]
    s_vel = data_dict[:subhalo]["S_Vel"]

    # Check that the requested halo index is within bounds
    n_halos = data_dict[:gc_data].header.n_groups_total

    (
        !iszero(n_halos) && !any(isempty, [n_subhalos_in_halo, g_vel, s_vel]) ||
        return zeros(typeof(1.0u"km*s^-1"), 3)
    )

    (
        0 < halo_idx <= n_halos ||
        throw(ArgumentError("computeVcm: There is only $(n_halos) FoF groups in \
        $(data_dict[:gc_data].path), so `halo_idx` = $(halo_idx) is out of bounds"))
    )

    # Select the halo velocity if `subhalo_rel_idx` == 0
    isPositive(subhalo_rel_idx) || return g_vel[:, halo_idx]

    # Check that the requested subhalo index is within bounds
    n_subfinds = n_subhalos_in_halo[halo_idx]

    if iszero(n_subfinds)

        (
            !logging[] ||
            @info("computeVcm: There are 0 subhalos in the FoF group $(halo_idx) from \
            $(data_dict[:gc_data].path), so the velocity will the halo velocity")
        )

        return g_vel[:, halo_idx]

    end

    (
        subhalo_rel_idx <= n_subfinds ||
        throw(ArgumentError("computeVcm: There is only $(n_subfinds) subhalos for the FoF \
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

    # Select the subhalo velocity
    return s_vel[:, subhalo_abs_idx]

end

"""
    computeVcm(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Velocity}

Read the velocity of the center of mass of a given subhalo.

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

  - The specified velocity.
"""
function computeVcm(data_dict::Dict, subhalo_abs_idx::Int)::Vector{<:Unitful.Velocity}

    # If there are no subfind data, return the origin
    if ismissing(data_dict[:gc_data].path) && !isSubfindActive(data_dict[:gc_data].path)
        return zeros(typeof(1.0u"km*s^-1"), 3)
    end

    s_vel = data_dict[:subhalo]["S_Vel"]

    # Check that the requested subhalo index is within bounds
    n_subgroups_total = data_dict[:gc_data].header.n_subgroups_total

    !iszero(n_subgroups_total) && !isempty(s_vel) || return zeros(typeof(1.0u"km*s^-1"), 3)

    (
        0 < subhalo_abs_idx <= n_subgroups_total ||
        throw(ArgumentError("computeCenter: There is only $(n_subgroups_total) subhalos in \
        $(data_dict[:gc_data].path), so `subhalo_abs_idx` = $(subhalo_abs_idx) is out of bounds"))
    )

    # Select the subhalo velocity
    return s_vel[:, subhalo_abs_idx]

end

"""
    computeVcm(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Velocity}

Compute the velocity of a characteristic center of the system.

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

  - The specified velocity.
"""
function computeVcm(data_dict::Dict, cm_type::Symbol)::Vector{<:Unitful.Velocity}

    if cm_type == :global_cm

        return computeGlobalVcm(data_dict)

    elseif cm_type ∈ keys(PARTICLE_INDEX)

        return computeComponentVcm(data_dict, cm_type)

    elseif cm_type == :zero

        return zeros(typeof(1.0u"km*s^-1"), 3)

    end

    throw(ArgumentError("computeVcm: `cm_type` can only be :global_cm, :zero or one of the keys of \
    `PARTICLE_INDEX` but I got :$(center_type)"))

end

"""
    computeComponentVcm(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Velocity}

Compute the velocity of the given component center of mass.

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
  - `component::Symbol`: Target component, it can be any of the keys of [`PARTICLE_INDEX`](@ref).

# Returns

  - The velocity of the center of mass.
"""
function computeComponentVcm(data_dict::Dict, component::Symbol)::Vector{<:Unitful.Velocity}

    # Load the necessary data
    velocities = data_dict[component]["VEL "]
    masses     = data_dict[component]["MASS"]

    # Check for missing data
    !any(isempty, [velocities, masses]) || return zeros(typeof(1.0u"km*s^-1"), 3)

    # Compute the total mass
    M = sum(masses)

    # Compute the velocity of the center of mass
    vcm = [sum(row .* masses) / M for row in eachrow(velocities)]

    return vcm

end

"""
    computeGlobalVcm(data_dict::Dict)::Vector{<:Unitful.Velocity}

Compute the velocity of the global center of mass.

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

  - The velocity of the global center of mass.
"""
function computeGlobalVcm(data_dict::Dict)::Vector{<:Unitful.Velocity}

    components = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["VEL "]), components)

    # Load the necessary data
    velocities = hcat([data_dict[component]["VEL "] for component in components]...)
    masses     = vcat([data_dict[component]["MASS"] for component in components]...)

    # Check for missing data
    !any(isempty, [velocities, masses]) || return zeros(typeof(1.0u"km*s^-1"), 3)

    # Compute the total mass
    M = sum(masses)

    # Compute the velocity of the center of mass
    vcm = [sum(row .* masses) / M for row in eachrow(velocities)]

    return vcm

end

@doc raw"""
    computeVcirc(
        data_dict::Dict;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

Compute the circular velocity of each particle of the given type, with respect to the origin.

The circular velocity of a particle is,

```math
v_\mathrm{circ} = \sqrt{\frac{\mathrm{G} \, M(r)}{r}} \, ,
```

where $r$ is the radial distance of the particle, and $M(r)$ is the total mass within a sphere of radius $r$.

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
  - `type::Symbol=:stars`: Target component.

# Returns

  - A tuple with two elements:

      + A vector with the radial distance of each particle to the origin.
      + A vector with the circular velocity of each particle.
"""
function computeVcirc(
    data_dict::Dict;
    type::Symbol=:stars,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

    # Compute the radial distance to each particle
    rs = computeDistance(data_dict[type]["POS "])

    # Check for missing data
    !isempty(rs) || return rs, Unitful.Velocity[]

    components = filter!(
        ts -> !isempty(data_dict[ts]["POS "]),
        [:stars, :gas, :halo, :black_hole],
    )

    # Concatenate the position and masses of all the cells and particles in the system
    distances = vcat([computeDistance(data_dict[component]["POS "]) for component in components]...)
    masses    = vcat([data_dict[component]["MASS"] for component in components]...)

    # Use the radial distances as bin edges for the mass histogram
    edges = [0.0u"kpc", rs...]

    # Compute to total mass within each particle radial distance
    M = similar(rs, eltype(masses))
    cumsum!(M, histogram1D(distances, masses, edges; empty_nan=false))

    # The mass histogram is a sorted array, so it is reverted to the unsorted order of `r`
    # to make `vcirc` the circular velocity of each particle in the order of the snapshot
    invpermute!(M, sortperm(rs))

    (
        !logging[] ||
        @info("computeVcirc: The circular velocity will be computed using $(components)")
    )

    vcirc = [iszero(r) ? 0.0u"km*s^-1" : sqrt(Unitful.G * m / r) for (m, r) in zip(M, rs)]

    return rs, vcirc

end

@doc raw"""
    computeVpolar(
        data_dict::Dict,
        component::Symbol;
        <keyword arguments>
    )::Vector{<:Unitful.Velocity}

Compute the cylindrical components of the velocity, $\mathbf{\vec{v}} = v_r \, \mathbf{e_r} + v_\theta \, \mathbf{e_\theta} + v_z \, \mathbf{e_z}$.

The speed in the radial direction expressed in Cartesian coordinates is

```math
v_r = \frac{x \, v_x + y \, v_y}{\sqrt(x^2 + y^2)} \, ,
```

in the tangential direction is

```math
v_\tau = \frac{x \, v_y - y \, v_x}{\sqrt(x^2 + y^2)} \, ,
```

and the speed in the z direction will be computes as

```math
v^*_z = v_z \, \mathrm{sign}(z) \, ,
```

in order to distinguish between inflows ($v^*_z < 0$) and outflows ($v^*_z > 0$).

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
  - `component::Symbol`: Which component will be calculated. The options are:

      + `:radial`     -> Stellar radial speed ($v_r$).
      + `:tangential` -> Stellar tangential speed ($v_\theta$).
      + `:zstar`      -> Stellar speed in the z direction, computed as $v_z \, \mathrm{sign}(z)$.
  - `type::Symbol=:stars`: Target cell/particle type.

# Returns

  - The chosen cylindricall component of the velocity.
"""
function computeVpolar(
    data_dict::Dict,
    component::Symbol;
    type::Symbol=:stars,
)::Vector{<:Unitful.Velocity}

    # Load the necessary data
    positions = data_dict[type]["POS "]
    velocities = data_dict[type]["VEL "]

    x = positions[1, :]
    y = positions[2, :]
    z = positions[3, :]

    vx = velocities[1, :]
    vy = velocities[2, :]
    vz = velocities[3, :]

    if component == :radial

        # Compute the radial component
        vp = @. (x * vx + y * vy) / sqrt(x^2 + y^2)

    elseif component == :tangential

        # Compute the tangential component
        vp = @. (x * vy - y * vx) / sqrt(x^2 + y^2)

    elseif component == :zstar

        # Compute the z component
        vp = @. vz * sign(z)

    else

        throw(ArgumentError("computeVpolar: `component` can only be :radial, :tangential \
        or :zstar, but I got :$(component)"))

    end

    return vp

end

"""
    computeAngularMomentum(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass},
    )::Vector{Vector{<:AngularMomentum}}

Compute the angular momentum of each cell/particle, with respect to the origin.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.

# Returns

  - The angular momentum of each cell/particle.
"""
function computeAngularMomentum(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass},
)::Vector{Vector{<:AngularMomentum}}

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return [0.0, 0.0, 1.0]

    iterator = zip(masses, eachcol(positions), eachcol(velocities))

    return map(x -> x[1] .* cross(x[2], x[3]), iterator)

end

"""
    computeTotalAngularMomentum(
        positions::Matrix{<:Unitful.Length},
        velocities::Matrix{<:Unitful.Velocity},
        masses::Vector{<:Unitful.Mass};
        <keyword arguments>
    )::Vector{<:Number}

Compute the total angular momentum of a group of cells/particles, with respect to the origin.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.
  - `normal::Bool=true`: If the result will be normalized.

# Returns

  - The angular momentum.
"""
function computeTotalAngularMomentum(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass};
    normal::Bool=true,
)::Vector{<:Number}

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return [0.0, 0.0, 1.0]

    iterator = zip(masses, eachcol(positions), eachcol(velocities))

    # Compute the total angular momentum
    L = mapreduce(x -> x[1] .* cross(x[2], x[3]), +, iterator)

    return normal ? normalize!(ustrip.(L)) : L

end

"""
    computeGlobalAngularMomentum(data_dict::Dict; <keyword arguments>)::Vector{<:Number}

Compute the total angular momentum with respect to the origin of the whole system.

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
  - `normal::Bool=true`: If the result will be normalized.

# Returns

  - The angular momentum.
"""
function computeGlobalAngularMomentum(data_dict::Dict; normal::Bool=true)::Vector{<:Number}

    components = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), components)

    # Concatenate the position, velocities, and masses of all the cells and particles in the system
    positions  = hcat([data_dict[component]["POS "] for component in components]...)
    velocities = hcat([data_dict[component]["VEL "] for component in components]...)
    masses     = vcat([data_dict[component]["MASS"] for component in components]...)

    # Check for missing data
    !any(isempty, [positions, velocities, masses]) || return [0.0, 0.0, 1.0]

    (
        !logging[] ||
        @info("computeGlobalAngularMomentum: The angular momentum will be computed using \
        $(components)")
    )

    return computeTotalAngularMomentum(positions, velocities, masses; normal)

end

@doc raw"""
    computeSpinParameter(
        positions::Matrix{<:Unitful.Length},
        masses::Vector{<:Unitful.Mass},
        velocities::Matrix{<:Unitful.Velocity};
        <keyword arguments>
    )::Float64

Compute the spin parameter for a system of cells/particles, with respect to the origin.

The spin parameter was originally defined by Peebles (1969) as,

```math
\lambda = \frac{J \, \sqrt{E}}{G \, M^{5/2}} \, ,
```

where $J$ is the norm of the total angular momentum, $M$ the total mass, $G$ the gravitational constant, and

```math
E = |E_P + E_k| \, ,
```

where $E_P$ is the total potencial energy and $E_k$ is the total kinetic energy (including thermal energy of the gas).

Due to the computational complexity of calculating $E_P$ for a large group of particles, Bullock et al. (2001) proposed an alternative definition of the spin parameter,

```math
\lambda = \frac{J}{\sqrt{2} \, M \, R \, V} \, ,
```

where $J$ is the norm of the total angular momentum inside a sphere of radius $R$ containing mass $M$, and

```math
V = \sqrt{\frac{G \, M}{R}} \, ,
```

is the circular velocity.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `velocities::Matrix{<:Unitful.Velocity}`: Velocities of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `masses::Vector{<:Unitful.Mass}`: Mass of every cell/particle.
  - `R::Unitful.Length=DISK_R`: Radius.

# Returns

  - The spin parameter.

# References

P. J. E. Peebles (1969). *Origin of the Angular Momentum of Galaxies*. Astrophysical Journal, **155**, 393. [doi:10.1086/149876](https://doi.org/10.1086/149876)

J. S. Bullock et al. (2001). *A Universal Angular Momentum Profile for Galactic Halos*. The Astrophysical Journal, **555(1)**, 240. [doi:10.1086/321477](https://doi.org/10.1086/321477)

J. Zjupa et al. (2017). *Angular momentum properties of haloes and their baryon content in the Illustris simulation*. Monthly Notices of the Royal Astronomical Society, **466(2)**, 1625–1647. [doi:10.1093/mnras/stw2945](https://doi.org/10.1093/mnras/stw2945)
"""
function computeSpinParameter(
    positions::Matrix{<:Unitful.Length},
    velocities::Matrix{<:Unitful.Velocity},
    masses::Vector{<:Unitful.Mass};
    R::Unitful.Length=DISK_R,
)::Float64

    (
        !any(isempty, [positions, velocities, masses]) ||
        throw(ArgumentError("computeSpinParameter: The spin parameter of an empty dataset \
        is undefined"))
    )

    # Find cells/particles within the virial radius
    idx = map(x -> x <= R, computeDistance(positions))

    # Compute the total mass within the virial radius
    M = sum(masses[idx]; init=0.0u"Msun")

    # Compute the norm of the total angular momentum
    J = norm(
        computeTotalAngularMomentum(
            positions[:, idx],
            velocities[:, idx],
            masses[idx];
            normal=false,
        )
    )

    return uconvert(Unitful.NoUnits, J / sqrt(2.0 * R * Unitful.G * M^3))

end

"""
    computeGlobalSpinParameter(data_dict::Dict; <keyword arguments>)::Float64

Compute the spin parameter of the whole system.

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
  - `R::Unitful.Length=DISK_R`: Radius.

# Returns

  - The spin parameter.
"""
function computeGlobalSpinParameter(data_dict::Dict; R::Unitful.Length=DISK_R)::Float64

    components = snapshotTypes(data_dict)

    filter!(ts -> !isempty(data_dict[ts]["POS "]), components)

    # Concatenate the position and masses of all the cells and particles in the system
    positions  = hcat([data_dict[component]["POS "] for component in components]...)
    velocities = hcat([data_dict[component]["VEL "] for component in components]...)
    masses     = vcat([data_dict[component]["MASS"] for component in components]...)

    (
        !logging[] ||
        @info("computeGlobalSpinParameter: The spin parameter will be computed using $(components)")
    )

    # Compute the total spin parameter
    return computeSpinParameter(positions, velocities, masses; R)

end

@doc raw"""
    computeCircularity(data_dict::Dict; <keyword arguments>)::Vector{Float64}

Compute the circularity of each particle of the given type, with respect to the origin and the $z$ direction [0, 0, 1].

The circularity of a particle is,

```math
\epsilon = j_z / j_\mathrm{circ} \, ,
```

where $j_z$ is the $z$ component of its specific angular momentum, and $j_\mathrm{circ}$ is the specific angular momentum of a circular orbit,

```math
j_\mathrm{circ} = r \, v_\mathrm{circ} = \sqrt{\mathrm{G} \, r \, M(r)} \, ,
```

where $r$ is the radial distance of the particle, and $M(r)$ is the total mass within a sphere of radius $r$.

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
  - `type::Symbol=:stars`: Target component.

# Returns

  - The circularity $\epsilon$ of each particle.
"""
function computeCircularity(data_dict::Dict; type::Symbol=:stars)::Vector{Float64}

    # Load the necessary data
    positions  = data_dict[type]["POS "]
    velocities = data_dict[type]["VEL "]

    # Check for missing data
    !any(isempty, [positions, velocities]) || return Float64[]

    # Compute the specific angular momentum in the z direction
    jzs = [x[1] * v[2] - x[2] * v[1] for (x, v) in zip(eachcol(positions), eachcol(velocities))]

    # Compute the circular velocities and the radial distances
    rs, vcircs = computeVcirc(data_dict; type)

    stellar_circularity = [
        any(iszero, [r, vcirc]) ? 0.0 : ustrip(Unitful.NoUnits, jz / (r * vcirc)) for
        (jz, r, vcirc) in zip(jzs, rs, vcircs)
    ]

    return stellar_circularity

end
