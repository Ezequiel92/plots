####################################################################################################
# Data analysis functions
####################################################################################################

####################################################################################################
# Signature for the plotSnapshot function in ./src/plotting/pipelines.jl
####################################################################################################
#
# A data analysis functions for plotSnapshot must take a dictionary with the following shape:
#
#   + :sim_data          -> ::Simulation (see the Simulation struct in ./src/constants/globals.jl).
#   + :snap_data         -> ::Snapshot (see the Snapshot struct in ./src/constants/globals.jl).
#   + :gc_data           -> ::GroupCatalog (see the GroupCatalog struct in ./src/constants/globals.jl).
#   + cell/particle type -> (block -> data of block, block -> data of block, ...).
#   + cell/particle type -> (block -> data of block, block -> data of block, ...).
#   + cell/particle type -> (block -> data of block, block -> data of block, ...).
#   + ...
#   + groupcat type      -> (block -> data of block, block -> data of block, ...).
#   + groupcat type      -> (block -> data of block, block -> data of block, ...).
#   + groupcat type      -> (block -> data of block, block -> data of block, ...).
#   + ...
#
# and return one or more vectors or matrices. It should return `nothing` if the input data has
# some problem that prevents computation (e.g. is empty).
#
# Expected signature:
#
#   da_function(data_dict, args...; kwargs...) -> (processed_data, ...)  or `nothing`
#
# where:
#
#   - data_dict::Dict
#   - processed_data::VecOrMat{<:Number}
#
####################################################################################################

"""
    daRotationCurve(
        data_dict::Dict,
        R::Unitful.Length;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

Compute a rotation curve.

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
  - `R::Unitful.Length`: Maximum radius.
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

  - A tuple with two elements:

      + A vector with the distances to each star.
      + A vector with the circular velocity of each star.
"""
function daRotationCurve(
    data_dict::Dict,
    R::Unitful.Length;
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Velocity}}

    filtered_dd = filterData(data_dict; filter_function)

    # Compute the circular velocities and the radial distances of each star
    r, vcirc = computeVcirc(filtered_dd)

    # Only leave the data within a sphere of radius `R`
    rangeCut!(r, vcirc, (0.0u"kpc", R))

    # Sort the arrays radialy
    idx = sortperm(r)

    return r[idx], vcirc[idx]

end

"""
    daKennicuttSchmidtLaw(
        data_dict::Dict,
        grid::CubicGrid,
        quantity::Symbol;
        <keyword arguments>
    )::Union{NTuple{2,Vector{<:Float64}},Nothing}

Compute the gas mass surface density and the SFR surface density, used in the Kennicutt-Schmidt law.

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
  - `grid::CubicGrid`: Cubic grid.
  - `quantity::Symbol=:molecular_mass`: Quantity for the x axis. The options are:

      + `:gas_mass`          -> Gas mass surface density. This one will be plotted with the results of Kennicutt (1998).
      + `:molecular_mass`    -> Molecular mass surface density. This one will be plotted with the results of Bigiel et al. (2008).
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006). This one will be plotted with the results of Bigiel et al. (2008).
      + `:neutral_mass`      -> Neutral mass surface density. This one will be plotted with the results of Bigiel et al. (2008).
  - `type::Symbol=:cells`: If the gas surface density will be calculated assuming the gas is in `:particles` or in Voronoi `:cells`.
  - `reduce::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection, averaging the value of neighboring pixels. It has to divide the size of `grid` exactly.
  - `stellar_ff::Function=filterNothing`: Filter function for the stars. It has to be a function with the signature:

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
  - `gas_ff::Function=filterNothing`: Filter function for the gas. It has to be a function with the signature:

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

  - A tuple with two elements:

      + A vector with log10(ΣH / M⊙ * kpc^-2).
      + A vector with log10(Σsfr / M⊙ * yr^-1 * kpc^-2).

    It returns `nothing` if any of the necessary quantities are missing.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function daKennicuttSchmidtLaw(
    data_dict::Dict,
    grid::CubicGrid,
    quantity::Symbol;
    type::Symbol=:cells,
    reduce_factor::Int=1,
    stellar_ff::Function=filterNothing,
    gas_ff::Function=filterNothing,
)::Union{NTuple{2,Vector{<:Float64}},Nothing}

    (
        quantity ∈ [:gas_mass, :molecular_mass, :br_molecular_mass, :neutral_mass] ||
        throw(ArgumentError("daKennicuttSchmidtLaw: `quantity` can only be :gas_mass, \
        :molecular_mass, :br_molecular_mass or :neutral_mass, but I got :$(quantity)"))
    )

    # Log factor to go from stellar surface density to SFR surface density
    # log10(Σsfr) = log10(Σ*) - log10Δt
    log10Δt = log10(ustrip(u"yr", AGE_RESOLUTION))

    _, _, stellar_density = daDensity2DProjection(
        data_dict,
        grid,
        :stellar_mass,
        :particles;
        reduce_factor,
        filter_function=stellar_ff,
    )

    _, _, gas_density = daDensity2DProjection(
        data_dict,
        grid,
        quantity,
        type;
        reduce_factor,
        filter_function=gas_ff,
    )

    x_axis = vec(gas_density)
    y_axis = vec(stellar_density)

    # Delete 0s and NaNs in the data vectors
    x_idxs = map(x -> isnan(x) || iszero(x), x_axis)
    y_idxs = map(x -> isnan(x) || iszero(x), y_axis)

    delete_idxs = x_idxs ∪ y_idxs

    deleteat!(x_axis, delete_idxs)
    deleteat!(y_axis, delete_idxs)

    !any(isempty, [x_axis, y_axis]) || return nothing

    return x_axis, y_axis .- log10Δt

end

"""
    daMolla2015(
        data_dict::Dict,
        grid::CircularGrid,
        quantity::Symbol;
        <keyword arguments>
    )::Union{
        Tuple{
            Vector{<:Unitful.Length},
            <:Union{Vector{<:SurfaceDensity},Vector{<:MassFlowDensity},Vector{Float64}}
        },
        Nothing,
    }

Compute a profile for the Milky Way, compatible with the experimental data in Mollá et al. (2015).

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
  - `grid::CircularGrid`: Circular grid.
  - `quantity::Symbol`: Quantity. The options are:

      + `:stellar_area_density`      -> Stellar area mass density.
      + `:molecular_area_density`    -> Molecular mass surface density.
      + `:br_molecular_area_density` -> Molecular mass surface density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_area_density`       -> Atomic hydrogen area mass density.
      + `:sfr_area_density`          -> Star formation rate area density, for the last `AGE_RESOLUTION`.
      + `:O_stellar_abundance`       -> Stellar abundance of oxygen, as ``12 + \\log_{10}(\\mathrm{O \\, / \\, H})``.
      + `:N_stellar_abundance`       -> Stellar abundance of nitrogen, as ``12 + \\log_{10}(\\mathrm{N \\, / \\, H})``.
      + `:C_stellar_abundance`       -> Stellar abundance of carbon, as ``12 + \\log_{10}(\\mathrm{C \\, / \\, H})``.
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

  - A tuple with two elements:

      + A vector with the position of each ring.
      + A vector with the `quantity` area density of each ring.

    It returns `nothing` if any of the necessary quantities are missing.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function*. Monthly Notices of the Royal Astronomical Society **451(4)**, 3693–3708. [doi:10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)
"""
function daMolla2015(
    data_dict::Dict,
    grid::CircularGrid,
    quantity::Symbol;
    filter_function::Function=filterNothing,
)::Union{
    Tuple{
        Vector{<:Unitful.Length},
        <:Union{Vector{<:SurfaceDensity},Vector{<:MassFlowDensity},Vector{Float64}},
    },
    Nothing,
}

    filtered_dd = filterData(data_dict; filter_function)

    if quantity == :stellar_area_density

        positions   = filtered_dd[:stars]["POS "]
        masses      = filtered_dd[:stars]["MASS"]
        norm_values = Number[]
        f           = identity
        density     = true

    elseif quantity == :sfr_area_density

        positions   = filtered_dd[:stars]["POS "]
        masses      = computeSFR(filtered_dd; age_resol=AGE_RESOLUTION)
        norm_values = Number[]
        f           = identity
        density     = true

    elseif quantity == :molecular_area_density

        positions   = filtered_dd[:gas]["POS "]
        masses      = computeMass(filtered_dd, :molecular)
        norm_values = Number[]
        f           = identity
        density     = true

    elseif quantity == :br_molecular_area_density

        positions   = filtered_dd[:gas]["POS "]
        masses      = computeMass(filtered_dd, :br_molecular)
        norm_values = Number[]
        f           = identity
        density     = true

    elseif quantity == :atomic_area_density

        positions   = filtered_dd[:gas]["POS "]
        masses      = computeMass(filtered_dd, :atomic)
        norm_values = Number[]
        f           = identity
        density     = true

    elseif quantity == :O_stellar_abundance

        positions   = filtered_dd[:stars]["POS "]
        masses      = computeElementMass(filtered_dd, :stars, :O) ./ ATOMIC_WEIGHTS[:O]
        norm_values = computeElementMass(filtered_dd, :stars, :H) ./ ATOMIC_WEIGHTS[:H]
        f           = x -> 12 .+ log10.(x)
        density     = false

    elseif quantity == :N_stellar_abundance

        positions   = filtered_dd[:stars]["POS "]
        masses      = computeElementMass(filtered_dd, :stars, :N) ./ ATOMIC_WEIGHTS[:N]
        norm_values = computeElementMass(filtered_dd, :stars, :H) ./ ATOMIC_WEIGHTS[:H]
        f           = x -> 12 .+ log10.(x)
        density     = false

    elseif quantity == :C_stellar_abundance

        positions   = filtered_dd[:stars]["POS "]
        masses      = computeElementMass(filtered_dd, :stars, :N) ./ ATOMIC_WEIGHTS[:N]
        norm_values = computeElementMass(filtered_dd, :stars, :H) ./ ATOMIC_WEIGHTS[:H]
        f           = x -> 12 .+ log10.(x)
        density     = false

    else

        throw(ArgumentError("daMolla2015: I don't recognize the quantity :$(quantity)"))

    end

    # Return `nothing` if any of the necessary quantities are missing
    !any(isempty, [positions, masses]) || return nothing

    density_profile = f(computeParticleProfile(positions, masses, grid; norm_values, total=true, density))

    return grid.grid, density_profile

end

"""
    daProfile(
        data_dict::Dict,
        quantity::Symbol,
        grid::CircularGrid;
        <keyword arguments>
    )::Union{Tuple{Vector{<:Unitful.Length},Vector{<:Number}},Nothing}

Compute a profile.

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
  - `quantity::Symbol`: Target quantity. The options are the same as for [`scatterQty`](@ref):

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
  - `grid::CircularGrid`: Circular grid.
  - `flat::Bool=true`: If the profile will be 2D, using rings, or 3D, using spherical shells.
  - `total::Bool=true`: If the sum (default) or the mean of `quantity` will be computed for each bin.
  - `cumulative::Bool=false`: If the profile will be accumulated or not.
  - `density::Bool=false`: If the profile will be of the density of `quantity`.
  - `fractions::Bool=false`: If a profile of the gas mass fractions will be calculated. It is only valid with `quantity` equal to :neutral_mass, :molecular_mass, :br_molecular_mass, :atomic_mass or :ionized_mass, and it forces `total` = true, `cumulative` = false, and `density` = false.
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

  - A tuple with two elements:

      + A vector with the position of each ring or spherical shells.
      + A vector with the value `quantity` in each each ring or spherical shells.

    It returns `nothing` if any of the necessary quantities are missing.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function daProfile(
    data_dict::Dict,
    quantity::Symbol,
    grid::CircularGrid;
    flat::Bool=true,
    total::Bool=true,
    cumulative::Bool=false,
    density::Bool=false,
    fractions::Bool=false,
    filter_function::Function=filterNothing,
)::Union{Tuple{Vector{<:Unitful.Length},Vector{<:Number}},Nothing}

    filtered_dd = filterData(data_dict; filter_function)

    # Get the cell/particle type
    type = first(keys(plotParams(quantity).request))

    # Read the positions and values
    positions = filtered_dd[type]["POS "]
    values    = scatterQty(filtered_dd, quantity)

    n_pos = size(positions, 2)
    n_val = length(values)

    # Check consistency in the number of positions and values
    (
        n_pos == n_val ||
        throw(ArgumentError("daProfile: `positions` and `values` should have the same number of \
        elements, but length(`positions`) = $(n_pos) != length(`values`) = $(n_val). \
        Check that the same cell/particle type was selected when using `plotParams` \
        and `scatterQty`."))
    )

    # Return `nothing` if any of the necessary quantities are missing
    !any(iszero, [n_pos, n_val]) || return nothing

    if fractions

        (
            quantity ∈ [
                :neutral_mass,
                :molecular_mass,
                :br_molecular_mass,
                :atomic_mass,
                :ionized_mass,
            ] ||
            throw(ArgumentError("daProfile: If `fractions`= true, quantity must be \
            :neutral_mass, :molecular_mass, :br_molecular_mass, :atomic_mass or :ionized_mass, but \
            I got `quantity` = :$(quantity)"))
        )

        total       = true
        cumulative  = false
        density     = false
        norm_values = scatterQty(filtered_dd, :gas_mass)

    else

        norm_values = Number[]

    end

    profile = computeParticleProfile(
        positions,
        values,
        grid;
        norm_values,
        flat,
        total,
        cumulative,
        density,
    )

    return grid.grid, profile

end

"""
    daBandProfile(
        data_dict::Dict,
        quantity::Symbol,
        grid::CircularGrid;
        <keyword arguments>
    )::Union{
        Tuple{Vector{<:Unitful.Length},Vector{<:Number},Vector{<:Number},Vector{<:Number}},
        Tuple{Vector{<:Unitful.Length},Vector{<:Number},Vector{<:Number}},
        Nothing,
    }

Compute the profile of a mean quantity with error bars or bands.

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
  - `quantity::Symbol`: Target quantity. The possibilities are:

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
  - `grid::CircularGrid`: Circular grid.
  - `flat::Bool=true`: If the profile will be 2D, using rings, or 3D, using spherical shells.
  - `error_bar::Bool=false`: If the returned values will be compatible with `errorbars!` or with `band!` (default).
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

  - A tuple with two elements:

      + A vector with the position of each ring or spherical shells.
      + A vector with the value `quantity` in each each ring or spherical shells.

    It returns `nothing` if any of the necessary quantities are missing.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function daBandProfile(
    data_dict::Dict,
    quantity::Symbol,
    grid::CircularGrid;
    flat::Bool=true,
    error_bar::Bool=false,
    filter_function::Function=filterNothing,
)::Union{
    Tuple{Vector{<:Unitful.Length},Vector{<:Number},Vector{<:Number},Vector{<:Number}},
    Tuple{Vector{<:Unitful.Length},Vector{<:Number},Vector{<:Number}},
    Nothing,
}

    filtered_dd = filterData(data_dict; filter_function)

    # Get the cell/particle type
    type = first(keys(plotParams(quantity).request))

    # Read the positions and values
    positions = filtered_dd[type]["POS "]
    values    = scatterQty(filtered_dd, quantity)

    n_pos = size(positions, 2)
    n_val = length(values)

    # Check consistency in the number of positions and values
    (
        n_pos == n_val ||
        throw(ArgumentError("daBandProfile: `positions` and `values` should have the same number \
        of elements, but length(`positions`) = $(n_pos) != length(`values`) = $(n_val). \
        Check that the same cell/particle type was selected when using `plotParams` \
        and `scatterQty`"))
    )

    # Return `nothing` if any of the necessary quantities are missing
    !any(iszero, [n_pos, n_val]) || return nothing

    mean, std = computeParticleBandProfile(positions, values, grid; flat)

    !error_bar || return grid.grid, mean, std, std

    return grid.grid, mean .- std, mean .+ std

end

"""
    daStellarHistory(
        data_dict::Dict;
        <keyword arguments>
    )::Union{Tuple{Vector{<:Unitful.Time},Vector{<:Number}},Nothing}

Compute the evolution of a given stellar `quantity` using the stellar ages at a given instant in time.

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
  - `quantity::Symbol=:sfr`: Target quantity. The options are:

      + `:sfr`                 -> The star formation rate.
      + `:ssfr`                -> The specific star formation rate.
      + `:stellar_mass`        -> Stellar mass.
      + `:stellar_metallicity` -> Mass fraction of all elements above He in the stars (solar units).
  - `n_bins::Int=100`: Number of bins (time intervals).
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

  - A tuple with two elements:

      + A vector with the physical times.
      + A vector with the values of `quantity` at each time.
"""
function daStellarHistory(
    data_dict::Dict;
    quantity::Symbol=:sfr,
    n_bins::Int=100,
    filter_function::Function=filterNothing,
)::Union{Tuple{Vector{<:Unitful.Time},Vector{<:Number}},Nothing}

    filtered_dd = filterData(data_dict; filter_function)

    birth_ticks = filtered_dd[:stars]["GAGE"]
    masses      = filtered_dd[:stars]["MASS"]

    # Return `nothing` if there are less than 3 stars
    !any(x -> x <= 2, length.([birth_ticks, masses])) || return nothing

    # Compute the stellar birth dates
    if filtered_dd[:sim_data].cosmological
        # Go from scale factor to physical time
        birth_times = computeTime(birth_ticks, filtered_dd[:snap_data].header)
    else
        birth_times = birth_ticks
    end

    # Compute the birth time range
    min, max = extrema(birth_times)

    # Compute the total stellar mass in each time bin
    grid = CircularGrid(max, n_bins; shift=min)

    stellar_masses = histogram1D(birth_times, masses, grid; empty_nan=false)

    # Compute the time axis
    bin_width = (max - min) / n_bins
    x_axis = collect(range(min + (bin_width * 0.5), length=n_bins, step=bin_width))

    # Compute the stellar quantity
    if quantity == :sfr

        y_axis = stellar_masses ./ bin_width

    elseif quantity == :ssfr

        accu_mass = cumsum(stellar_masses)
        y_axis = (stellar_masses ./ bin_width) ./ accu_mass

    elseif quantity == :stellar_mass

        y_axis = cumsum(stellar_masses)

    elseif quantity == :stellar_metallicity

        metallicities = computeMetalMass(data_dict, :stars)

        # Return `nothing` if any of the necessary quantities are missing
        !isempty(metallicities) || return nothing

        stellar_metallicities = histogram1D(birth_times, metallicities, grid; empty_nan=false)

        Z = cumsum(stellar_metallicities) ./ cumsum(stellar_masses)
        y_axis = Z ./ SOLAR_METALLICITY

    else

        throw(ArgumentError("daStellarHistory: `quantity` can only be :sfr, :ssfr, \
        :stellar_metallicity or :stellar_mass, but I got :$(quantity)"))

    end

    return x_axis, y_axis

end

"""
    daLineHistogram(
        data_dict::Dict,
        quantity::Symbol,
        grid::LinearGrid;
        <keyword arguments>
    )::Union{Tuple{Vector{<:Number},Vector{<:Number}},Nothing}

Compute a 1D histogram of a given `quantity`, normalized to the maximum number of counts.

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
  - `grid::LinearGrid`: Linear grid.
  - `type::Symbol`: Type of cell/particle.
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
  - `norm::Int=0`: Number of count that will be use to normalize the histogram. If left as 0, the histogram will be normalize with the maximum bin count.

# Returns

  - A tuple with two elements:

      + A vector with the value corresponding to each bin.
      + A vector with the counts, normalized to the maximum value.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function daLineHistogram(
    data_dict::Dict,
    quantity::Symbol,
    grid::LinearGrid,
    type::Symbol;
    filter_function::Function=filterNothing,
    norm::Int=0,
)::Union{Tuple{Vector{<:Number},Vector{<:Number}},Nothing}

    # Compute the values
    scatter_qty = scatterQty(data_dict, quantity)

    if isempty(scatter_qty)

        (
            !logging[] ||
            @warn("daLineHistogram: There is no data for the quantity :$(quantity), \\
            even before filtering")

        )

        return nothing

    end

    # Filter after computing the values, to preserve quantities that depends on
    # global properties (e.g. global gravitational potential)
    idxs = filter_function(data_dict)[type]

    values = scatter_qty[idxs]

    if isempty(values)

        (
            !logging[] ||
            @warn("daLineHistogram: After filtering, there is no data left for the quantity \
            :$(quantity)")

        )

        return nothing

    end

    # Compute the quantity histogram
    counts = histogram1D(values, grid)

    if logging[]

        clean_values = filter(x -> !isnan(x) && !isinf(x), values)

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

        @info(
            "\nHistogram statistics \
            \n  Simulation: $(basename(data_dict[:sim_data].path)) \
            \n  Snapshot:   $(data_dict[:snap_data].global_index) \
            \n  Quantity:   $(quantity) \
            \n  Type:       $(type) \
            \n  Max bin:    $(grid.grid[argmax(counts)]) \
            \n  Max count:  $(maximum(counts)) \
            \n  Min - Max:  $(min_max_v) \
            \n  Mean:       $(mean_v) \
            \n  Median:     $(meadian_v) \
            \n  Mode:       $(mode_v)"
        )

    end

    # Normalize the counts
    norm_counts = isPositive(norm) ? counts ./ norm : counts ./ maximum(counts)

    return grid.grid, norm_counts

end

"""
    daDensity2DProjection(
        data_dict::Dict,
        grid::CubicGrid,
        quantity::Symbol,
        field_type::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Union{Matrix{Float64},Vector{Float64}}}

Project a 3D density field into a given plane.

!!! note

    If the source of the field are particles, a simple 2D histogram is used. If they are Voronoi cells instead, the density of the cells that cross the line of sight of each pixel are added up.

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
  - `grid::CubicGrid`: Cubic grid.
  - `quantity::Symbol`: Which density will be calculated. The options are:

      + `:stellar_mass`      -> Stellar density.
      + `:gas_mass`          -> Gas density.
      + `:hydrogen_mass`     -> Hydrogen density.
      + `:dm_mass`           -> Dark matter density.
      + `:bh_mass`           -> Black hole density.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) density.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) density.
      + `:ionized_mass`      -> Ionized hydrogen (``\\mathrm{HII}``) density.
      + `:neutral_mass`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) density.
      + `:stellar_gas_mass`  -> Stellar gas mass (according to out SF model).
  - `field_type::Symbol`: If the source of the field are `:particles` or Voronoi `:cells`.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection. If `reduce_grid` = :square, the new values will be computed averaging the values of neighboring pixels. `reduce_factor` has to divide the size of `grid` exactly. If `reduce_grid` = :circular, the new values will be computed averaging the values of the pixels the fall within each of the `reduce_factor` concentric rings.
  - `reduce_grid::Symbol=:square`: Type of grid to reduce the resolution of the result. The options are:

      + `:square`    -> The density distribution will be reduced into a regular square grid, with a resolution `reduce_factor` times lower than `grid`. This emulates the way the surface densities are measured in observations. `reduce_factor` = 1 means no reduction in resolution.
      + `:circular` -> The density distribution will be reduced into a flat circular grid, formed by a series of `reduce_factor` concentric rings. This emulates the traditonal way the Kennicutt-Schmidt law is measured in simulations. `reduce_factor` = 1 means that the result will be a single point, the opposite of the `reduce_grid` = :square case.
  - `projection_plane::Symbol=:xy`: Projection plane. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
  - `m_unit::Unitful.Units=u"Msun"`: Mass unit.
  - `l_unit::Unitful.Units=u"kpc"`: Length unit.
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

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the ``\\log_{10}`` of the density at each point of the 2D grid.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function daDensity2DProjection(
    data_dict::Dict,
    grid::CubicGrid,
    quantity::Symbol,
    field_type::Symbol;
    reduce_factor::Int=1,
    reduce_grid::Symbol=:square,
    projection_plane::Symbol=:xy,
    m_unit::Unitful.Units=u"Msun",
    l_unit::Unitful.Units=u"kpc",
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Union{Matrix{Float64},Vector{Float64}}}

    filtered_dd = filterData(data_dict; filter_function)

    # Set the cell/particle type
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
        throw(ArgumentError("daDensity2DProjection: I don't recognize the quantity :$(quantity)"))
    end

    # For comological simulations with comoving units, correct
    # the density so it is always in physical units
    if !PHYSICAL_UNITS && filtered_dd[:sim_data].cosmological
        # Correction factor for the area
        # A [physical units] = A [comoving units] * a0^2
        physical_factor = filtered_dd[:snap_data].scale_factor^2
    else
        physical_factor = 1.0
    end

    # Load the cell/particle positions
    positions = filtered_dd[component]["POS "]

    # Compute the masses of the target quantity
    masses = scatterQty(filtered_dd, quantity)

    # If any of the necessary quantities are missing return an empty density field
    if any(isempty, [masses, positions])
        return grid.x_ticks, grid.y_ticks, fill(NaN, (grid.n_bins, grid.n_bins))
    end

    if field_type == :cells

        # Compute the volume of each cell
        cell_volumes = filtered_dd[component]["MASS"] ./ filtered_dd[component]["RHO "]

        # Compute the densities of the target quantity
        densities = ustrip.(m_unit * l_unit^-3, masses ./ cell_volumes)

        # Load the volume and area of the voxels
        voxel_volume = ustrip(l_unit^3, grid.bin_volume)
        voxel_area   = ustrip(l_unit^2, grid.bin_area)

        # Allocate memory
        physical_grid = Matrix{Float64}(undef, 3, grid.n_bins^3)

        # Compute the tree for a nearest neighbor search
        kdtree = KDTree(ustrip.(l_unit, positions))

        # Reshape the grid to conform to the way `nn` expect the matrix to be structured
        for i in eachindex(grid.grid)
            physical_grid[1, i] = ustrip(l_unit, grid.grid[i][1])
            physical_grid[2, i] = ustrip(l_unit, grid.grid[i][2])
            physical_grid[3, i] = ustrip(l_unit, grid.grid[i][3])
        end

        # Find the nearest cell to each voxel
        idxs, _ = nn(kdtree, physical_grid)

        # Allocate memory
        mass_grid = similar(grid.grid, Float64)

        # Compute the mass in each voxel
        for i in eachindex(grid.grid)
            mass_grid[i] = densities[idxs[i]] * voxel_volume
        end

        # Project `mass_grid` to the target plane
        if projection_plane == :xy
            dims = 3
        elseif projection_plane == :xz
            # Project across dimension 1 to keep it consistent with :xz for `field_type` = :particles
            dims = 1
        elseif projection_plane == :yz
            # Project across dimension 2 to keep it consistent with :yz for `field_type` = :particles
            dims = 2
        else
            throw(ArgumentError("daDensity2DProjection: The argument `projection_plane` must be \
            :xy, :xz or :yz, but I got :$(projection_plane)"))
        end

        density = dropdims(sum(mass_grid; dims) ./ voxel_area; dims)

    elseif field_type == :particles

        # Project the particles to the given plane
        if projection_plane == :xy
            pos_2D = positions[[1, 2], :]
        elseif projection_plane == :xz
            pos_2D = positions[[1, 3], :]
        elseif projection_plane == :yz
            pos_2D = positions[[2, 3], :]
        else
            throw(ArgumentError("daDensity2DProjection: The argument `projection_plane` must be \
            :xy, :xz or :yz, but I got :$(projection_plane)"))
        end

        # Compute the 2D histogram
        density = ustrip.(
            m_unit * l_unit^-2,
            histogram2D(pos_2D, masses, flattenGrid(grid); empty_nan=false) ./ grid.bin_area,
        )

    else

        throw(ArgumentError("daDensity2DProjection: The argument `field_type` must be :cells or \
        :particles, but I got :$(field_type)"))

    end

    if reduce_grid == :square

        # Reduce the resolution of the result into a new square grid
        # `reduce_factor` here is the factor by wich the number of rows and columns will be reduced
        density = reduceResolution(density ./ physical_factor, reduce_factor)
        x_axis  = reduceTicks(grid.x_ticks, reduce_factor)
        y_axis  = reduceTicks(grid.y_ticks, reduce_factor)

        # The transpose and reverse operation are used to conform to
        # the way `heatmap!` expect the matrix to be structured
        # Depending on the `field_type` and `projection_plane`, different operations
        # are applied to keep the axis consistent between cells and particles
        if field_type == :particles || projection_plane == :xy
            density = reverse!(transpose(density), dims=2)
        elseif projection_plane == :yz
            reverse!(density, dims=1)
        end

    elseif reduce_grid == :circular

        # Reduce the resolution of the result into a circular grid
        # `reduce_factor` here is the number of bins for the circular grid
        density = projectIntoCircularGrid(density ./ physical_factor, reduce_factor)
        x_axis  = [grid.physical_size * (2 * i - 1) / (4 * reduce_factor) for i in 1:reduce_factor]
        y_axis  = x_axis

    else

        throw(ArgumentError("daDensity2DProjection: `reduce_grid` can only be :square or \
        :circular, but I got :$( reduce_grid)"))

    end

    # Set bins with a value of 0 to NaN
    replace!(x -> iszero(x) ? NaN : x, density)

    # Apply log10 to enhance the contrast
    z_axis = log10.(density)

    if logging[]

        log_z_axis = filter(!isnan, z_axis)

        if isempty(log_z_axis)

            min_max_z = (NaN, NaN)
            mean_z    = NaN
            meadian_z = NaN
            mode_z    = NaN

        else

            min_max_z = extrema(log_z_axis)
            mean_z    = mean(log_z_axis)
            meadian_z = median(log_z_axis)
            mode_z    = mode(log_z_axis)

        end

        # Print the density range
        @info(
            "\nDensity range - log₁₀(Σ [$(m_unit * l_unit^-2)]) \
            \n  Simulation: $(basename(filtered_dd[:sim_data].path)) \
            \n  Snapshot:   $(filtered_dd[:snap_data].global_index) \
            \n  Quantity:   $(quantity) \
            \n  Field type: $(field_type) \
            \n  Plane:      $(projection_plane) \
            \n  Min - Max:  $(min_max_z) \
            \n  Mean:       $(mean_z) \
            \n  Median:     $(meadian_z) \
            \n  Mode:       $(mode_z)"
        )

    end

    return x_axis, y_axis, z_axis

end

"""
    daGasSFR2DProjection(
        data_dict::Dict,
        grid::CubicGrid,
        field_type::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Union{Matrix{Float64},Vector{Float64}}}

Project the 3D gas SFR field into a given plane.

!!! note

    If the source of the field are particles, a simple 2D histogram is used. If they are Voronoi cells instead, the SFR of the cells that cross the line of sight of each pixel are added up.

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
  - `grid::CubicGrid`: Cubic grid.
  - `field_type::Symbol`: If the source of the field are `:particles` or Voronoi `:cells`.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection. If `reduce_grid` = :square, the new values will be computed adding up the values of neighboring pixels. `reduce_factor` has to divide the size of `grid` exactly. If `reduce_grid` = :circular, the new values will be computed adding up the values of the pixels the fall within each of the `reduce_factor` concentric rings.
  - `reduce_grid::Symbol=:square`: Type of grid to reduce the resolution of the result. The options are:

      + `:square`    -> The density distribution will be reduced into a regular square grid, with a resolution `reduce_factor` times lower than `grid`. This emulates the way the surface densities are measured in observations. `reduce_factor` = 1 means no reduction in resolution.
      + `:circular` -> The density distribution will be reduced into a flat circular grid, formed by a series of `reduce_factor` concentric rings. This emulates the traditonal way the Kennicutt-Schmidt law is measured in simulations. `reduce_factor` = 1 means that the result will be a single point, the opposite of the `reduce_grid` = :square case.
  - `projection_plane::Symbol=:xy`: Projection plane. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
  - `m_unit::Unitful.Units=u"Msun"`: Mass unit.
  - `l_unit::Unitful.Units=u"kpc"`: Length unit.
  - `t_unit::Unitful.Units=u"yr"`: Time unit.
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

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the ``\\log_{10}`` of the gas SFR at each point of the 2D grid.
"""
function daGasSFR2DProjection(
    data_dict::Dict,
    grid::CubicGrid,
    field_type::Symbol;
    reduce_factor::Int=1,
    reduce_grid::Symbol=:square,
    projection_plane::Symbol=:xy,
    m_unit::Unitful.Units=u"Msun",
    l_unit::Unitful.Units=u"kpc",
    t_unit::Unitful.Units=u"yr",
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Union{Matrix{Float64},Vector{Float64}}}

    filtered_dd = filterData(data_dict; filter_function)

    # Load the cell/particle positions
    positions = filtered_dd[:gas]["POS "]

    # Load the gas SFR
    sfrs = filtered_dd[:gas]["SFR "]

    # If any of the necessary quantities are missing return an empty SFR field
    if any(isempty, [positions, sfrs])
        return grid.x_ticks, grid.y_ticks, fill(NaN, (grid.n_bins, grid.n_bins))
    end

    if field_type == :cells

        # Compute the volume of each cell
        cell_volumes = filtered_dd[:gas]["MASS"] ./ filtered_dd[:gas]["RHO "]

        # Compute the gas SFR densities
        sfr_densities = ustrip.(m_unit * t_unit^-1 * l_unit^-3, sfrs ./ cell_volumes)

        # Allocate memory
        physical_grid = Matrix{Float64}(undef, 3, grid.n_bins^3)

        # Compute the tree for a nearest neighbor search
        kdtree = KDTree(ustrip.(l_unit, positions))

        # Reshape the grid to conform to the way `nn` expect the matrix to be structured
        for i in eachindex(grid.grid)
            physical_grid[1, i] = ustrip(l_unit, grid.grid[i][1])
            physical_grid[2, i] = ustrip(l_unit, grid.grid[i][2])
            physical_grid[3, i] = ustrip(l_unit, grid.grid[i][3])
        end

        # Find the nearest cell to each voxel
        idxs, _ = nn(kdtree, physical_grid)

        # Allocate memory
        sfr_grid = similar(grid.grid, Float64)

        # Load the volume of the voxels
        voxel_volume = ustrip(l_unit^3, grid.bin_volume)

        # Compute the gas SFR in each voxel
        for i in eachindex(grid.grid)
            sfr_grid[i] = sfr_densities[idxs[i]] * voxel_volume
        end

        # Project `sfr_grid` to the target plane
        if projection_plane == :xy
            dims = 3
        elseif projection_plane == :xz
            # Project across dimension 1 to keep it consistent with :xz for `field_type` = :particles
            dims = 1
        elseif projection_plane == :yz
            # Project across dimension 2 to keep it consistent with :yz for `field_type` = :particles
            dims = 2
        else
            throw(ArgumentError("daGasSFR2DProjection: The argument `projection_plane` must be \
            :xy, :xz or :yz, but I got :$(projection_plane)"))
        end

        sfr = dropdims(sum(sfr_grid; dims); dims)

    elseif field_type == :particles

        # Project the particles to the given plane
        if projection_plane == :xy
            pos_2D = positions[[1, 2], :]
        elseif projection_plane == :xz
            pos_2D = positions[[1, 3], :]
        elseif projection_plane == :yz
            pos_2D = positions[[2, 3], :]
        else
            throw(ArgumentError("daGasSFR2DProjection: The argument `projection_plane` must be \
            :xy, :xz or :yz, but I got :$(projection_plane)"))
        end

        # Compute the 2D histogram
        sfr = ustrip.(
            m_unit * t_unit^-1,
            histogram2D(pos_2D, sfrs, flattenGrid(grid); empty_nan=false),
        )

    else

        throw(ArgumentError("daGasSFR2DProjection: The argument `field_type` must be :cells or \
        :particles, but I got :$(field_type)"))

    end

    if reduce_grid == :square

        # Reduce the resolution of the result into a new square grid
        # `reduce_factor` here is the factor by wich the number of rows and columns will be reduced
        sfr = reduceResolution(sfr, reduce_factor; total=true)
        x_axis  = reduceTicks(grid.x_ticks, reduce_factor)
        y_axis  = reduceTicks(grid.y_ticks, reduce_factor)

        # The transpose and reverse operation are used to conform to
        # the way `heatmap!` expect the matrix to be structured
        # Depending on the `field_type` and `projection_plane`, different operations
        # are applied to keep the axis consistent between cells and particles
        if field_type == :particles || projection_plane == :xy
            sfr = reverse!(transpose(sfr), dims=2)
        elseif projection_plane == :yz
            reverse!(sfr, dims=1)
        end

    elseif reduce_grid == :circular

        # Reduce the resolution of the result into a circular grid
        # `reduce_factor` here is the number of bins for the circular grid
        sfr = projectIntoCircularGrid(sfr, reduce_factor; total=true)
        x_axis  = [grid.physical_size * (2 * i - 1) / (4 * reduce_factor) for i in 1:reduce_factor]
        y_axis  = x_axis

    else

        throw(ArgumentError("daGasSFR2DProjection: `reduce_grid` can only be :square or :circular, \
        but I got :$( reduce_grid)"))

    end

    # Set bins with 0 or Inf to NaN
    replace!(x -> (iszero(x) || isinf(x)) ? NaN : x, sfr)

    # Apply log10 to enhance the contrast
    z_axis = log10.(sfr)

    if logging[]

        log_z_axis = filter(!isnan, z_axis)

        if isempty(log_z_axis)

            min_max_z = (NaN, NaN)
            mean_z    = NaN
            meadian_z = NaN
            mode_z    = NaN

        else

            min_max_z = extrema(log_z_axis)
            mean_z    = mean(log_z_axis)
            meadian_z = median(log_z_axis)
            mode_z    = mode(log_z_axis)

        end

        # Print the SFR range
        @info(
            "\nGas SFR range - log₁₀(SFR [$(m_unit * t_unit^-1)]) \
            \n  Simulation: $(basename(filtered_dd[:sim_data].path)) \
            \n  Snapshot:   $(filtered_dd[:snap_data].global_index) \
            \n  Field type: $(field_type) \
            \n  Plane:      $(projection_plane) \
            \n  Min - Max:  $(min_max_z) \
            \n  Mean:       $(mean_z) \
            \n  Median:     $(meadian_z) \
            \n  Mode:       $(mode_z)"
        )

    end

    return x_axis, y_axis, z_axis

end

"""
    daMetallicity2DProjection(
        data_dict::Dict,
        grid::CubicGrid,
        component::Symbol,
        field_type::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Union{Matrix{Float64},Vector{Float64}}}

Project the 3D metallicity field to a given plane.

!!! note

    The metallicity in each pixel is the total metal mass divided by the total gas mass, in the column given by that pixel. By default, the total metallicity (`element` = :all) is given in solar units.

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
  - `grid::CubicGrid`: Cubic grid.
  - `component::Symbol`: Target component. It can be either `:stars` or `:gas`.
  - `field_type::Symbol`: If the source of the field are `:particles` or Voronoi `:cells`.
  - `element::Symbol=:all`: Target element. The possibilities are the keys of [`ELEMENT_INDEX`](@ref). Set it to :all if you want the total metallicity.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection. If `reduce_grid` = :square, the new values will be computed adding the values of mass in each neighboring pixel. `reduce_factor` has to divide the size of `grid` exactly. If `reduce_grid` = :circular, the new values will be computed adding up the values of mass of the pixels the fall within each of the `reduce_factor` concentric rings.
  - `reduce_grid::Symbol=:square`: Type of grid to reduce the resolution of the result. The options are:

      + `:square`    -> The density distribution will be reduced into a regular square grid, with a resolution `reduce_factor` times lower than `grid`. This emulates the way the surface densities are measured in observations. `reduce_factor` = 1 means no reduction in resolution.
      + `:circular` -> The density distribution will be reduced into a flat circular grid, formed by a series of `reduce_factor` concentric rings. This emulates the traditonal way the Kennicutt-Schmidt law is measured in simulations. `reduce_factor` = 1 means that the result will be a single point, the opposite of the `reduce_grid` = :square case.
  - `projection_plane::Symbol=:xy`: Projection plane. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
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

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the ``\\log_{10}`` of the metallicity at each point of the 2D grid.
"""
function daMetallicity2DProjection(
    data_dict::Dict,
    grid::CubicGrid,
    component::Symbol,
    field_type::Symbol;
    element::Symbol=:all,
    reduce_factor::Int=1,
    reduce_grid::Symbol=:square,
    projection_plane::Symbol=:xy,
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Union{Matrix{Float64},Vector{Float64}}}

    filtered_dd = filterData(data_dict; filter_function)

    (
        component ∈ [:gas, :stars] ||
        throw(ArgumentError("daMetallicity2DProjection: The argument `component` must be :gas \
        or :stars, but I got :$(component)"))
    )

    (
        element ∈ [:all, keys(ELEMENT_INDEX)...] ||
        throw(ArgumentError("daMetallicity2DProjection: The argument `element` can only be :all \
        or one of the keys of `ELEMENT_INDEX` (see `./src/constants/globals.jl`), \
        but I got :$(quantity)"))
    )

    # Load the cell/particle positions
    positions = filtered_dd[component]["POS "]

    # If the necessary quantities are missing return an empty metallicity field
    if isempty(positions)
        return grid.x_ticks, grid.y_ticks, fill(NaN, (grid.n_bins, grid.n_bins))
    end

    if field_type == :cells

        # Load the cell densities
        cell_densities = filtered_dd[component]["RHO "]

        # Compute the volume of each cell
        cell_volumes = filtered_dd[component]["MASS"] ./ cell_densities

        if element == :all

            # Compute the metal mass density in each cell
            metal_densities = computeMetalMass(filtered_dd, component) ./ cell_volumes

        else

            # Compute the `element` mass density in each cell
            metal_densities  = computeElementMass(filtered_dd, component, element) ./ cell_volumes

            # Compute the `hydrogen` mass density in each cell
            hydrogen_density = computeElementMass(filtered_dd, component, :H) ./ cell_volumes

        end

        # Allocate memory
        physical_grid = Matrix{Float64}(undef, 3, grid.n_bins^3)

        # Reshape the grid to conform to the way `nn` expect the matrix to be structured
        for i in eachindex(grid.grid)
            physical_grid[1, i] = ustrip(u"kpc", grid.grid[i][1])
            physical_grid[2, i] = ustrip(u"kpc", grid.grid[i][2])
            physical_grid[3, i] = ustrip(u"kpc", grid.grid[i][3])
        end

        # Compute the tree for a nearest neighbor search
        kdtree = KDTree(ustrip.(u"kpc", positions))

        # Find the nearest cell to each voxel
        idxs, _ = nn(kdtree, physical_grid)

        # Allocate memory
        norm_mass_grid  = similar(grid.grid, Float64)
        metal_mass_grid = similar(grid.grid, Float64)

        # Compute the corresponding masses in each voxel
        for i in eachindex(grid.grid)

            if element == :all
                norm_mass_grid[i]   = ustrip.(u"Msun", cell_densities[idxs[i]] * grid.bin_volume)
            else
                norm_mass_grid[i]   = ustrip.(u"Msun", hydrogen_density[idxs[i]] * grid.bin_volume)
            end

            metal_mass_grid[i] = ustrip.(u"Msun", metal_densities[idxs[i]] * grid.bin_volume)

        end

        # Project `metal_mass_grid` to the target plane
        if projection_plane == :xy
            dims = 3
        elseif projection_plane == :xz
            # Project across dimension 1 to keep it consistent with :xz for `field_type` = :particles
            dims = 1
        elseif projection_plane == :yz
            # Project across dimension 2 to keep it consistent with :yz for `field_type` = :particles
            dims = 2
        else
            throw(ArgumentError("daMetallicity2DProjection: The argument `projection_plane` must \
            be :xy, :xz or :yz, but I got :$(projection_plane)"))
        end

        metal_mass = dropdims(sum(metal_mass_grid; dims); dims)
        norm_mass  = dropdims(sum(norm_mass_grid; dims); dims)

    elseif field_type == :particles

        # Project the particles to the given plane
        if projection_plane == :xy
            pos_2D = positions[[1, 2], :]
        elseif projection_plane == :xz
            pos_2D = positions[[1, 3], :]
        elseif projection_plane == :yz
            pos_2D = positions[[2, 3], :]
        else
            throw(ArgumentError("daMetallicity2DProjection: The argument `projection_plane` must \
            be :xy, :xz or :yz, but I got :$(projection_plane)"))
        end

        # Compute the metal mass 2D histogram
        metal_mass = histogram2D(
            pos_2D,
            computeMetalMass(filtered_dd, component),
            flattenGrid(grid);
            empty_nan=true,
        )

        # Compute the normalization mass 2D histogram
        norm_mass = histogram2D(
            pos_2D,
            filtered_dd[component]["MASS"],
            flattenGrid(grid);
            empty_nan=true,
        )

    else

        throw(ArgumentError("daMetallicity2DProjection: The argument `field_type` must be :cells \
        or :particles, but I got :$(field_type)"))

    end

    if reduce_grid == :square

        # Reduce the resolution of the result into a new square grid
        # `reduce_factor` here is the factor by wich the number of rows and columns will be reduced
        metal_mass = reduceResolution(metal_mass, reduce_factor; total=true)
        norm_mass  = reduceResolution(norm_mass, reduce_factor; total=true)

        x_axis = reduceTicks(grid.x_ticks, reduce_factor)
        y_axis = reduceTicks(grid.y_ticks, reduce_factor)

        # Compute the metallicity in each bin
        metallicity = uconvert.(Unitful.NoUnits, metal_mass ./ norm_mass)

        # The transpose and reverse operation are used to conform to
        # the way `heatmap!` expect the matrix to be structured
        # Depending on the `field_type` and `projection_plane`, different operations
        # are applied to keep the axis consistent between cells and particles
        if field_type == :particles || projection_plane == :xy
            metallicity = reverse!(transpose(metallicity), dims=2)
        elseif projection_plane == :yz
            reverse!(metallicity, dims=1)
        end

    elseif reduce_grid == :circular

        # Reduce the resolution of the result into a circular grid
        # `reduce_factor` here is the number of bins for the circular grid
        metal_mass = projectIntoCircularGrid(metal_mass, reduce_factor; total=true)
        norm_mass  = projectIntoCircularGrid(norm_mass, reduce_factor; total=true)

        x_axis = [grid.physical_size * (2 * i - 1) / (4 * reduce_factor) for i in 1:reduce_factor]
        y_axis = x_axis

        # Compute the metallicity in each bin
        metallicity = uconvert.(Unitful.NoUnits, metal_mass ./ norm_mass)

    else

        throw(ArgumentError("daMetallicity2DProjection: `reduce_grid` can only be :square or \
        :circular, but I got :$(reduce_grid)"))

    end

    # Set bins with 0 or Inf to NaN
    replace!(x -> (iszero(x) || isinf(x)) ? NaN : x, metallicity)

    # Apply log10 to enhance the contrast
    if element == :all
        z_axis = log10.(metallicity ./ SOLAR_METALLICITY)
    else
        # Add 12 so the result is by convention 12 + log10(X / H)
        z_axis = 12 .+ log10.(metallicity)
    end

    if logging[]

        clean_z_axis = filter(!isnan, z_axis)

        if isempty(clean_z_axis)

            min_max_z = (NaN, NaN)
            mean_z    = NaN
            meadian_z = NaN
            mode_z    = NaN

        else

            min_max_z = extrema(clean_z_axis)
            mean_z    = mean(clean_z_axis)
            meadian_z = median(clean_z_axis)
            mode_z    = mode(clean_z_axis)

        end

        @info(
            "\nMetallicity range - $(element == :all ? "log₁₀(Z [Z⊙])" : "12 + log₁₀($(element) / H)") \
            \n  Simulation: $(basename(filtered_dd[:sim_data].path)) \
            \n  Snapshot:   $(filtered_dd[:snap_data].global_index) \
            \n  Component:  $(component) \
            \n  Field type: $(field_type) \
            \n  Plane:      $(projection_plane)\
            \n  Min - Max:  $(min_max_z) \
            \n  Mean:       $(mean_z) \
            \n  Median:     $(meadian_z) \
            \n  Mode:       $(mode_z)"
        )

    end

    return x_axis, y_axis, z_axis

end

"""
    daTemperature2DProjection(
        data_dict::Dict,
        grid::CubicGrid,
        field_type::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Union{Matrix{Float64},Vector{Float64}}}

Project the 3D temperature field to a given plane.

!!! note

    The temperature in each pixel is the mean temperature of the column given by that pixel. By default, ``K`` is used as unit of temperature, so the output will be ``\\log_{10}(T \\, [\\mathrm{K}])``.

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
  - `grid::CubicGrid`: Cubic grid.
  - `field_type::Symbol`: If the source of the field are `:particles` or Voronoi `:cells`.
  - `reduce_factor::Int=1`: Factor by which the resolution of the result will be reduced. This will be applied after the density projection. If `reduce_grid` = :square, the new values will be computed averaging the values of neighboring pixels. `reduce_factor` has to divide the size of `grid` exactly. If `reduce_grid` = :circular, the new values will be computed averaging the values of the pixels the fall within each of the `reduce_factor` concentric rings.
  - `reduce_grid::Symbol=:square`: Type of grid to reduce the resolution of the result. The options are:

      + `:square`    -> The density distribution will be reduced into a regular square grid, with a resolution `reduce_factor` times lower than `grid`. This emulates the way the surface densities are measured in observations. `reduce_factor` = 1 means no reduction in resolution.
      + `:circular` -> The density distribution will be reduced into a flat circular grid, formed by a series of `reduce_factor` concentric rings. This emulates the traditonal way the Kennicutt-Schmidt law is measured in simulations. `reduce_factor` = 1 means that the result will be a single point, the opposite of the `reduce_grid` = :square case.
  - `projection_plane::Symbol=:xy`: Projection plane. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
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

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the values of temperature at each grid point.
"""
function daTemperature2DProjection(
    data_dict::Dict,
    grid::CubicGrid,
    field_type::Symbol;
    reduce_factor::Int=1,
    reduce_grid::Symbol=:square,
    projection_plane::Symbol=:xy,
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Union{Matrix{Float64},Vector{Float64}}}

    filtered_dd = filterData(data_dict; filter_function)

    # Load the cell/particle temperatures
    temperatures = ustrip.(u"K", filtered_dd[:gas]["TEMP"])

    # Load the cell/particle positions
    positions = filtered_dd[:gas]["POS "]

    # If any of the necessary quantities are missing return an empty temperature field
    if any(isempty, [temperatures, positions])
        return grid.x_ticks, grid.y_ticks, fill(NaN, (grid.n_bins, grid.n_bins))
    end

    # Set the units
    m_unit = u"Msun"
    l_unit = u"kpc"

    if field_type == :cells

        # Load the gas densities
        densities = ustrip.(m_unit * l_unit^-3, filtered_dd[:gas]["RHO "])

        # Load the volume and area of the voxels
        voxel_volume = ustrip(l_unit^3, grid.bin_volume)

        # Allocate memory
        physical_grid = Matrix{Float64}(undef, 3, grid.n_bins^3)

        # Reshape the grid to conform to the way `nn` expect the matrix to be structured
        for i in eachindex(grid.grid)
            physical_grid[1, i] = ustrip(u"kpc", grid.grid[i][1])
            physical_grid[2, i] = ustrip(u"kpc", grid.grid[i][2])
            physical_grid[3, i] = ustrip(u"kpc", grid.grid[i][3])
        end

        # Compute the tree for a nearest neighbor search
        kdtree = KDTree(ustrip.(u"kpc", positions))

        # Find the nearest cell to each voxel
        idxs, _ = nn(kdtree, physical_grid)

        # Allocate memory
        temperature_grid = similar(grid.grid, Float64)
        mass_grid        = similar(grid.grid, Float64)

        # Compute the temperature and mass of each voxel
        for i in eachindex(grid.grid)
            temperature_grid[i] = temperatures[idxs[i]]
            mass_grid[i]        = densities[idxs[i]] * voxel_volume
        end

        # Use the mass of each voxel as a weight
        weighted_temperature = temperature_grid .* mass_grid

        # Project `temperature_grid` to the target plane
        if projection_plane == :xy
            dims = 3
        elseif projection_plane == :xz
            # Project across dimension 1 to keep it consistent with :xz for `field_type` = :particles
            dims = 1
        elseif projection_plane == :yz
            # Project across dimension 2 to keep it consistent with :yz for `field_type` = :particles
            dims = 2
        else
            throw(ArgumentError("daTemperature2DProjection: The argument `projection_plane` must \
            be :xy, :xz or :yz, but I got :$(projection_plane)"))
        end

        normalization = dropdims(sum(mass_grid; dims); dims)
        temperature = dropdims(sum(weighted_temperature; dims); dims) ./ normalization

    elseif field_type == :particles

        # Project the particles to the given plane
        if projection_plane == :xy
            pos_2D = positions[[1, 2], :]
        elseif projection_plane == :xz
            pos_2D = positions[[1, 3], :]
        elseif projection_plane == :yz
            pos_2D = positions[[2, 3], :]
        else
            throw(ArgumentError("daTemperature2DProjection: The argument `projection_plane` must \
            be :xy, :xz or :yz, but I got :$(projection_plane)"))
        end

        # Compute the 2D histogram
        temperature = histogram2D(
            pos_2D,
            temperatures,
            flattenGrid(grid);
            total=false,
            empty_nan=false,
        )

    else

        throw(ArgumentError("daTemperature2DProjection: The argument `field_type` must be :cells \
        or :particles, but I got :$(field_type)"))

    end

    if reduce_grid == :square

        # Reduce the resolution of the result into a new square grid
        # `reduce_factor` here is the factor by wich the number of rows and columns will be reduced
        temperature = reduceResolution(temperature, reduce_factor)
        x_axis  = reduceTicks(grid.x_ticks, reduce_factor)
        y_axis  = reduceTicks(grid.y_ticks, reduce_factor)

        # The transpose and reverse operation are used to conform to
        # the way `heatmap!` expect the matrix to be structured
        # Depending on the `field_type` and `projection_plane`, different operations
        # are applied to keep the axis consistent between cells and particles
        if field_type == :particles || projection_plane == :xy
            temperature = reverse!(transpose(temperature), dims=2)
        elseif projection_plane == :yz
            reverse!(temperature, dims=1)
        end

    elseif reduce_grid == :circular

        # Reduce the resolution of the result into a circular grid
        # `reduce_factor` here is the number of bins for the circular grid
        temperature = projectIntoCircularGrid(temperature, reduce_factor)
        x_axis  = [grid.physical_size * (2 * i - 1) / (4 * reduce_factor) for i in 1:reduce_factor]
        y_axis  = x_axis

    else

        throw(ArgumentError("daTemperature2DProjection: `reduce_grid` can only be :square or \
        :circular, but I got :$( reduce_grid)"))

    end

    # Set bins with a value of 0 to NaN
    replace!(x -> iszero(x) ? NaN : x, temperature)

    # Apply log10 to enhance the contrast
    z_axis = log10.(temperature)

    if logging[]

        clean_z_axis = filter(!isinf, z_axis)

        if isempty(clean_z_axis)

            min_max_T = (NaN, NaN)
            mean_T    = NaN
            meadian_T = NaN
            mode_T    = NaN

        else

            min_max_T = extrema(clean_z_axis)
            mean_T    = mean(clean_z_axis)
            meadian_T = median(clean_z_axis)
            mode_T    = mode(clean_z_axis)

        end

        # Print the temperature range
        @info(
            "\nTemperature range - log₁₀(T [K]) \
            \n  Simulation: $(basename(filtered_dd[:sim_data].path)) \
            \n  Snapshot:   $(filtered_dd[:snap_data].global_index) \
            \n  Field type: $(field_type) \
            \n  Type:       $(type) \
            \n  Min - Max:  $(min_max_T) \
            \n  Mean:       $(mean_T) \
            \n  Median:     $(meadian_T) \
            \n  Mode:       $(mode_T)"
        )

    end

    return x_axis, y_axis, z_axis

end

"""
    daScatterDensity(
        data_dict::Dict,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Number},Vector{<:Number},Matrix{Float64}}

Turn a scatter plot into a 2D histogram.

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
  - `x_range::Union{NTuple{2,<:Number},Nothing}=nothing`: x axis range for the histogram grid. If set to `nothing`, the extrema of the values will be used.
  - `y_range::Union{NTuple{2,<:Number},Nothing}=nothing`: y axis range for the histogram grid. If set to `nothing`, the extrema of the values will be used.
  - `x_log::Union{Unitful.Units,Nothing}=nothing`: Desired unit of `x_quantity`, if you want to use log10(`x_quantity`) for the x axis.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Desired unit of `y_quantity`, if you want to use log10(`y_quantity`) for the y axis.
  - `n_bins::Int=100`: Number of bins per side of the grid.
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

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the counts.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function daScatterDensity(
    data_dict::Dict,
    x_quantity::Symbol,
    y_quantity::Symbol;
    x_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    y_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    x_log::Union{Unitful.Units,Nothing}=nothing,
    y_log::Union{Unitful.Units,Nothing}=nothing,
    n_bins::Int=100,
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Number},Vector{<:Number},Matrix{Float64}}

    filtered_dd = filterData(data_dict; filter_function)

    # Compute the values of the quantities
    x_values = scatterQty(filtered_dd, x_quantity)
    y_values = scatterQty(filtered_dd, y_quantity)

    # If any of the necessary quantities are missing return an empty histogram
    if any(isempty, [x_values, y_values])
        return 1:n_bins, 1:n_bins, fill(NaN, (n_bins, n_bins))
    end

    (
        length(x_values) == length(y_values) ||
        throw(ArgumentError("daScatterDensity: :$(x_quantity) and :$(y_quantity) have a diferent \
        number of values. They should be the same"))
    )

    # If requested, apply log10 to the x axis data, ignoring 0 values
    if !isnothing(x_log)
        null_x_idxs = findall(iszero, x_values)
        x_values    = log10.(deleteat!(ustrip.(x_log, x_values), null_x_idxs))
        y_values    = deleteat!(y_values, null_x_idxs)
    end

    # If requested, apply log10 to the y axis data, ignoring 0 values
    if !isnothing(y_log)
        null_y_idxs = findall(iszero, y_values)
        x_values    = deleteat!(x_values, null_y_idxs)
        y_values    = log10.(deleteat!(ustrip.(y_log, y_values), null_y_idxs))
    end

    # If there is no specified range, use the extrema of the x values
    if isnothing(x_range)
        x_range = extrema(x_values)
    end

    # If there is no specified range, use the extrema of the y values
    if isnothing(y_range)
        y_range = extrema(y_values)
    end

    # Compute the bin half width for each axis
    x_bin_h_width = 0.5 * (x_range[2] - x_range[1]) / n_bins
    y_bin_h_width = 0.5 * (y_range[2] - y_range[1]) / n_bins

    # Compute the center value of each bin for each axis
    x_axis = collect(range(x_range[1] + x_bin_h_width; length=n_bins, step=2 * x_bin_h_width))
    y_axis = collect(range(y_range[1] + y_bin_h_width; length=n_bins, step=2 * y_bin_h_width))

    # Compute the 2D histogram
    counts = Float64.(histogram2D(
        permutedims(hcat(x_values, y_values), (2, 1)),
        collect(range(x_range[1], x_range[2]; length=n_bins + 1)),
        collect(range(y_range[1], y_range[2]; length=n_bins + 1)),
    ))

    # Set bins with a value of 0 to NaN
    replace!(x -> iszero(x) ? NaN : x, counts)

    # The transpose and reverse operation are used to conform to the way heatmap! expect the matrix to be structured,
    # and log10 is used to enhance the contrast
    z_axis = reverse!(transpose(log10.(counts)), dims=2)

    return x_axis, y_axis, z_axis

end

"""
    daScatterWeightedDensity(
        data_dict::Dict,
        x_quantity::Symbol,
        y_quantity::Symbol,
        z_quantity::Symbol,
        z_unit::Unitful.Units;
        <keyword arguments>
    )::Tuple{Vector{<:Number},Vector{<:Number},Matrix{Float64}}

Turn a scatter plot into a 2D histogram, weighted by `z_quantity`.

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
  - `x_range::Union{NTuple{2,<:Number},Nothing}=nothing`: x axis range for the histogram grid. If set to `nothing`, the extrema of the values will be used.
  - `y_range::Union{NTuple{2,<:Number},Nothing}=nothing`: y axis range for the histogram grid. If set to `nothing`, the extrema of the values will be used.
  - `x_log::Union{Unitful.Units,Nothing}=nothing`: Desired unit of `x_quantity`, if you want to use log10(`x_quantity`) for the x axis.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Desired unit of `y_quantity`, if you want to use log10(`y_quantity`) for the y axis.
  - `total::Bool=true`: If the sum (default) or the mean of `z_quantity` will be used as the value of each pixel.
  - `n_bins::Int=100`: Number of bins per side of the grid.
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

  - A tuple with three elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the counts.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function daScatterWeightedDensity(
    data_dict::Dict,
    x_quantity::Symbol,
    y_quantity::Symbol,
    z_quantity::Symbol,
    z_unit::Unitful.Units;
    x_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    y_range::Union{NTuple{2,<:Number},Nothing}=nothing,
    x_log::Union{Unitful.Units,Nothing}=nothing,
    y_log::Union{Unitful.Units,Nothing}=nothing,
    total::Bool=true,
    n_bins::Int=100,
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Number},Vector{<:Number},Matrix{Float64}}

    filtered_dd = filterData(data_dict; filter_function)

    # Compute the values of the quantities
    x_values = scatterQty(filtered_dd, x_quantity)
    y_values = scatterQty(filtered_dd, y_quantity)
    z_values = scatterQty(filtered_dd, z_quantity)

    # If any of the necessary quantities are missing return an empty histogram
    if any(isempty, [x_values, y_values, z_values])
        return collect(1:n_bins), collect(1:n_bins), fill(NaN, (n_bins, n_bins))
    end

    (
        allequal(length, [x_values, y_values, z_values]) ||
        throw(ArgumentError("daScatterWeightedDensity: :$(x_quantity), :$(y_quantity), \
        and :$(z_quantity) have a diferent number of values. They should be the same"))
    )

    # Delete NaN values
    nan_x_idxs = map(isnan, x_values)
    nan_y_idxs = map(isnan, y_values)
    nan_z_idxs = map(isnan, z_values)
    nan_idxs   = nan_x_idxs ∪ nan_y_idxs ∪ nan_z_idxs
    deleteat!(x_values, nan_idxs)
    deleteat!(y_values, nan_idxs)
    deleteat!(z_values, nan_idxs)

    # If requested, apply log10 to the x axis data, ignoring 0 values
    if !isnothing(x_log)

        null_x_idxs = map(iszero, x_values)
        x_values    = log10.(deleteat!(ustrip.(x_log, x_values), null_x_idxs))
        deleteat!(y_values, null_x_idxs)
        deleteat!(z_values, null_x_idxs)

    end

    # If requested, apply log10 to the y axis data, ignoring 0 values
    if !isnothing(y_log)

        null_y_idxs = findall(iszero, y_values)
        y_values    = log10.(deleteat!(ustrip.(y_log, y_values), null_y_idxs))
        deleteat!(x_values, null_y_idxs)
        deleteat!(z_values, null_y_idxs)

    end

    # If there is no specified range, use the extrema of the x values
    if isnothing(x_range)
        x_range = extrema(x_values)
    end

    # If there is no specified range, use the extrema of the y values
    if isnothing(y_range)
        y_range = extrema(y_values)
    end

    # Compute the bin half width for each axis
    x_bin_h_width = 0.5 * (x_range[2] - x_range[1]) / n_bins
    y_bin_h_width = 0.5 * (y_range[2] - y_range[1]) / n_bins

    # Compute the center value of each bin for each axis
    x_axis = collect(range(x_range[1] + x_bin_h_width; length=n_bins, step=2 * x_bin_h_width))
    y_axis = collect(range(y_range[1] + y_bin_h_width; length=n_bins, step=2 * y_bin_h_width))

    # Compute the 2D histogram
    values = log10.(ustrip.(z_unit, histogram2D(
        permutedims(hcat(x_values, y_values), (2, 1)),
        z_values,
        collect(range(x_range[1], x_range[2]; length=n_bins + 1)),
        collect(range(y_range[1], y_range[2]; length=n_bins + 1));
        total,
    )))

    if logging[]

        clean_c = filter(!isnan, values)

        if isempty(clean_c)

            min_max_c = (NaN, NaN)
            mean_c    = NaN
            meadian_c = NaN
            mode_c    = NaN

        else

            min_max_c = extrema(clean_c)
            mean_c    = mean(clean_c)
            meadian_c = median(clean_c)
            mode_c    = mode(clean_c)

        end

        # Print the color range
        @info(
            "\nColor range \
            \n  Simulation: $(basename(filtered_dd[:sim_data].path)) \
            \n  Snapshot:   $(filtered_dd[:snap_data].global_index) \
            \n  Quantity:   $(z_quantity) \
            \n  Min - Max:  $(min_max_c) \
            \n  Mean:       $(mean_c) \
            \n  Median:     $(meadian_c) \
            \n  Mode:       $(mode_c)"
        )

    end

    # Set bins with a value of 0 to NaN
    replace!(x -> iszero(x) ? NaN : x, values)

    # The transpose and reverse operation are used to conform to the way heatmap! expect the matrix to be structured,
    # and log10 is used to enhance the contrast
    z_axis = reverse!(transpose(values), dims=2)

    return x_axis, y_axis, z_axis

end

"""
    daVelocityField(
        data_dict::Dict,
        grid::SquareGrid,
        component::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{<:Number},Matrix{<:Number}}

Compute a 2D mean velocity field.

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
  - `grid::SquareGrid`: Square grid.
  - `component::Symbol`: For which cell/particle type the velocity field will be computed. The possibilities are the keys of [`PARTICLE_INDEX`](@ref).
  - `projection_plane::Symbol=:xy`: Projection plane. The options are `:xy`, `:xz`, and `:yz`. The disk is generally oriented to have its axis of rotation parallel to the z axis.
  - `velocity_units::Bool=false`: If the velocity will be given as an `Unitful.Quantity` with units or as a `Flot64` (in which case the underlying unit is ``\\mathrm{km} \\, \\mathrm{s}^{-1}``).
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

  - A tuple with four elements:

      + A vector with the x coordinates of the grid.
      + A vector with the y coordinates of the grid.
      + A matrix with the mean velocity in the x direction at each grid point.
      + A matrix with the mean velocity in the y direction at each grid point.
"""
function daVelocityField(
    data_dict::Dict,
    grid::SquareGrid,
    component::Symbol;
    projection_plane::Symbol=:xy,
    velocity_units::Bool=false,
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{<:Unitful.Length},Matrix{<:Number},Matrix{<:Number}}

    filtered_dd = filterData(data_dict; filter_function)

    positions  = filtered_dd[component]["POS "]
    velocities = filtered_dd[component]["VEL "]

    # If any of the necessary quantities are missing return an empty velocity field
    if any(isempty, [positions, velocities])
        return grid.x_ticks, grid.y_ticks, zeros(size(grid.grid)), zeros(size(grid.grid))
    end

    # Project the cell/particles to the chosen plane
    if projection_plane == :xy

        pos_2D = positions[[1, 2], :]

        # Compute the components of the mean velocity
        vx = histogram2D(pos_2D, vec(velocities[1, :]), grid; total=false)
        vy = histogram2D(pos_2D, vec(velocities[2, :]), grid; total=false)

    elseif projection_plane == :xz

        pos_2D = positions[[1, 3], :]

        # Compute the components of the mean velocity
        vx = histogram2D(pos_2D, vec(velocities[1, :]), grid; total=false)
        vy = histogram2D(pos_2D, vec(velocities[3, :]), grid; total=false)

    elseif projection_plane == :yz

        pos_2D = positions[[2, 3], :]

        # Compute the components of the mean velocity
        vx = histogram2D(pos_2D, vec(velocities[2, :]), grid; total=false)
        vy = histogram2D(pos_2D, vec(velocities[3, :]), grid; total=false)

    else

        throw(ArgumentError("daVelocityField: The argument `projection_plane` must be \
        :xy, :xz or :yz, but I got :$(projection_plane)"))

    end

    # The transpose and reverse operation are used to conform to the way arrows! expect the matrix to be structured
    vx = collect(reverse!(transpose(vx), dims=2))
    vy = collect(reverse!(transpose(vy), dims=2))

    if !velocity_units
        vx = ustrip.(u"km*s^-1", vx)
        vy = ustrip.(u"km*s^-1", vy)
    end

    return grid.x_ticks, grid.y_ticks, vx, vy

end

"""
    daIntegrateGalaxy(
        data_dict::Dict,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute two global quantities of the simulation.

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
      + `:ode_stellar_gas_P`           -> Gas pressure, for the gas that form the stars.
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

  - A tuple with two elements:

      + A single element vector with the value of `x_quantity`.
      + A single element vector with the value of `y_quantity`.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function daIntegrateGalaxy(
    data_dict::Dict,
    x_quantity::Symbol,
    y_quantity::Symbol;
    filter_function::Function=filterNothing,
)::NTuple{2,Vector{<:Number}}

    filtered_dd = filterData(data_dict; filter_function)

    return [integrateQty(filtered_dd, x_quantity)], [integrateQty(filtered_dd, y_quantity)]

end

"""
    daScatterGalaxy(
        data_dict::Dict,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute two quantities for every cell/particle in the simulation.

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
  - `x_log::Union{Unitful.Units,Nothing}=nothing`: Desired unit of `x_quantity`, if you want to use log10(`x_quantity`) for the x axis.
  - `y_log::Union{Unitful.Units,Nothing}=nothing`: Desired unit of `y_quantity`, if you want to use log10(`y_quantity`) for the y axis.
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

  - A tuple with two elements:

      + A vector with the values of `x_quantity`.
      + A vector with the values of `y_quantity`.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function daScatterGalaxy(
    data_dict::Dict,
    x_quantity::Symbol,
    y_quantity::Symbol;
    x_log::Union{Unitful.Units,Nothing}=nothing,
    y_log::Union{Unitful.Units,Nothing}=nothing,
    filter_function::Function=filterNothing,
)::NTuple{2,Vector{<:Number}}

    filtered_dd = filterData(data_dict; filter_function)

    x_values = scatterQty(filtered_dd, x_quantity)
    y_values = scatterQty(filtered_dd, y_quantity)

    (
        length(x_values) == length(y_values) ||
        throw(ArgumentError("daScatterGalaxy: :$(x_quantity) and :$(y_quantity) have a diferent \
        number of values. They should be the same"))
    )

    if isnothing(x_log)
        x_idxs = Int64[]
    else
        x_idxs = map(iszero, x_values)
    end

    if isnothing(y_log)
        y_idxs = Int64[]
    else
        y_idxs = map(iszero, y_values)
    end

    delete_idxs = x_idxs ∪ y_idxs

    deleteat!(x_values, delete_idxs)
    deleteat!(y_values, delete_idxs)

    if isnothing(x_log)
        x_axis = x_values
    else
        x_axis = log10.(ustrip.(x_log, x_values))
    end

    if isnothing(y_log)
        y_axis = y_values
    else
        y_axis = log10.(ustrip.(y_log, y_values))
    end

    if any(isempty, [x_axis, y_axis])
        return Float64[], Float64[]
    end

    return x_axis, y_axis

end

"""
    daGasFractions(
        data_dict::Dict,
        quantity::Symbol,
        edges::Vector{<:Number};
        <keyword arguments>
    )::Union{NTuple{2,Vector{<:Number}},Nothing}

Compute the values for a bar plot of the gas fractions, where the bins are a given gas `quantity`.

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
  - `quantity::Symbol`: Target quantity for the bins. The possibilities are:

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
      + `:gas_metallicity`             -> Mass fraction of all elements above He in the gas (solar units).
      + `:X_gas_abundance`             -> Gas abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. The possibilities are the keys of [`ELEMENT_INDEX`](@ref).
      + `:gas_radial_distance`         -> Distance of every gas cell to the origin.
      + `:gas_xy_distance`             -> Projected distance of every gas cell to the origin.
      + `:gas_sfr`                     -> SFR associated to each gas particle/cell within the code.
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

  - A tuple with two elements:

      + A vector with the positions of each bar.
      + A vector with the height of each bar.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function daGasFractions(
    data_dict::Dict,
    quantity::Symbol,
    edges::Vector{<:Number};
    filter_function::Function=filterNothing,
)::Union{NTuple{2,Vector{<:Number}},Nothing}

    filtered_dd = filterData(data_dict; filter_function)

    # Compute the number of bins for the gas quantity
    n_bins = length(edges) - 1

    # Compute the number of bars per bin
    n_bars = 3

    # Compute the gas quantities
    gas_qty  = scatterQty(filtered_dd, quantity)

    # If any of the necessary quantities are missing return nothing
    !isempty(gas_qty) || return nothing

    # Compute the mass of each gas phase
    ionized_mass   = computeMass(filtered_dd, :ionized)
    atomic_mass    = computeMass(filtered_dd, :atomic)
    molecular_mass = computeMass(filtered_dd, :molecular)

    # If any of the necessary quantities are missing return nothing
    !any(isempty, [ionized_mass, atomic_mass, molecular_mass]) || return nothing

    # Allocate memory
    percents = Vector{Float64}(undef, n_bins * n_bars)

    # Loop over each bin
    for i in 1:n_bins

        # Compute the indices of the gas cells inside the current bin
        qty_idx = findall(q->edges[i] < q <= edges[i + 1], gas_qty)

        if isempty(qty_idx)

            percents[((i - 1) * n_bars + 1):(i * n_bars)] .= 0.0

        else

            total_masses = [
                sum(molecular_mass[qty_idx]),
                sum(atomic_mass[qty_idx]),
                sum(ionized_mass[qty_idx]),
            ]

            # Compute the mass fraction of each phase inside the current bin
            fractions = uconvert.(Unitful.NoUnits, (total_masses ./ sum(total_masses)) .* 100)

            percents[((i - 1) * n_bars + 1):(i * n_bars)] .= fractions

        end

    end

    return repeat(1:n_bins, inner=n_bars), percents

end

"""
    daStellarMetallictyHistogram(data_dict::Dict)::Union{Tuple{Vector{Float64}},Nothing}

Compute the stellar metallicity (in [`SOLAR_METALLICITY`](@ref) units), for an histogram.

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

  - A Tuple with one elements:

      + A Vector with the stellar metallicites.
"""
function daStellarMetallictyHistogram(data_dict::Dict)::Union{Tuple{Vector{Float64}},Nothing}

    metallicity = scatterQty(data_dict, :stellar_metallicity)

    !isempty(metallicity) || return nothing

    return (metallicity,)

end

"""
    daStellarBTHistogram(data_dict::Dict)::Union{Tuple{Vector{<:Unitful.Time}},Nothing}

Compute the stellar birth times, for an histogram.

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

  - A Tuple with one elements:

      + A Vector with the birth times.
"""
function daStellarBTHistogram(data_dict::Dict)::Union{Tuple{Vector{<:Unitful.Time}},Nothing}

    birth_ticks = data_dict[:stars]["GAGE"]

    !isempty(birth_ticks) || return nothing

    if data_dict[:sim_data].cosmological
        # Go from scale factor to physical time
        birth_times = computeTime(birth_ticks, data_dict[:snap_data].header)
    else
        birth_times = birth_ticks
    end

    return (birth_times,)

end

####################################################################################################
# Signature for the plotTimeSeries function in ./src/plotting/pipelines.jl
####################################################################################################
#
# A data analysis functions for plotTimeSeries must take a Simulation struct, and return two
# vectors. It should return `nothing` if the input data has some problem that prevents computation
# (e.g. is empty).
#
# Expected signature:
#
#   da_function(sim_data, args...; kw_args...) -> (processed_data_x, processed_data_y)
#
# where:
#
#   - sim_data::Simulation, see the Simulation struct in ./src/constants/globals.jl
#   - processed_data_x::Vector{<:Number}
#   - processed_data_y::Vector{<:Number}
#
####################################################################################################

"""
    daEvolution(
        sim_data::Simulation,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute the time series of two quantities.

# Arguments

  - `sim_data::Simulation`: Information about the simulation in a [`Simulation`](@ref) object.
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
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `extra_filter::Function=filterNothing`: Filter function that will be applied after the one given by `filter_mode`.
  - `ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}()`: Request dictionary for the `extra_filter` filter function.
  - `smooth::Int=0`: The result of [`integrateQty`](@ref) will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.
  - `cumulative::Bool=false`: If the `y_quantity` will be accumulated or not.
  - `fraction::Bool=false`: If the `y_quantity` will be represented as a fraction of the last value. If `cumulative` = true, this will apply to the accumulated values.
  - `scaling::Function=identity`: Function to scale the x-axis (only relevant if `smooth` != 0). The bins will be computed accordingly. The options are the scaling functions accepted by Makie.jl: log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.

# Returns

  - A Tuple with two elements:

      + A Vector with the time series of `x_quantity`.
      + A Vector with the time series of `y_quantity`.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function daEvolution(
    sim_data::Simulation,
    x_quantity::Symbol,
    y_quantity::Symbol;
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    extra_filter::Function=filterNothing,
    ff_request::Dict{Symbol,Vector{String}}=Dict{Symbol,Vector{String}}(),
    smooth::Int=0,
    cumulative::Bool=false,
    fraction::Bool=false,
    scaling::Function=identity,
)::NTuple{2,Vector{<:Number}}

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        mergeRequests(plotParams(x_quantity).request, plotParams(y_quantity).request, ff_request),
    )

    # Iterate over each snapshot in the slice
    iterator = eachrow(DataFrame(sim_data.table[sim_data.slice, :]))

    # Allocate memory
    x_axis = Vector{Number}(fill(NaN, length(iterator)))
    y_axis = Vector{Number}(fill(NaN, length(iterator)))

    for (slice_index, sim_table_data) in pairs(iterator)

        global_index  = sim_table_data[1]
        scale_factor  = sim_table_data[3]
        redshift      = sim_table_data[4]
        physical_time = sim_table_data[5]
        lookback_time = sim_table_data[6]
        snapshot_path = sim_table_data[7]
        groupcat_path = sim_table_data[8]

        # Skip missing snapshots
        !ismissing(snapshot_path) || continue

        # Get the snapshot header
        snapshot_header = readSnapHeader(snapshot_path)

        # Get the group catalog header
        groupcat_header = readGroupCatHeader(groupcat_path)

        # Construct the metadata dictionary
        metadata = Dict(
            :sim_data => sim_data,
            :snap_data => Snapshot(
                snapshot_path,
                global_index,
                slice_index,
                physical_time,
                lookback_time,
                scale_factor,
                redshift,
                snapshot_header,
            ),
            :gc_data => GroupCatalog(groupcat_path, groupcat_header),
        )

        # Read the data in the snapshot
        data_dict = merge(
            metadata,
            readSnapshot(snapshot_path, request),
            readGroupCatalog(groupcat_path, snapshot_path, request),
        )

        # Filter the data
        filterData!(data_dict; filter_function)

        # Translate the data
        translateData!(data_dict, translation)

        # Rotate the data
        rotateData!(data_dict, rotation)

        # Filter the data again
        filterData!(data_dict; filter_function=extra_filter)

        # Compute the value for the x axis
        x_axis[slice_index] = integrateQty(data_dict, x_quantity)

        # Compute the value for the y axis
        y_axis[slice_index] = integrateQty(data_dict, y_quantity)

    end

    if cumulative
        cumsum!(y_axis, y_axis)
    end

    if fraction
        y_axis ./= y_axis[end]
    end

    if iszero(smooth)
        return x_axis, y_axis
    else
        return smoothWindow(x_axis, y_axis, smooth; scaling)
    end

end

"""
    daVirialAccretion(
        sim_data::Simulation;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute the evolution of the accreted mass into the virial radius.

# Arguments

  - `sim_data::Simulation`: Information about the simulation in a [`Simulation`](@ref) object.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted. Only valid if `tracers` = true. The options are:

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
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `halo_idx::Int=1`: Index of the target halo (FoF group). Starts at 1.
  - `tracers::Bool=false`: If tracers will be use to compute the mass accretion. If false, `filter_mode` will be ignored.
  - `smooth::Int=0`: The time series will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.

# Returns

  - A Tuple with two elements:

      + A Vector with the physical times.
      + A Vector with the accreted mass at each time.
"""
function daVirialAccretion(
    sim_data::Simulation;
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    halo_idx::Int=1,
    tracers::Bool=false,
    smooth::Int=0,
)::NTuple{2,Vector{<:Number}}

    filter_function, translation, _, request = selectFilter(
        filter_mode,
        Dict(
            :gas         => ["ID  ", "MASS"],
            :stars       => ["ID  ", "MASS"],
            :black_hole  => ["ID  ", "MASS"],
            :group       => ["G_R_Crit200", "G_M_Crit200"],
            :tracer      => ["PAID", "TRID"],
        ),
    )

    # Read the metadata table for the simulation
    simulation_dataframe = DataFrame(sim_data.table[sim_data.slice, :])

    # Delete missing snapshots
    filter!(row -> !ismissing(row[:snapshot_paths]), simulation_dataframe)

    # Iterate over each snapshot in the slice
    iterator = eachrow(simulation_dataframe)

    # Check that there are at least 2 snapshots left
    (
        length(iterator) >= 2 ||
        throw(ArgumentError("daVirialAccretion: The given slice, $(sim_data.slice), selected for \
        less than two snapshots. I need at least two snapshots to compute the a time series of \
        gas accretion. The full simulation table is:\n$(sim_data.table)"))
    )

    ################################################################################################
    # First element of the iteration over the snapshots
    ################################################################################################

    sim_table_data = iterator[1]

    snapshot_path = sim_table_data[7]
    groupcat_path = sim_table_data[8]

    # Get the snapshot header
    snapshot_header = readSnapHeader(snapshot_path)

    # Get the group catalog header
    groupcat_header = readGroupCatHeader(groupcat_path)

    # Construct the metadata dictionary
    metadata = Dict(
        :sim_data => sim_data,
        :snap_data => Snapshot(
            snapshot_path,
            sim_table_data[1],
            1,
            sim_table_data[5],
            sim_table_data[6],
            sim_table_data[3],
            sim_table_data[4],
            snapshot_header,
        ),
        :gc_data => GroupCatalog(groupcat_path, groupcat_header),
    )

    # Read the data in the snapshot
    past_dd = merge(
        metadata,
        readSnapshot(snapshot_path, request),
        readGroupCatalog(groupcat_path, snapshot_path, request),
    )

    if tracers
        # Filter the data
        filterData!(past_dd; filter_function)
        # Translate the data
        translateData!(past_dd, translation)
    end

    # Allocate memory fo the mass axis
    Δm = Vector{Unitful.Mass}(undef, length(iterator) - 1)

    ################################################################################################
    # Iteration over the snapshots
    ################################################################################################

    for (slice_index, sim_table_data) in pairs(iterator[2:end])

        global_index  = sim_table_data[1]
        scale_factor  = sim_table_data[3]
        redshift      = sim_table_data[4]
        physical_time = sim_table_data[5]
        lookback_time = sim_table_data[6]
        snapshot_path = sim_table_data[7]
        groupcat_path = sim_table_data[8]

        # Get the snapshot header
        snapshot_header = readSnapHeader(snapshot_path)

        # Get the group catalog header
        groupcat_header = readGroupCatHeader(groupcat_path)

        # Construct the metadata dictionary
        metadata = Dict(
            :sim_data => sim_data,
            :snap_data => Snapshot(
                snapshot_path,
                global_index,
                slice_index + 1,
                physical_time,
                lookback_time,
                scale_factor,
                redshift,
                snapshot_header,
            ),
            :gc_data => GroupCatalog(groupcat_path, groupcat_header),
        )

        # Read the data in the snapshot
        present_dd = merge(
            metadata,
            readSnapshot(snapshot_path, request),
            readGroupCatalog(groupcat_path, snapshot_path, request),
        )

        if tracers
            # Filter the data
            filterData!(present_dd; filter_function)
            # Translate the data
            translateData!(present_dd, translation)
        end

        if tracers

            Δm[slice_index], _, _ = computeVirialAccretion(present_dd, past_dd; halo_idx)

        else

            if isempty(past_dd[:group]["G_M_Crit200"])
                m_past = 0.0u"Msun"
            else
                m_past = past_dd[:group]["G_M_Crit200"][halo_idx]
            end

            if isempty(present_dd[:group]["G_M_Crit200"])
                m_present = 0.0u"Msun"
            else
                m_present = present_dd[:group]["G_M_Crit200"][halo_idx]
            end

            Δm[slice_index] = m_present -  m_past

        end

        past_dd = present_dd

    end

    # Compure the time ticks
    t  = sim_data.table[sim_data.slice, :physical_times]

    # Compure the time axis
    Δt = deltas(t)[2:end]

    if iszero(smooth)
        return t[2:end], Δm ./ Δt
    else
        return smoothWindow(t[2:end], Δm ./ Δt, smooth)
    end

end

"""
    daDiscAccretion(
        sim_data::Simulation;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute the evolution of the accreted mass into the disc.

# Arguments

  - `sim_data::Simulation`: Information about the simulation in a [`Simulation`](@ref) object.
  - `filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all`: Which cells/particles will be plotted. The options are:

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
              + `(halo_idx, subhalo_rel_idx)` -> Sets the principal axis of the stars in `subhalo_rel_idx::Int` subhalo (of the `halo_idx::Int` halo), as the new coordinate system.
              + `(halo_idx, 0)`               -> Sets the principal axis of the stars in the `halo_idx::Int` halo, as the new coordinate system.
              + `subhalo_abs_idx`             -> Sets the principal axis of the stars in the `subhalo_abs_idx::Int` subhalo as the new coordinate system.
  - `max_r::Unitful.Length=DISK_R`: Radius of the cylinder.
  - `max_z::Unitful.Length=5.0u"kpc"`: Half height of the cylinder.
  - `smooth::Int=0`: The time series will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.

# Returns

  - A Tuple with two elements:

      + A Vector with the physical times.
      + A Vector with the accreted mass at each time.
"""
function daDiscAccretion(
    sim_data::Simulation;
    filter_mode::Union{Symbol,Dict{Symbol,Any}}=:all,
    max_r::Unitful.Length=DISK_R,
    max_z::Unitful.Length=5.0u"kpc",
    smooth::Int=0,
)::NTuple{2,Vector{<:Number}}

    filter_function, translation, rotation, request = selectFilter(
        filter_mode,
        Dict(
            :gas         => ["ID  ", "MASS"],
            :stars       => ["ID  ", "MASS"],
            :black_hole  => ["ID  ", "MASS"],
            :tracer      => ["PAID", "TRID"],
        ),
    )

    # Read the metadata table for the simulation
    simulation_dataframe = DataFrame(sim_data.table[sim_data.slice, :])

    # Delete missing snapshots
    filter!(row -> !ismissing(row[:snapshot_paths]), simulation_dataframe)

    # Iterate over each snapshot in the slice
    iterator = eachrow(simulation_dataframe)

    # Check that there are at least 2 snapshots left
    (
        length(iterator) >= 2 ||
        throw(ArgumentError("daDiscAccretion: The given slice, $(sim_data.slice), selected for \
        less than two snapshots. I need at least two snapshots to compute the a time series of \
        gas accretion. The full simulation table is:\n$(sim_data.table)"))
    )

    ################################################################################################
    # First element of the iteration over the snapshots
    ################################################################################################

    sim_table_data = iterator[1]

    snapshot_path = sim_table_data[7]
    groupcat_path = sim_table_data[8]

    # Get the snapshot header
    snapshot_header = readSnapHeader(snapshot_path)

    # Get the group catalog header
    groupcat_header = readGroupCatHeader(groupcat_path)

    # Construct the metadata dictionary
    metadata = Dict(
        :sim_data => sim_data,
        :snap_data => Snapshot(
            snapshot_path,
            sim_table_data[1],
            1,
            sim_table_data[5],
            sim_table_data[6],
            sim_table_data[3],
            sim_table_data[4],
            snapshot_header,
        ),
        :gc_data => GroupCatalog(groupcat_path, groupcat_header),
    )

    # Read the data in the snapshot
    past_dd = merge(
        metadata,
        readSnapshot(snapshot_path, request),
        readGroupCatalog(groupcat_path, snapshot_path, request),
    )

    # Filter the data
    filterData!(past_dd; filter_function)

    # Translate the data
    translateData!(past_dd, translation)

    # Rotate the data
    rotateData!(past_dd, rotation)

    # Allocate memory fo the mass axis
    Δm = Vector{Unitful.Mass}(undef, length(iterator) - 1)

    ################################################################################################
    # Iteration over the snapshots
    ################################################################################################

    for (slice_index, sim_table_data) in pairs(iterator[2:end])

        global_index  = sim_table_data[1]
        scale_factor  = sim_table_data[3]
        redshift      = sim_table_data[4]
        physical_time = sim_table_data[5]
        lookback_time = sim_table_data[6]
        snapshot_path = sim_table_data[7]
        groupcat_path = sim_table_data[8]

        # Get the snapshot header
        snapshot_header = readSnapHeader(snapshot_path)

        # Get the group catalog header
        groupcat_header = readGroupCatHeader(groupcat_path)

        # Construct the metadata dictionary
        metadata = Dict(
            :sim_data => sim_data,
            :snap_data => Snapshot(
                snapshot_path,
                global_index,
                slice_index + 1,
                physical_time,
                lookback_time,
                scale_factor,
                redshift,
                snapshot_header,
            ),
            :gc_data => GroupCatalog(groupcat_path, groupcat_header),
        )

        # Read the data in the snapshot
        present_dd = merge(
            metadata,
            readSnapshot(snapshot_path, request),
            readGroupCatalog(groupcat_path, snapshot_path, request),
        )

        # Filter the data
        filterData!(present_dd; filter_function)

        # Translate the data
        translateData!(present_dd, translation)

        # Rotate the data
        rotateData!(present_dd, rotation)

        Δm[slice_index], _, _ = computeDiscAccretion(present_dd, past_dd; max_r, max_z)

        past_dd = present_dd

    end

    # Compure the time ticks
    t  = sim_data.table[sim_data.slice, :physical_times]

    # Compure the time axis
    Δt = deltas(t)[2:end]

    if iszero(smooth)
        return t[2:end], Δm ./ Δt
    else
        return smoothWindow(t[2:end], Δm ./ Δt, smooth)
    end

end

"""
    daSFRtxt(
        sim_data::Simulation,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute the stellar mass or SFR evolution using the data in the `sfr.txt` file.

# Arguments

  - `sim_data::Simulation`: Information about the simulation in a [`Simulation`](@ref) object.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:scale_factor`  -> Scale factor.
      + `:redshift`      -> Redshift.
      + `:physical_time` -> Physical time since the Big Bang.
      + `:lookback_time` -> Physical time left to reach the last snapshot.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass` -> Stellar mass.
      + `:sfr`          -> The star formation rate.
  - `smooth::Int=0`: The result will be smooth out using `smooth` bins. Set it to 0 if you want no smoothing.

# Returns

  - A Tuple with two elements:

      + A Vector with the time series of `x_quantity`.
      + A Vector with the time series of `y_quantity`.
"""
function daSFRtxt(
    sim_data::Simulation,
    x_quantity::Symbol,
    y_quantity::Symbol;
    smooth::Int=0,
)::NTuple{2,Vector{<:Number}}

    snapshot_paths = filter(!ismissing, sim_data.table[!, 7])

    (
        !isempty(snapshot_paths) ||
        throw(ArgumentError("daSFRtxt: I couldn't find any snapshots in $(sim_data.path), \
        and I need at least one for unit conversion"))
    )

    # Find the path to one snapshot
    snapshot_path = first(snapshot_paths)

    # Read its header
    header = readSnapHeader(snapshot_path)

    # Read the data in the `sfr.txt` file
    sfr_txt_data = readSfrFile(joinpath(sim_data.path, SFR_REL_PATH), snapshot_path)
    time_ticks = sfr_txt_data[1]

    if x_quantity == :scale_factor

        (
            sim_data.cosmological || !logging[] ||
            @warn("daSFRtxt: For non-cosmological simulations `x_quantity` can only be \
            :physical_time")
        )

        x_axis = time_ticks

    elseif x_quantity == :redshift

        if sim_data.cosmological

            x_axis = (1.0 ./ time_ticks) .- 1.0

        else

            (
                !logging[] ||
                @warn("daSFRtxt: For non-cosmological simulations `x_quantity` can only be \
                :physical_time")
            )

            x_axis = time_ticks

        end

    elseif x_quantity == :physical_time

        if sim_data.cosmological
            x_axis = computeTime(time_ticks, header)
        else
            x_axis = time_ticks
        end

    elseif x_quantity == :lookback_time

        if sim_data.cosmological
            physical_time = computeTime(time_ticks, header)
        else
            physical_time = time_ticks
        end

        x_axis = last(physical_time) .- physical_time

    else

        throw(ArgumentError("daSFRtxt: `x_quantity` can only be :scale_factor, :redshift, \
        :physical_time or :lookback_time, but I got :$(x_quantity)"))

    end

    if y_quantity == :stellar_mass

        y_axis = sfr_txt_data[6]

    elseif y_quantity == :sfr

        y_axis = sfr_txt_data[3]

    else

        throw(ArgumentError("daSFRtxt: `y_quantity` can only be :stellar_mass or :sfr, \
        but I got :$(y_quantity)"))

    end

    # Apply smoothing if required
    if !iszero(smooth)
        x_axis, y_axis = smoothWindow(x_axis, y_axis, smooth)
    end

    return x_axis, y_axis

end

"""
    daCPUtxt(
        sim_data::Simulation,
        target::String,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute the evolution of a measured quantity in the `cpu.txt` file, for a given `target` process.

# Arguments

  - `sim_data::Simulation`: Information about the simulation in a [`Simulation`](@ref) object.
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

# Returns

  - A Tuple with two elements:

      + A Vector with the time series of `x_quantity`.
      + A Vector with the time series of `y_quantity`.
"""
function daCPUtxt(
    sim_data::Simulation,
    target::String,
    x_quantity::Symbol,
    y_quantity::Symbol;
    smooth::Int=0,
)::NTuple{2,Vector{<:Number}}

    snapshot_paths = filter(!ismissing, sim_data.table[!, 7])

    (
        !isempty(snapshot_paths) ||
        throw(ArgumentError("daCPUtxt: I couldn't find any snapshots in $(sim_data.path), \
        and I need at least one for unit conversion"))
    )

    # Find the path to one snapshot
    snapshot_path = first(snapshot_paths)

    # Read its header
    header = readSnapHeader(snapshot_path)

    # Read the data in the `sfr.txt` file
    cpu_txt_data = readCpuFile(joinpath(sim_data.path, CPU_REL_PATH), [target])[target]

    if x_quantity == :time_step

        x_axis = cpu_txt_data[:, 1]

    elseif x_quantity == :physical_time

        if sim_data.cosmological
            x_axis = computeTime(cpu_txt_data[:, 2], header)
        else
            x_axis = cpu_txt_data[:, 2] .* internalUnits("CLKT", snapshot_path)
        end

    elseif x_quantity == :clock_time_s

        x_axis = cpu_txt_data[:, 3] .* u"s"

    elseif x_quantity == :clock_time_percent

        x_axis = cpu_txt_data[:, 4]

    elseif x_quantity == :tot_clock_time_s

        x_axis = cpu_txt_data[:, 5] .* u"s"

    elseif x_quantity == :tot_clock_time_percent

        x_axis = cpu_txt_data[:, 6]

    else

        throw(ArgumentError("daCPUtxt: I don't recognize the x_quantity :$(x_quantity)"))

    end

    if y_quantity == :time_step

        y_axis = cpu_txt_data[:, 1]

    elseif y_quantity == :physical_time

        if sim_data.cosmological
            y_axis = computeTime(cpu_txt_data[:, 2], header)
        else
            y_axis = cpu_txt_data[:, 2] .* internalUnits("CLKT", snapshot_path)
        end

    elseif y_quantity == :clock_time_s

        y_axis = cpu_txt_data[:, 3] .* u"s"

    elseif y_quantity == :clock_time_percent

        y_axis = cpu_txt_data[:, 4]

    elseif y_quantity == :tot_clock_time_s

        y_axis = cpu_txt_data[:, 5] .* u"s"

    elseif y_quantity == :tot_clock_time_percent

        y_axis = cpu_txt_data[:, 6]

    else

        throw(ArgumentError("daCPUtxt: I don't recognize the y_quantity :$(y_quantity)"))

    end

    # Apply smoothing if required
    if !iszero(smooth)
        x_axis, y_axis = smoothWindow(x_axis, y_axis, smooth)
    end

    return x_axis, y_axis

end

"""
    daVSFLaw(
        data_dict::Dict,
        grid::CubicGrid,
        quantity::Symbol;
        <keyword arguments>
    )::Union{NTuple{2,Vector{<:Float64}},Nothing}

Compute the gas mass density and the SFR density, used in the volumetric star formation (VSF) law.

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
  - `grid::CubicGrid`: Cubic grid.
  - `quantity::Symbol`: Quantity for the x axis. The options are:

      + `:gas_mass`          -> Gas density.
      + `:hydrogen_mass`     -> Hydrogen density.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) density.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) density.
      + `:ionized_mass`      -> Ionized hydrogen (``\\mathrm{HII}``) density.
      + `:neutral_mass`      -> Neutral hydrogen (``\\mathrm{HI + H_2}``) density.
  - `type::Symbol=:cells`: If the gas surface density will be calculated assuming the gas is in `:particles` or in Voronoi `:cells`.
  - `stellar_ff::Function=filterNothing`: Filter function for the stars. It has to be a function with the signature:

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
  - `gas_ff::Function=filterNothing`: Filter function for the gas. It has to be a function with the signature:

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

  - A tuple with two elements:

      + A vector with log10(ρH / M⊙ * pc^-2).
      + A vector with log10(ρsfr / M⊙ * yr^-1 * kpc^-2).

    It returns `nothing` if any of the necessary quantities are missing.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function daVSFLaw(
    data_dict::Dict,
    grid::CubicGrid,
    quantity::Symbol;
    type::Symbol=:cells,
    stellar_ff::Function=filterNothing,
    gas_ff::Function=filterNothing,
)::Union{NTuple{2,Vector{<:Float64}},Nothing}

    (
        quantity ∈ [:gas_mass, :molecular_mass, :br_molecular_mass, :neutral_mass] ||
        throw(ArgumentError("daVSFLaw: `quantity` can only be :gas_mass, \
        :molecular_mass, :br_molecular_mass or :neutral_mass, but I got :$(quantity)"))
    )

    # Log factor to go from stellar surface density to SFR surface density
    # log10(Σsfr) = log10(Σ*) - log10Δt
    log10Δt = log10(ustrip(u"yr", AGE_RESOLUTION))

    stellar_density = density3DProjection(
        data_dict,
        grid,
        :stellar_mass,
        :particles;
        m_unit=u"Msun",
        l_unit=u"kpc",
        filter_function=stellar_ff,
    )

    gas_density = density3DProjection(
        data_dict,
        grid,
        quantity,
        type;
        m_unit=u"Msun",
        l_unit=u"pc",
        filter_function=gas_ff,
    )

    x_axis = vec(gas_density)
    y_axis = vec(stellar_density)

    # Delete NaNs in the data vectors
    x_idxs = map(isnan, x_axis)
    y_idxs = map(isnan, y_axis)

    delete_idxs = x_idxs ∪ y_idxs

    deleteat!(x_axis, delete_idxs)
    deleteat!(y_axis, delete_idxs)

    !any(isempty, [x_axis, y_axis]) || return nothing

    return log10.(x_axis), log10.(y_axis) .- log10Δt

end

@doc raw"""
    daClumpingFactor(
        data_dict::Dict,
        quantity::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Volume},Vector{Float64}}

Compute the clumping factor (``C_\rho``), for the number density of `quantity`, at different volume scales.

```math
C_\rho = \frac{\langle n^2 \rangle}{\langle n \rangle^2} \, ,
```

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
  - `quantity::Symbol`: The number density of which quantity will be used. The options are:

      + `:gas`          -> Gas number density.
      + `:molecular`    -> Molecular hydrogen number density.
      + `:br_molecular` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic`       -> Atomic hydrogen number density.
      + `:ionized`      -> Ionized hydrogen number density.
      + `:neutral`      -> Neutral hydrogen number density.
  - `nn::Int=32`: Number of neighbors.
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

  - A tuple with two elements:

      + A vector with the volumes.
      + A vector with the clumping factors.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function daClumpingFactor(
    data_dict::Dict,
    quantity::Symbol;
    nn::Int=32,
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Volume},Vector{Float64}}

    filtered_dd = filterData(data_dict; filter_function)

    if quantity == :gas

        number_densities = scatterQty(filtered_dd, :gas_number_density)

    elseif quantity == :molecular

        number_densities = scatterQty(filtered_dd, :molecular_number_density)

    elseif quantity == :br_molecular

        number_densities = scatterQty(filtered_dd, :br_molecular_number_density)

    elseif quantity == :atomic

        number_densities = scatterQty(filtered_dd, :atomic_number_density)

    elseif quantity == :ionized

        number_densities = scatterQty(filtered_dd, :ionized_number_density)

    elseif quantity == :neutral

        number_densities = scatterQty(filtered_dd, :neutral_number_density)

    else

        throw(ArgumentError("daClumpingFactor: `quantity` can only be :gas, :molecular, \
        :br_molecular, :atomic, :ionized or :neutral, but I got :$(quantity)"))

    end

    # Load the position of each cell/particle
    positions = ustrip.(u"kpc", filtered_dd[:gas]["POS "])

    # Compute the volume of each cell/particle
    cell_volumes = filtered_dd[:gas]["MASS"] ./ filtered_dd[:gas]["RHO "]

    if any(isempty, [positions, cell_volumes, number_densities])
        return Float64[], Unitful.Volume[]
    end

    # Compute the tree for a nearest neighbor search
    kdtree = KDTree(positions)

    # Find the `nn` nearest cells/particles to each cell/particle
    idxs, _ = knn(kdtree, positions, nn, true)

    # Allocate memory
    Cρ = similar(number_densities, Float64)
    V  = similar(cell_volumes)

    (
        allequal(length, [V, Cρ, idxs]) ||
        throw(DomainError("daClumpingFactor: The lists of numer densities, volumes, and nearest \
        neighbor indices don't have the same lengths. This should not happen!"))
    )

    for (i, idx) in pairs(idxs)

        V[i]  = sum(cell_volumes[idx]; init=0.0u"kpc^3")
        Cρ[i] = computeClumpingFactor(number_densities[idx])

    end

    return V, Cρ

end

@doc raw"""
    daClumpingFactorProfile(
        data_dict::Dict,
        quantity::Symbol,
        grid::CircularGrid;
        <keyword arguments>
    )::Tuple{Vector{<:Unitful.Length},Vector{Float64}}

Compute a clumping factor (``C_\rho``) profile, for the number density of `quantity`.

```math
C_\rho = \frac{\langle n^2 \rangle}{\langle n \rangle^2} \, ,
```

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
  - `quantity::Symbol`: The number density of which quantity will be used. The options are:

      + `:gas`          -> Gas number density.
      + `:molecular`    -> Molecular hydrogen number density.
      + `:br_molecular` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic`       -> Atomic hydrogen number density.
      + `:ionized`      -> Ionized hydrogen number density.
      + `:neutral`      -> Neutral hydrogen number density.
  - `grid::CircularGrid`: Circular grid.
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

  - A tuple with two elements:

      + A vector with the central position of each bin.
      + A vector with the clumping factors.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function daClumpingFactorProfile(
    data_dict::Dict,
    quantity::Symbol,
    grid::CircularGrid;
    filter_function::Function=filterNothing,
)::Tuple{Vector{<:Unitful.Length},Vector{Float64}}

    filtered_dd = filterData(data_dict; filter_function)

    if quantity == :gas

        number_densities = scatterQty(filtered_dd, :gas_number_density)

    elseif quantity == :molecular

        number_densities = scatterQty(filtered_dd, :molecular_number_density)

    elseif quantity == :br_molecular

        number_densities = scatterQty(filtered_dd, :br_molecular_number_density)

    elseif quantity == :atomic

        number_densities = scatterQty(filtered_dd, :atomic_number_density)

    elseif quantity == :ionized

        number_densities = scatterQty(filtered_dd, :ionized_number_density)

    elseif quantity == :neutral

        number_densities = scatterQty(filtered_dd, :neutral_number_density)

    else

        throw(ArgumentError("daClumpingFactor: `quantity` can only be :gas, :molecular, \
        :br_molecular, :atomic, :ionized or :neutral, but I got :$(quantity)"))

    end

    # Load the position of each cell/particle
    positions = filtered_dd[:gas]["POS "]

    # Compute the radial distance of each cell/particle
    distances = computeDistance(positions[1:2, :]; center=grid.center[1:2])

    # Find which cells/particles fall within each bin of `grid`
    n_profile = listHistogram1D(distances, number_densities, grid)

    # Allocate memory
    Cρ = similar(grid.grid, Float64)

    for (i, n) in pairs(n_profile)

        Cρ[i] = computeClumpingFactor(n)

    end

    return grid.grid, Cρ

end
