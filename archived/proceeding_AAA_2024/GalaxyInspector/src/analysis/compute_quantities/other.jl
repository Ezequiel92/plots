####################################################################################################
# Computation of derived quantities
####################################################################################################

@doc raw"""
    computeEqQuotient(data_dict::Dict, type::Symbol)::Vector{Float64}

Compute the equilibrium quotient for the molecular or ionized equations of the SF model.

Molecular equation

From

```math
	0 = \frac{f_a}{\tau_\mathrm{cond}} - \eta_\mathrm{diss} \, \psi - \psi \, ,
```
and using

```math
    \tau_\mathrm{cond} = \frac{C_\mathrm{cond}}{(Z + Z_\mathrm{eff}) \, \rho_\mathrm{cell} \, (1 - f_s)} \, ,
```
and
```math
	\psi = \frac{f_m}{\tau_\mathrm{star}} \, .
```

We get

```math
	(\eta_\mathrm{diss} + 1) \, \frac{f_m^0}{\tau_\mathrm{star}} = \frac{f_a^0}{C_\mathrm{cond}} \, (Z + Z_\mathrm{eff}) \, \rho_\mathrm{cell} \, (1 - f_s^0) \, .
```

Which can be rewritten as

```math
	\frac{f_a^0}{f_m^0} \, (1 - f_s^0) = \frac{\eta_\mathrm{diss} + 1}{\tau_\mathrm{star}} \, \frac{C_\mathrm{cond}}{(Z + Z_\mathrm{eff}) \, \rho_\mathrm{cell}} \, ,
```
where we use the notation $f_X^0$ to indicate equilibrium fractions.

This last expression shows the value of $f_a^0$, $f_m^0$, and $f_s^0$ (as a function of the parameters of the model) that make the molecular equation have a $0$ derivative (equilibrium point for that particular equation).

Ionized equation

From

```math
	0 = - \frac{f_i}{\tau_\mathrm{rec}} + \eta_\mathrm{ion} \, \psi + R \, \psi \, ,
```
and using

```math
    \tau_\mathrm{rec} = \frac{C_\mathrm{rec}}{f_i \, \rho_\mathrm{cell}} \, ,
```
and
```math
	\psi = \frac{f_m}{\tau_\mathrm{star}} \, .
```

We get

```math
	0 = - \frac{f_i^2 \, \rho_\mathrm{cell}}{C_\mathrm{rec}} + (\eta_\mathrm{ion} + R) \, \frac{f_m}{\tau_\mathrm{star}} \, .
```

Which can be rewritten as

```math
	\frac{(f_i^0)^2}{f_m^0} = \frac{\eta_\mathrm{ion} + R}{\tau_\mathrm{star}} \, \frac{C_\mathrm{rec}}{\rho_\mathrm{cell}} \, .
```

As before, this last expression shows the value of $f_i^0$ and $f_m^0$ (as a function of the parameters of the model) that make the ionized equation have a $0$ derivative (equilibrium point for that particular equation).

# Arguments

  - `data_dict::Dict`: Scale factors.
  - `type::Symbol`: If the :molecular or the :ionized equation will be used.

# Returns

  - A vector with the equilibrium quotient for the molecular or ionized equations.
"""
function computeEqQuotient(data_dict::Dict, type::Symbol)::Vector{Float64}

    dg = data_dict[:gas]

    !any(isempty, [dg["FRAC"], dg["RHO "]]) || return Float64[]

    # Allocate memory
    eq_quotients = fill(NaN, length(dg["RHO "]))

    if type == :molecular

        !any(isempty, [dg["ETAD"], dg["PARZ"]]) || return Float64[]

        iterator = zip(
            dg["ETAD"],
            dg["FRAC"][2, :],
            dg["FRAC"][3, :],
            dg["FRAC"][4, :],
            τ_star.(dg["RHO "]),
            τ_cond.(dg["RHO "], dg["PARZ"]),
        )

        for (i, (ηd, fa, fm, fs, τS, τC)) in enumerate(iterator)

            !(isnan(fa) || iszero(fa) || isone(fs) || iszero(fm)) || continue

            mol_ls = (fa / fm) * (1 - fs)
            mol_rs = uconvert(Unitful.NoUnits, ((ηd + 1) * τC) / τS)

            eq_quotients[i] = log10(mol_ls / mol_rs)

        end

    elseif type == :ionized

        !any(isempty, [dg["ETAI"], dg["PARR"]]) || return Float64[]

        iterator = zip(
            dg["ETAI"],
            dg["PARR"],
            dg["FRAC"][1, :],
            dg["FRAC"][3, :],
            GalaxyInspector.τ_star.(dg["RHO "]),
            GalaxyInspector.τ_rec.(dg["RHO "]),
        )

        for (i, (ηi, R, fi, fm, τS, τR)) in enumerate(iterator)

            !(isnan(fi) || iszero(fi) || iszero(fm)) || continue

            ion_ls = (fi * fi) / fm
            ion_rs = uconvert(Unitful.NoUnits, ((ηi + R) * τR) / τS)

            eq_quotients[i] = log10(ion_ls / ion_rs)

        end

    else

        throw(ArgumentError("computeEqQuotient: `type` can only be :molecular or :ionized, \
        but I got :$(type)"))

    end

    return eq_quotients

end

@doc raw"""
    computeTime(
        scale_factors::Vector{<:Real},
        header::SnapshotHeader;
        <keyword arguments>
    )::Vector{<:Unitful.Time}

Compute the physical time corresponding to each of the `scale_factors`.

To get the physical time $t$ from the scale factor `a`, one does the integral:

```math
t = \frac{1}{H_0} \int_0^a \frac{\mathrm{d}a'}{a' \, \sqrt{\mathcal{E}(a')}} \, ,
```

where

```math
\mathcal{E}(a) = \Omega_\Lambda + \Omega_m \, a^{-3} + \Omega_r \, a^{-4} + \Omega_K \, a^{-2} \, .
```

# Arguments

  - `scale_factors::Vector{<:Real}`: Scale factors.
  - `header::SnapshotHeader`: A header of the simulation, containing the cosmological parameters.
  - `a0::Float64=0.0`: Initial scale factor.

# Returns

  - A vector with the physical times.
"""
function computeTime(
    scale_factors::Vector{<:Real},
    header::SnapshotHeader;
    a0::Float64=0.0,
)::Vector{<:Unitful.Time}

    f = x -> energyIntegrand(x, header)

    return [quadgk(f, a0, a)[1] * u"Gyr" for a in scale_factors]

end

@doc raw"""
    computeTime(a::Real, header::SnapshotHeader; <keyword arguments>)::Unitful.Time

Compute the physical time corresponding to the scale factor `a`.

To get the physical time $t$ from the scale factor `a`, one does the integral:

```math
t = \frac{1}{H_0} \int_0^a \frac{\mathrm{d}a'}{a' \, \sqrt{\mathcal{E}(a')}} \, ,
```

where

```math
\mathcal{E}(a) = \Omega_\Lambda + \Omega_m \, a^{-3} + \Omega_r \, a^{-4} + \Omega_K \, a^{-2} \, .
```

# Arguments

  - `a::Real`: Scale factor.
  - `header::SnapshotHeader`: A header of the simulation, containing the cosmological parameters.
  - `a0::Float64=0.0`: Initial scale factor.

# Returns

  - The physical time.
"""
function computeTime(a::Real, header::SnapshotHeader; a0::Float64=0.0)::Unitful.Time

    return computeTime([a], header; a0)[1]

end

"""
    computeTimeTicks(
        paths::Vector{<:Union{Missing,String}},
    )::Tuple{Vector{Float64},Vector{Float64},Vector{<:Unitful.Time},Vector{<:Unitful.Time}}

Compute the different times stamps associated with each snapshot in `paths`.

# Arguments

  - `paths::Vector{<:Union{Missing,String}}`: Paths to the snapshots.

# Returns

  - A tuple with four elements:

      + A vector with the scale factors.
      + A vector with the redshifts.
      + A vector with the physical times (physical time since the Big Bang).
      + A vector with the lookback times (physical time left to reach the last snapshot).
"""
function computeTimeTicks(
    paths::Vector{<:Union{Missing,String}},
)::Tuple{Vector{Float64},Vector{Float64},Vector{<:Unitful.Time},Vector{<:Unitful.Time}}

    snapshot_paths = filter(!ismissing, paths)

    !isempty(snapshot_paths) || return [NaN], [NaN], [NaN*u"s"], [NaN*u"s"]

    first_snapshot = first(snapshot_paths)

    if isCosmological(first_snapshot)

        # For cosmological simulations, the time field in the Header of the snapshot is the scale factor
        scale_factors = [readTime(path) for path in snapshot_paths]
        redshifts = @. (1.0 / scale_factors) - 1.0
        physical_times = computeTime(scale_factors, readSnapHeader(first_snapshot))
        lookback_times = last(physical_times) .- physical_times

    else

        # Compute the factor for internal units of time
        u_time = internalUnits("CLKT", first_snapshot)

        # a = 1.0 for non-cosmological simulations
        scale_factors = ones(length(snapshot_paths))
        # z = 0.0 for non-cosmological simulations
        redshifts = zeros(length(snapshot_paths))
        # For non-cosmological simulations, the time in the snapshot is the physical time
        physical_times = [readTime(path) * u_time for path in snapshot_paths]
        lookback_times = last(physical_times) .- physical_times

    end

    return scale_factors, redshifts, physical_times, lookback_times

end

"""
    computeTemperature(
        internal_energy::Vector{<:SpecificEnergy},
        electron_fraction::Vector{Float32},
    )::Vector{<:Unitful.Temperature}

Compute the gas temperature.

# Arguments

  - `internal_energy::Vector{<:SpecificEnergy}`: Specific internal energy of every gas cell/particle.
  - `electron_fraction::Vector{Float32}`: Number fraction of electrons in every gas cell/particle.

# Returns

  - The temperature of each gas cell/particle.
"""
function computeTemperature(
    internal_energy::Vector{<:SpecificEnergy},
    electron_fraction::Vector{Float32},
)::Vector{<:Unitful.Temperature}

    # xH := mass_fraction_of_hydrogen
    xH = HYDROGEN_MASSFRAC

    # yHe := number_of_helium_atoms / number_of_hydrogen_atoms
    # Take the mass fraction of metals as negligible
    yHe = @. (1.0 - xH) / (4.0 * xH)

    # electron_fraction := number_of_electrons / number_of_hydrogen_atoms
    # μ := total_mass / (total_number_of_particles * proton_mass)
    #   ≈ number_of_protons / total_number_of_particles
    # For the total mass, take the mass of electrons as negligible
    μ = @. (1.0 + 4.0 * yHe) / (1.0 + yHe + electron_fraction)

    # T = (adiabatic_index - 1) * internal_energy_per_unit_mass *
    #     (total_mass / total_number_of_particles) / boltzmann_constant
    return @. 0.6667 * internal_energy * μ * Unitful.mp / Unitful.k

end

"""
    computeStellarAge(data_dict::Dict)::Vector{<:Unitful.Time}

Compute the age of the stars.

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

  - The stellar ages.
"""
function computeStellarAge(data_dict::Dict)::Vector{<:Unitful.Time}

    birth_ticks = data_dict[:stars]["GAGE"]

    !isempty(birth_ticks) || return Unitful.Time[]

    if data_dict[:sim_data].cosmological
        # Go from scale factor to physical time
        birth_times = computeTime(birth_ticks, data_dict[:snap_data].header)
    else
        birth_times = birth_ticks
    end

    return data_dict[:snap_data].physical_time .- birth_times

end

"""
    computeSFR(
        data_dict::Dict;
        <keyword arguments>
    )::Vector{<:Unitful.MassFlow}

Compute the star formation rate of each stellar particle.

For stellar particles younger that `age_resol`, the SFR is its mass divided by `age_resol`. It is defined as 0 for older particles.

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
  - `age_resol::Unitful.Time=AGE_RESOLUTION`: Age resolution for the SFR.

# Returns

  - The star formation rate of each stellar particle.
"""
function computeSFR(
    data_dict::Dict;
    age_resol::Unitful.Time=AGE_RESOLUTION,
)::Vector{<:Unitful.MassFlow}

    # Compute the stellar ages
    ages = computeStellarAge(data_dict)

    !isempty(ages) || return Unitful.MassFlow[]

    # Allocate memory
    sfr = zeros(typeof(1.0u"Msun*yr^-1"), length(ages))

    # Find the stellar particles younger than `age_resol`
    idxs = map(x -> x <= age_resol, ages)

    # Compute the SFR
    sfr[idxs] .= data_dict[:stars]["MASS"][idxs] ./ age_resol

    return sfr

end

@doc raw"""
    computeClumpingFactor(density::Vector{<:Number})::Float64

Compute the clumping factor,

```math
C_\rho = \frac{\langle \rho^2 \rangle}{\langle \rho \rangle^2} \, .
```

# Arguments

  - `density::Vector{<:Number}`: The density of the cells/particles.

# Returns

  - The clumping factor.
"""
function computeClumpingFactor(density::Vector{<:Number})::Float64

    !isempty(density) || return NaN

    μ, var = mean_and_var(density)

    return 1.0 + uconvert(Unitful.NoUnits, var / μ^2)

end

@doc raw"""
    computeDepletionTime(
        mass::Vector{<:Unitful.Mass},
        sfr::Vector{<:Unitful.MassFlow},
    )::Vector{<:Unitful.Time}

Compute the depletion time,

```math
t_\mathrm{ff} = \frac{M_\mathrm{gas}}{\dot{M}_\star} \, .
```

# Arguments

  - `mass::Vector{<:Unitful.Mass}`: The gas mass of the cells/particles.
  - `sfr::Vector{<:Unitful.MassFlow}`: The SFR associated to each cell/particle.

# Returns

  - The depletion time.
"""
function computeDepletionTime(
    mass::Vector{<:Unitful.Mass},
    sfr::Vector{<:Unitful.MassFlow},
)::Vector{<:Unitful.Time}

    if any(isempty, [mass, sfr])

        (
            !logging[] ||
            @warn("computeDepletionTime: There is missing data for the gas, so the result will be \
            an empty vector")
        )

        return Unitful.Time[]

    end

    return @. mass / sfr

end

@doc raw"""
    computeEfficiencyFF(
        density::Vector{<:Unitful.Density},
        mass::Vector{<:Unitful.Mass},
        sfr::Vector{<:Unitful.MassFlow},
    )::Vector{Float64}

Compute the star formation efficiency per free-fall time, according to the definition in equation 1 of Krumholz (2012),

```math
\epsilon_\mathrm{ff} = \frac{t_\mathrm{ff}}{t_\mathrm{dep}} \, .
```
where

```math
t_\mathrm{ff} = \sqrt{\frac{3 \, \pi}{32 \, G \, \rho}} \, ,
```
is the free-fall time, and

```math
t_\mathrm{ff} = \frac{M_\mathrm{H_2}}{\dot{M}_\star} \, ,
```
is the depletion time.

# Arguments

  - `density::Vector{<:Unitful.Density}`: The molecular hydrogen (``\\mathrm{H_2}``) density of the cells/particles.
  - `mass::Vector{<:Unitful.Mass}`: The gas mass of the cells/particles.
  - `sfr::Vector{<:Unitful.MassFlow}`: The SFR associated to each cell/particle.

# Returns

  - The star formation efficiency per free-fall time.

# References

M. R. Krumholz et al. (2011). *A UNIVERSAL, LOCAL STAR FORMATION LAW IN GALACTIC CLOUDS, NEARBY GALAXIES, HIGH-REDSHIFT DISKS, AND STARBURSTS*. The Astrophysical Journal, **745(1)**, 69. [doi:10.1088/0004-637X/745/1/69](https://doi.org/10.1088/0004-637X/745/1/69)
"""
function computeEfficiencyFF(
    density::Vector{<:Unitful.Density},
    mass::Vector{<:Unitful.Mass},
    sfr::Vector{<:Unitful.MassFlow},
)::Vector{Float64}

    if any(isempty, [density, mass, sfr])

        (
            !logging[] ||
            @warn("computeEfficiencyFF: There is missing data for the gas, so the result will be \
            an empty vector")
        )

        return Float64[]

    end

    ϵff = Vector{Float64}(undef, length(density))

    for i in eachindex(ϵff)

        if density[i] < THRESHOLD_DENSITY

            ϵff[i] = 0.0

        else

            # Compute the free-fall time
            tff = sqrt(3π / (32 * Unitful.G * density[i]))

            # Compute the depletion time
            tdep = mass[i] / sfr[i]

            ϵff[i] = uconvert(Unitful.NoUnits, tff / tdep)

        end

    end

    return ϵff

end

"""
    integrateQty(data_dict::Dict, quantity::Symbol)::Number

Compute an integrated quantity for the whole system in `data_dict`.

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

# Returns

  - The velue of `quantity` for the whole system in `data_dict`.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function integrateQty(data_dict::Dict, quantity::Symbol)::Number

    if quantity == :stellar_mass

        integrated_qty = sum(computeMass(data_dict, :stars); init=0.0u"Msun")

    elseif quantity == :gas_mass

        integrated_qty = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

    elseif quantity == :hydrogen_mass

        integrated_qty = sum(computeMass(data_dict, :hydrogen); init=0.0u"Msun")

    elseif quantity == :dm_mass

        integrated_qty = sum(computeMass(data_dict, :dark_matter); init=0.0u"Msun")

    elseif quantity == :bh_mass

        integrated_qty = sum(computeMass(data_dict, :black_holes); init=0.0u"Msun")

    elseif quantity == :molecular_mass

        integrated_qty = sum(computeMass(data_dict, :molecular); init=0.0u"Msun")

    elseif quantity == :br_molecular_mass

        integrated_qty = sum(computeMass(data_dict, :br_molecular); init=0.0u"Msun")

    elseif quantity == :atomic_mass

        integrated_qty = sum(computeMass(data_dict, :atomic); init=0.0u"Msun")

    elseif quantity == :ionized_mass

        integrated_qty = sum(computeMass(data_dict, :ionized); init=0.0u"Msun")

    elseif quantity == :neutral_mass

        integrated_qty = sum(computeMass(data_dict, :neutral); init=0.0u"Msun")

    elseif quantity == :stellar_gas_mass

        integrated_qty = sum(computeMass(data_dict, :stellar); init=0.0u"Msun")

    elseif quantity == :stellar_number

        integrated_qty = length(data_dict[:stars]["MASS"])

    elseif quantity == :gas_number

        integrated_qty = length(data_dict[:gas]["MASS"])

    elseif quantity == :dm_number

        integrated_qty = length(data_dict[:halo]["MASS"])

    elseif quantity == :bh_number

        integrated_qty = length(data_dict[:black_hole]["MASS"])

    elseif quantity == :molecular_fraction

        molecular_mass = sum(computeMass(data_dict, :molecular); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = molecular_mass / gas_mass
        end

    elseif quantity == :br_molecular_fraction

        molecular_mass = sum(computeMass(data_dict, :br_molecular); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = molecular_mass / gas_mass
        end

    elseif quantity == :atomic_fraction

        atomic_mass = sum(computeMass(data_dict, :atomic); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = atomic_mass / gas_mass
        end

    elseif quantity == :ionized_fraction

        ionized_mass = sum(computeMass(data_dict, :ionized); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = ionized_mass / gas_mass
        end

    elseif quantity == :neutral_fraction

        neutral_mass = sum(computeMass(data_dict, :neutral); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = neutral_mass / gas_mass
        end

    elseif quantity == :molecular_neutral_fraction

        molecular_mass = sum(computeMass(data_dict, :molecular); init=0.0u"Msun")
        neutral_mass = sum(computeMass(data_dict, :neutral); init=0.0u"Msun")

        if iszero(neutral_mass)
            integrated_qty = NaN
        else
            integrated_qty = molecular_mass / neutral_mass
        end

    elseif quantity == :ionized_neutral_fraction

        ionized_mass = sum(computeMass(data_dict, :ionized); init=0.0u"Msun")
        neutral_mass = sum(computeMass(data_dict, :atomic); init=0.0u"Msun")

        if iszero(neutral_mass)
            integrated_qty = NaN
        else
            integrated_qty = ionized_mass / neutral_mass
        end

    elseif quantity == :gas_mass_density

        density = computeVolumeDensity(data_dict, :gas)

        if isempty(density)
            integrated_qty = 0.0u"Msun * kpc^-3"
        else
            integrated_qty = mean(density)
        end

    elseif quantity == :stellar_gas_fraction

        stellar_gas_mass = sum(computeMass(data_dict, :stellar); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = neutral_mass / stellar_gas_mass
        end

    elseif quantity == :stellar_area_density

        integrated_qty = sum(computeMass(data_dict, :stars); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :gas_area_density

        integrated_qty = sum(computeMass(data_dict, :gas); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :molecular_area_density

        integrated_qty = sum(computeMass(data_dict, :molecular); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :br_molecular_area_density

        integrated_qty = sum(computeMass(data_dict, :br_molecular); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :atomic_area_density

        integrated_qty = sum(computeMass(data_dict, :atomic); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :ionized_area_density

        integrated_qty = sum(computeMass(data_dict, :ionized); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :neutral_area_density

        integrated_qty = sum(computeMass(data_dict, :neutral); init=0.0u"Msun") / area(DISK_R)

    elseif quantity == :sfr_area_density

        sfr = sum(computeSFR(data_dict; age_resol=AGE_RESOLUTION); init=0.0u"Msun*yr^-1")

        integrated_qty = sfr / area(DISK_R)

    elseif quantity == :gas_td

        mass = computeMass(data_dict, :gas)
        td   = computeDepletionTime(mass, data_dict[:gas]["SFR "])

        if isempty(td)
            integrated_qty = NaN
        else
            integrated_qty = mean(td)
        end

    elseif quantity == :molecular_td

        mass = computeMass(data_dict, :molecular)
        td   = computeDepletionTime(mass, data_dict[:gas]["SFR "])

        if isempty(td)
            integrated_qty = NaN
        else
            integrated_qty = mean(td)
        end

    elseif quantity == :br_molecular_td

        mass = computeMass(data_dict, :br_molecular)
        td   = computeDepletionTime(mass, data_dict[:gas]["SFR "])

        if isempty(td)
            integrated_qty = NaN
        else
            integrated_qty = mean(td)
        end

    elseif quantity == :atomic_td

        mass = computeMass(data_dict, :atomic)
        td   = computeDepletionTime(mass, data_dict[:gas]["SFR "])

        if isempty(td)
            integrated_qty = NaN
        else
            integrated_qty = mean(td)
        end

    elseif quantity == :ionized_td

        mass = computeMass(data_dict, :ionized)
        td   = computeDepletionTime(mass, data_dict[:gas]["SFR "])

        if isempty(td)
            integrated_qty = NaN
        else
            integrated_qty = mean(td)
        end

    elseif quantity == :neutral_td

        mass = computeMass(data_dict, :neutral)
        td   = computeDepletionTime(mass, data_dict[:gas]["SFR "])

        if isempty(td)
            integrated_qty = NaN
        else
            integrated_qty = mean(td)
        end

    elseif quantity == :gas_metallicity

        metal_mass = sum(computeMetalMass(data_dict, :gas); init=0.0u"Msun")
        gas_mass = sum(computeMass(data_dict, :gas); init=0.0u"Msun")

        if iszero(gas_mass)
            integrated_qty = NaN
        else
            integrated_qty = (metal_mass / gas_mass) / SOLAR_METALLICITY
        end

    elseif quantity == :stellar_metallicity

        metal_mass = sum(computeMetalMass(data_dict, :stars); init=0.0u"Msun")
        stellar_mass = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun")

        if iszero(stellar_mass)
            integrated_qty = NaN
        else
            integrated_qty = (metal_mass / stellar_mass) / SOLAR_METALLICITY
        end

    elseif quantity ∈ GAS_ABUNDANCE

        element_symbol = Symbol(first(split(string(quantity), "_")))

        abundance = 12 + log10(computeGlobalAbundance(data_dict, :gas, element_symbol))
        integrated_qty = isinf(abundance) ? NaN : abundance

    elseif quantity ∈ STELLAR_ABUNDANCE

        element_symbol = Symbol(first(split(string(quantity), "_")))

        abundance = 12 + log10(computeGlobalAbundance(data_dict, :stars, element_symbol))
        integrated_qty = isinf(abundance) ? NaN : abundance

    elseif quantity == :stellar_specific_am

        positions = data_dict[:stars]["POS "]
        velocities = data_dict[:stars]["VEL "]
        masses = data_dict[:stars]["MASS"]

        if any(isempty, [positions, velocities, masses])
            integrated_qty = NaN
        else
            J = norm(computeTotalAngularMomentum(positions, velocities, masses; normal=false))
            integrated_qty = J / sum(masses)
        end

    elseif quantity == :gas_specific_am

        positions = data_dict[:gas]["POS "]
        velocities = data_dict[:gas]["VEL "]
        masses = computeMass(data_dict, :gas)

        if any(isempty, [positions, velocities, masses])
            integrated_qty = NaN
        else
            J = norm(computeTotalAngularMomentum(positions, velocities, masses; normal=false))
            integrated_qty = J / sum(masses)
        end

    elseif quantity == :dm_specific_am

        positions = data_dict[:halo]["POS "]
        velocities = data_dict[:halo]["VEL "]
        masses = data_dict
        masses = [:halo]["MASS"]

        if any(isempty, [positions, velocities, masses])
            integrated_qty = NaN
        else
            J = norm(computeTotalAngularMomentum(positions, velocities, masses; normal=false))
            integrated_qty = J / sum(masses)
        end

    elseif quantity == :sfr

        # Get the global index (index in the context of the whole simulation) of the current snapshot
        present_idx = data_dict[:snap_data].global_index

        if present_idx == 1

            integrated_qty = 0.0u"Msun*yr^-1"

        else

            # Get the physical times
            times = data_dict[:sim_data].table[:, 5]
            # Compute the time between snapshots
            Δt = times[present_idx] - times[present_idx - 1]

            integrated_qty = sum(computeSFR(data_dict; age_resol=Δt); init=0.0u"Msun*yr^-1")

        end

    elseif quantity == :ssfr

        # Get the global index (index in the context of the whole simulation) of the current snapshot
        present_idx = data_dict[:snap_data].global_index

        # Compute the total stellar mass
        stellar_mass = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun")

        if present_idx == 1 || iszero(stellar_mass)

            integrated_qty = 0.0u"yr^-1"

        else

            # Get the physical times
            times = data_dict[:sim_data].table[:, 5]
            # Compute the time between snapshots
            Δt = times[present_idx] - times[present_idx - 1]

            integrated_qty = sum(
                computeSFR(data_dict; age_resol=Δt);
                init=0.0u"Msun*yr^-1",
            ) / stellar_mass

        end

    elseif quantity == :observational_sfr

        integrated_qty = sum(computeSFR(data_dict; age_resol=AGE_RESOLUTION); init=0.0u"Msun*yr^-1")

    elseif quantity == :observational_ssfr

        sfr = sum(computeSFR(data_dict; age_resol=AGE_RESOLUTION); init=0.0u"Msun*yr^-1")
        stellar_mass = sum(data_dict[:stars]["MASS"]; init=0.0u"Msun")

        if iszero(stellar_mass)
            integrated_qty = 0.0u"yr^-1"
        else
            integrated_qty = sfr / stellar_mass
        end

    elseif quantity == :stellar_eff

        ϵffs = computeEfficiencyFF(
            data_dict[:stars]["RHOC"] .* u"mp",
            data_dict[:stars]["GMAS"],
            data_dict[:stars]["GSFR"],
        )

        # Filter out zeros and NaNs
        filter!(eff -> !isnan(eff) && !iszero(eff), ϵffs)

        if isempty(ϵffs)
            integrated_qty = NaN
        else
            integrated_qty = mean(ϵffs)
        end

    elseif quantity == :gas_eff

        mass = computeMass(data_dict, :gas)

        ϵffs = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

        # Filter out zeros and NaNs
        filter!(eff -> !isnan(eff) && !iszero(eff), ϵffs)

        if isempty(ϵffs)
            integrated_qty = NaN
        else
            integrated_qty = mean(ϵffs)
        end

    elseif quantity == :molecular_eff

        mass = computeMass(data_dict, :molecular)

        ϵffs = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

        # Filter out zeros and NaNs
        filter!(eff -> !isnan(eff) && !iszero(eff), ϵffs)

        if isempty(ϵffs)
            integrated_qty = NaN
        else
            integrated_qty = mean(ϵffs)
        end

    elseif quantity == :br_molecular_eff

        mass = computeMass(data_dict, :br_molecular)

        ϵffs = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

        # Filter out zeros and NaNs
        filter!(eff -> !isnan(eff) && !iszero(eff), ϵffs)

        if isempty(ϵffs)
            integrated_qty = NaN
        else
            integrated_qty = mean(ϵffs)
        end

    elseif quantity == :atomic_eff

        mass = computeMass(data_dict, :atomic)

        ϵffs = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

        # Filter out zeros and NaNs
        filter!(eff -> !isnan(eff) && !iszero(eff), ϵffs)

        if isempty(ϵffs)
            integrated_qty = NaN
        else
            integrated_qty = mean(ϵffs)
        end

    elseif quantity == :ionized_eff

        mass = computeMass(data_dict, :ionized)

        ϵffs = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

        # Filter out zeros and NaNs
        filter!(eff -> !isnan(eff) && !iszero(eff), ϵffs)

        if isempty(ϵffs)
            integrated_qty = NaN
        else
            integrated_qty = mean(ϵffs)
        end

    elseif quantity == :neutral_eff

        mass = computeMass(data_dict, :neutral)

        ϵffs = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

        # Filter out zeros and NaNs
        filter!(eff -> !isnan(eff) && !iszero(eff), ϵffs)

        if isempty(ϵffs)
            integrated_qty = NaN
        else
            integrated_qty = mean(ϵffs)
        end

    elseif quantity == :scale_factor

        integrated_qty = data_dict[:sim_data].table[data_dict[:snap_data].global_index, 3]

    elseif quantity == :redshift

        integrated_qty = data_dict[:sim_data].table[data_dict[:snap_data].global_index, 4]

    elseif quantity == :physical_time

        integrated_qty = data_dict[:sim_data].table[data_dict[:snap_data].global_index, 5]

    elseif quantity == :lookback_time

        integrated_qty = data_dict[:sim_data].table[data_dict[:snap_data].global_index, 6]

    elseif quantity == :ode_gas_it

        odit = data_dict[:gas]["ODIT"]
        filter!(!isnan, odit)

        if isempty(odit)
            integrated_qty = NaN
        else
            integrated_qty = mean(odit)
        end

    elseif quantity == :ode_gas_accu_it

        acit = data_dict[:gas]["ACIT"]
        filter!(!isnan, acit)

        if isempty(acit)
            integrated_qty = NaN
        else
            integrated_qty = mean(acit)
        end

    elseif quantity == :ode_gas_tau_s

        τS = data_dict[:gas]["TAUS"]
        filter!(!isnan, τS)

        if isempty(τS)
            integrated_qty = NaN
        else
            integrated_qty = mean(τS)
        end

    elseif quantity == :ode_gas_eta_d

        ηd = data_dict[:gas]["ETAD"]
        filter!(!isnan, ηd)

        if isempty(ηd)
            integrated_qty = NaN
        else
            integrated_qty = mean(ηd)
        end

    elseif quantity == :ode_gas_eta_i

        ηi = data_dict[:gas]["ETAI"]
        filter!(!isnan, ηi)

        if isempty(ηi)
            integrated_qty = NaN
        else
            integrated_qty = mean(ηi)
        end

    elseif quantity == :ode_gas_r

        R = data_dict[:gas]["PARR"]
        filter!(!isnan, R)

        if isempty(R)
            integrated_qty = NaN
        else
            integrated_qty = mean(R)
        end

    elseif quantity == :ode_gas_cold_mf

        cold_fraction = data_dict[:gas]["COLF"]
        gas_mass = data_dict[:gas]["MASS"]

        filter!(!isnan, cold_fraction)
        filter!(!isnan, gas_mass)

        if isempty(cold_fraction) || isempty(gas_mass)
            integrated_qty = NaN
        else
            cold_mass = cold_fraction .* gas_mass
            integrated_qty = sum(cold_mass) ./ sum(gas_mass)
        end

    elseif quantity == :ode_stellar_it

        odit = data_dict[:stars]["ODIT"]
        filter!(!isnan, odit)

        if isempty(odit)
            integrated_qty = NaN
        else
            integrated_qty = mean(odit)
        end

    elseif quantity == :ode_stellar_accu_it

        acit = data_dict[:stars]["ACIT"]
        filter!(!isnan, acit)

        if isempty(acit)
            integrated_qty = NaN
        else
            integrated_qty = mean(acit)
        end

    elseif quantity == :ode_stellar_tau_s

        τS = data_dict[:stars]["TAUS"]
        filter!(!isnan, τS)

        if isempty(τS)
            integrated_qty = NaN
        else
            integrated_qty = mean(τS)
        end

    elseif quantity == :ode_stellar_eta_d

        ηd = data_dict[:stars]["ETAD"]
        filter!(!isnan, ηd)

        if isempty(ηd)
            integrated_qty = NaN
        else
            integrated_qty = mean(ηd)
        end

    elseif quantity == :ode_stellar_eta_i

        ηi = data_dict[:stars]["ETAI"]
        filter!(!isnan, ηi)

        if isempty(ηi)
            integrated_qty = NaN
        else
            integrated_qty = mean(ηi)
        end

    elseif quantity == :ode_stellar_r

        R = data_dict[:stars]["PARR"]
        filter!(!isnan, R)

        if isempty(R)
            integrated_qty = NaN
        else
            integrated_qty = mean(R)
        end

    elseif quantity == :ode_stellar_cold_mf

        cold_fraction = data_dict[:stars]["COLF"]
        gas_mass = data_dict[:stars]["GMAS"]

        filter!(!isnan, cold_fraction)
        filter!(!isnan, gas_mass)

        if isempty(cold_fraction) || isempty(gas_mass)
            integrated_qty = NaN
        else
            cold_mass = cold_fraction .* gas_mass
            integrated_qty = sum(cold_mass) ./ sum(gas_mass)
        end

    elseif quantity == :ode_stellar_gas_rho

        ρ = data_dict[:stars]["RHOC"]
        filter!(!isnan, ρ)

        if isempty(ρ)
            integrated_qty = NaN
        else
            integrated_qty = mean(ρ .* u"mp")
        end

    elseif quantity == :ode_stellar_gas_Z

        Z = data_dict[:stars]["PARZ"]
        gas_mass = data_dict[:stars]["GMAS"]

        filter!(!isnan, Z)
        filter!(!isnan, gas_mass)

        if isempty(Z) || isempty(gas_mass)
            integrated_qty = NaN
        else
            metal_mass = Z .* gas_mass
            integrated_qty = (sum(metal_mass) ./ sum(gas_mass)) ./ SOLAR_METALLICITY
        end

    elseif quantity == :ode_stellar_gas_mass

        gm = data_dict[:stars]["GMAS"]
        filter!(!isnan, gm)

        if isempty(gm)
            integrated_qty = NaN
        else
            integrated_qty = sum(gm)
        end

    elseif quantity == :ode_stellar_gas_sfr

        gsfr = data_dict[:stars]["GSFR"]
        filter!(!isnan, gsfr)

        if isempty(gsfr)
            integrated_qty = NaN
        else
            integrated_qty = sum(gsfr)
        end

    elseif quantity == :ode_stellar_gas_P

        P = data_dict[:stars]["GPRE"]
        filter!(!isnan, P)

        if isempty(P)
            integrated_qty = NaN
        else
            integrated_qty = mean(P)
        end

    else

        throw(ArgumentError("integrateQty: I don't recognize the quantity :$(quantity)"))

    end

    return integrated_qty

end

"""
    scatterQty(data_dict::Dict, quantity::Symbol)::Vector{<:Number}

Compute a quantity for each cell/particle in `data_dict`.

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

# Returns

  - The values of `quantity` for every cell/particle.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function scatterQty(data_dict::Dict, quantity::Symbol)::Vector{<:Number}

    if quantity == :stellar_mass

        scatter_qty = computeMass(data_dict, :stars)

    elseif quantity == :gas_mass

        scatter_qty = computeMass(data_dict, :gas)

    elseif quantity == :hydrogen_mass

        scatter_qty = computeMass(data_dict, :hydrogen)

    elseif quantity == :dm_mass

        scatter_qty = computeMass(data_dict, :dark_matter)

    elseif quantity == :bh_mass

        scatter_qty = computeMass(data_dict, :black_holes)

    elseif quantity == :molecular_mass

        scatter_qty = computeMass(data_dict, :molecular)

    elseif quantity == :br_molecular_mass

        scatter_qty = computeMass(data_dict, :br_molecular)

    elseif quantity == :atomic_mass

        scatter_qty = computeMass(data_dict, :atomic)

    elseif quantity == :ionized_mass

        scatter_qty = computeMass(data_dict, :ionized)

    elseif quantity == :neutral_mass

        scatter_qty = computeMass(data_dict, :neutral)

    elseif quantity == :stellar_gas_mass

        scatter_qty = computeMass(data_dict, :stellar)

    elseif quantity == :molecular_fraction

        scatter_qty = computeFraction(data_dict, :molecular)

    elseif quantity == :br_molecular_fraction

        scatter_qty = computeFraction(data_dict, :br_molecular)

    elseif quantity == :atomic_fraction

        scatter_qty = computeFraction(data_dict, :atomic)

    elseif quantity == :ionized_fraction

        scatter_qty = computeFraction(data_dict, :ionized)

    elseif quantity == :neutral_fraction

        scatter_qty = computeFraction(data_dict, :neutral)

    elseif quantity == :molecular_neutral_fraction

        fm = computeFraction(data_dict, :molecular)
        fn = computeFraction(data_dict, :neutral)

        scatter_qty = fm ./ fn

    elseif quantity == :ionized_neutral_fraction

        fi = computeFraction(data_dict, :ionized)
        fn = computeFraction(data_dict, :neutral)

        scatter_qty = fi ./ fn

    elseif quantity == :stellar_gas_fraction

        scatter_qty = computeFraction(data_dict, :stellar)

    elseif quantity == :mol_eq_quotient

        scatter_qty = computeEqQuotient(data_dict, :molecular)

    elseif quantity == :ion_eq_quotient

        scatter_qty = computeEqQuotient(data_dict, :ionized)

    elseif quantity == :gas_mass_density

        scatter_qty = computeVolumeDensity(data_dict, :gas)

    elseif quantity == :hydrogen_mass_density

        scatter_qty = computeVolumeDensity(data_dict, :hydrogen)

    elseif quantity == :gas_number_density

        scatter_qty = computeNumberDensity(data_dict, :gas)

    elseif quantity == :molecular_number_density

        scatter_qty = computeNumberDensity(data_dict, :molecular)

    elseif quantity == :br_molecular_number_density

        scatter_qty = computeNumberDensity(data_dict, :br_molecular)

    elseif quantity == :atomic_number_density

        scatter_qty = computeNumberDensity(data_dict, :atomic)

    elseif quantity == :ionized_number_density

        scatter_qty = computeNumberDensity(data_dict, :ionized)

    elseif quantity == :neutral_number_density

        scatter_qty = computeNumberDensity(data_dict, :neutral)

    elseif quantity == :gas_td

        scatter_qty = computeDepletionTime(computeMass(data_dict, :gas), data_dict[:gas]["SFR "])

    elseif quantity == :molecular_td

        scatter_qty = computeDepletionTime(
            computeMass(data_dict, :molecular),
            data_dict[:gas]["SFR "],
        )

    elseif quantity == :br_molecular_td

        scatter_qty = computeDepletionTime(
            computeMass(data_dict, :br_molecular),
            data_dict[:gas]["SFR "],
        )

    elseif quantity == :atomic_td

        scatter_qty = computeDepletionTime(
            computeMass(data_dict, :atomic),
            data_dict[:gas]["SFR "],
        )

    elseif quantity == :ionized_td

        scatter_qty = computeDepletionTime(
            computeMass(data_dict, :ionized),
            data_dict[:gas]["SFR "],
        )

    elseif quantity == :neutral_td

        scatter_qty = computeDepletionTime(
            computeMass(data_dict, :neutral),
            data_dict[:gas]["SFR "],
        )

    elseif quantity == :gas_metallicity

        scatter_qty = setPositive(data_dict[:gas]["GZ  "]) ./ SOLAR_METALLICITY

    elseif quantity == :stellar_metallicity

        scatter_qty = setPositive(data_dict[:stars]["GZ2 "]) ./ SOLAR_METALLICITY

    elseif quantity ∈ GAS_ABUNDANCE

        element_symbol = Symbol(first(split(string(quantity), "_")))

        abundances = computeAbundance(
            data_dict,
            :gas,
            element_symbol;
            solar=false,
        )

        if isempty(abundances)
            scatter_qty = Float64[]
        else
            scatter_qty = 12 .+ log10.(abundances)
            replace!(x -> isinf(x) ? NaN : x, scatter_qty)
        end

    elseif quantity ∈ STELLAR_ABUNDANCE

        element_symbol = Symbol(first(split(string(quantity), "_")))

        abundances = computeAbundance(
            data_dict,
            :stars,
            element_symbol;
            solar=false,
        )

        if isempty(abundances)
            scatter_qty = Float64[]
        else
            scatter_qty = 12 .+ log10.(abundances)
            replace!(x -> isinf(x) ? NaN : x, scatter_qty)
        end

    elseif quantity == :stellar_radial_distance

        scatter_qty = computeDistance(data_dict[:stars]["POS "])

    elseif quantity == :gas_radial_distance

        scatter_qty = computeDistance(data_dict[:gas]["POS "])

    elseif quantity == :dm_radial_distance

        scatter_qty = computeDistance(data_dict[:halo]["POS "])

    elseif quantity == :stellar_xy_distance

        if isempty(data_dict[:stars]["POS "])
            scatter_qty = eltype(data_dict[:stars]["POS "])[]
        else
            scatter_qty = computeDistance(data_dict[:stars]["POS "][1:2, :])
        end

    elseif quantity == :gas_xy_distance

        if isempty(data_dict[:gas]["POS "])
            scatter_qty = eltype(data_dict[:gas]["POS "])[]
        else
            scatter_qty = computeDistance(data_dict[:gas]["POS "][1:2, :])
        end

    elseif quantity == :dm_xy_distance

        if isempty(data_dict[:halo]["POS "])
            scatter_qty = eltype(data_dict[:halo]["POS "])[]
        else
            scatter_qty = computeDistance(data_dict[:halo]["POS "][1:2, :])
        end

    elseif quantity == :gas_sfr

        scatter_qty = data_dict[:gas]["SFR "]

    elseif quantity == :stellar_circularity

        (
            !logging[] ||
            @info("scatterQty: The stellar circularity depends on the positions and velocities of \
            all cell/particles. So, after filtering, the result for a given star will change")
        )

        scatter_qty = computeCircularity(data_dict)

    elseif quantity == :stellar_vcirc

       (
            !logging[] ||
            @info("scatterQty: The stellar circular velocity depends on the positions and \
            velocities of all cell/particles. So, after filtering, the result for a given star \
            will change")
       )

        _, scatter_qty = computeVcirc(data_dict)

    elseif quantity == :stellar_vradial

        scatter_qty = computeVpolar(data_dict, :radial)

    elseif quantity == :stellar_vtangential

        scatter_qty = computeVpolar(data_dict, :tangential)

    elseif quantity == :stellar_vzstar

        scatter_qty = computeVpolar(data_dict, :zstar)

    elseif quantity == :stellar_age

        scatter_qty = computeStellarAge(data_dict)

    elseif quantity == :sfr

        # Get the global index (index in the context of the whole simulation) of the current snapshot
        present_idx = data_dict[:snap_data].global_index

        if present_idx == 1

            scatter_qty = zeros(typeof(1.0u"Msun*yr^-1"), length(data_dict[:stars]["MASS"]))

        else

            # Get the physical times
            times = data_dict[:sim_data].table[:, 5]
            # Compute the time between snapshots
            Δt = times[present_idx] - times[present_idx - 1]

            scatter_qty = computeSFR(data_dict; age_resol=Δt)

        end

    elseif quantity == :ssfr

        # Get the global index (index in the context of the whole simulation) of the current snapshot
        present_idx = data_dict[:snap_data].global_index

        # Load the stellar masses
        stellar_masses = data_dict[:stars]["MASS"]

        if present_idx == 1 || iszero(stellar_mass)

            scatter_qty = zeros(typeof(1.0u"yr^-1"), length(stellar_masses))

        else

            # Get the physical times
            times = data_dict[:sim_data].table[:, 5]
            # Compute the time between snapshots
            Δt = times[present_idx] - times[present_idx - 1]

            scatter_qty = computeSFR(data_dict; age_resol=Δt) ./ stellar_masses

        end

    elseif quantity == :observational_sfr

        scatter_qty = computeSFR(data_dict; age_resol=AGE_RESOLUTION)

    elseif quantity == :observational_ssfr

        sfr = computeSFR(data_dict; age_resol=AGE_RESOLUTION)
        stellar_masses = data_dict[:stars]["MASS"]

        scatter_qty = sfr ./ stellar_masses

    elseif quantity == :stellar_eff

        scatter_qty = computeEfficiencyFF(
            data_dict[:stars]["RHOC"] .* u"mp",
            data_dict[:stars]["GMAS"],
            data_dict[:stars]["GSFR"],
        )

    elseif quantity == :gas_eff

        mass = computeMass(data_dict, :gas)

        scatter_qty = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

    elseif quantity == :molecular_eff

        mass = computeMass(data_dict, :molecular)

        scatter_qty = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

    elseif quantity == :br_molecular_eff

        mass = computeMass(data_dict, :br_molecular)

        scatter_qty = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

    elseif quantity == :atomic_eff

        mass = computeMass(data_dict, :atomic)

        scatter_qty = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

    elseif quantity == :ionized_eff

        mass = computeMass(data_dict, :ionized)

        scatter_qty = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

    elseif quantity == :neutral_eff

        mass = computeMass(data_dict, :neutral)

        scatter_qty = computeEfficiencyFF(data_dict[:gas]["RHO "], mass, data_dict[:gas]["SFR "])

    elseif quantity == :temperature

        scatter_qty = log10.(ustrip.(u"K", data_dict[:gas]["TEMP"]))
        replace!(x -> isinf(x) ? NaN : x, scatter_qty)

    elseif quantity == :pressure

        scatter_qty = data_dict[:gas]["PRES"]

    elseif quantity == :ode_gas_it

        scatter_qty = data_dict[:gas]["ODIT"]

    elseif quantity == :ode_gas_accu_it

        scatter_qty = data_dict[:gas]["ACIT"]

    elseif quantity == :ode_gas_tau_s

        scatter_qty = data_dict[:gas]["TAUS"]

    elseif quantity == :ode_gas_eta_d

        scatter_qty = data_dict[:gas]["ETAD"]

    elseif quantity == :ode_gas_eta_i

        scatter_qty = data_dict[:gas]["ETAI"]

    elseif quantity == :ode_gas_r

        scatter_qty = data_dict[:gas]["PARR"]

    elseif quantity == :ode_gas_cold_mf

        scatter_qty = data_dict[:gas]["COLF"]

    elseif quantity == :ode_stellar_it

        scatter_qty = data_dict[:stars]["ODIT"]

    elseif quantity == :ode_stellar_accu_it

        scatter_qty = data_dict[:stars]["ACIT"]

    elseif quantity == :ode_stellar_tau_s

        scatter_qty = data_dict[:stars]["TAUS"]

    elseif quantity == :ode_stellar_eta_d

        scatter_qty = data_dict[:stars]["ETAD"]

    elseif quantity == :ode_stellar_eta_i

        scatter_qty = data_dict[:stars]["ETAI"]

    elseif quantity == :ode_stellar_r

        scatter_qty = data_dict[:stars]["PARR"]

    elseif quantity == :ode_stellar_cold_mf

        scatter_qty = data_dict[:stars]["COLF"]

    elseif quantity == :ode_stellar_gas_rho

        scatter_qty = data_dict[:stars]["RHOC"] .* u"mp"

    elseif quantity == :ode_stellar_gas_Z

        scatter_qty = data_dict[:stars]["PARZ"] ./ SOLAR_METALLICITY

    elseif quantity == :ode_stellar_gas_mass

        scatter_qty = data_dict[:stars]["GMAS"]

    elseif quantity == :ode_stellar_gas_sfr

        scatter_qty = data_dict[:stars]["GSFR"]

    elseif quantity == :ode_stellar_gas_P

        scatter_qty = data_dict[:stars]["GPRE"]

    else

        throw(ArgumentError("scatterQty: I don't recognize the quantity :$(quantity)"))

    end

    if isempty(scatter_qty)
        return Number[]
    end

    return scatter_qty

end
