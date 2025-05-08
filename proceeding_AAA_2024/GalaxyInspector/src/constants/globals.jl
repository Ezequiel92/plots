####################################################################################################
# Constants and data structures
####################################################################################################

################
# Configuration
################

"""
If physical units (lengths) will be used, instead of comoving units (lengths).
"""
PHYSICAL_UNITS = false

"""
Relative path, within the simulation directory, of the `sfr.txt` file.
"""
const SFR_REL_PATH = "output/sfr.txt"

"""
Relative path, within the simulation directory, of the `cpu.txt` file.
"""
const CPU_REL_PATH = "output/cpu.txt"

"""
If logging messages will printed out.
"""
const logging = Ref(false)

########################
# Characteristic scales
########################

"""
Characteristic radius.
"""
const DISK_R = 40.0u"kpc"

"""
Characteristic box size.
"""
const BOX_L = 65.0u"kpc"

"""
Characteristic stellar age limit for the SFR and sSFR.
"""
const AGE_RESOLUTION = 200.0u"Myr"

######################
# Cell/particle types
######################

"""
Code index for each type of cell/particle.

# References

See for example Gadget2 [User's Guide](https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf), or Gadget4 [documentation](https://wwwmpa.mpa-garching.mpg.de/gadget4/).
"""
const LONG_PARTICLE_INDEX = Dict(
    :gas        => 0,
    :halo       => 1,
    :disk       => 2,
    :bulge      => 3,
    :stars      => 4,
    :black_hole => 5,
    :tracer     => 6,
)

"""
Human readable name for each type of cell/particle.
"""
const LONG_PARTICLE_NAMES = Dict(
    :gas        => "Gas cells",
    :halo       => "HR DM particles",
    :disk       => "IR DM particles",
    :bulge      => "LR DM particles",
    :stars      => "Stellar particles",
    :black_hole => "Black hole particles",
    :tracer     => "Tracer particles",
)

"""
Human readable name for each type of cell/particle.
"""
const ISOLATED_PARTICLE_NAMES = Dict(
    :gas   => "Gas cells",
    :halo  => "DM particles",
    :disk  => "Stellar disk",
    :bulge => "Stellar bulge",
)

"""
Human readable name for each morphological component.
"""
const MORPHOLOGICAL_COMPONENTS = Dict(
    :disk       => "Disk",
    :bulge      => "Bulge",
    :thin_disk  => "Thin disk",
    :thick_disk => "Thick disk",
)

############################
# Constants of the SF model
############################

# Star formation efficiency
const ϵff  = 1.0
# Recombination coefficient
const αH   = 2.6e-13u"cm^3 * s^-1"
# Formation rate coefficient of H₂ on dust grain (at solar metallicity)
const Rsun = 3.5e-17u"cm^3 * s^-1"
# Solar metallicity
const Zsun = 0.0127
# Effective metallicity
const Zeff = 1e-3 * Zsun
# Clumping factor
const Cρ   = 100.0

const C_star = sqrt(3π / 32u"G") / ϵff
const C_rec  = u"mp" / αH
const C_cond = (u"mp" * Zsun) / (2 * Rsun * Cρ)

#######################################
# Reference values from the literature
#######################################

"""
Internal unit of length used in IllustrisTNG, equivalent to ``1.0  \\, \\mathrm{kpc}``.
See the documentation [here](https://www.tng-project.org/data/docs/specifications/)
"""
const ILLUSTRIS_L_UNIT = 3.085678e21u"cm"

"""
Internal unit of mass used in IllustrisTNG, equivalent to ``10^{10} \\, \\mathrm{M_\\odot}``.
See the documentation [here](https://www.tng-project.org/data/docs/specifications/)
"""
const ILLUSTRIS_M_UNIT = 1.989e43u"g"

"""
Internal unit of velocity used in IllustrisTNG, equivalent to ``1.0 \\, \\mathrm{km \\, s^{-1}}``.
See the documentation [here](https://www.tng-project.org/data/docs/specifications/)
"""
const ILLUSTRIS_V_UNIT = 1.0e5u"cm*s^-1"

"""
Cosmological threshold density above which the gas cells/particles can turn into stars.

This value corresponds to `CritOverDensity` ``= 57.7 \\, [\\mathrm{cm^{-3}}]`` in the `param.txt` file (used only in cosmological simulations). Which is converted to internal units within the code using `OverDensThresh` = `CritOverDensity` * `OmegaBaryon` * 3 * `Hubble` * `Hubble` / (8 * `M_PI` * `G`)`. Then, to go to physical units again one has to do: `OverDensThresh` * `UnitDensity_in_cgs` * `cf_a3inv` * `HubbleParam` * `HubbleParam`.

Using the unit factors,

`UnitLength_in_cm`         = ``3.085678 \\times 10^{24}``

`UnitMass_in_g`            = ``1.989 \\times 10^{43}``

`UnitVelocity_in_cm_per_s` = ``100000``

The derived units,

`UnitTime_in_s`      = `UnitLength_in_cm` * `UnitVelocity_in_cm_per_s`^-1 = ``3.08568 \\times 10^{19}``

`UnitDensity_in_cgs` = `UnitMass_in_g` * `UnitLength_in_cm^-3`            = ``6.76991 \\times 10^{-31}``

The parameters,

`OmegaBaryon`       = ``0.048``

`HubbleParam`       = ``0.6777``

`PROTONMASS`        = ``1.67262178 \\times 10^{-24}``

`HYDROGEN_MASSFRAC` = ``0.76``

`GRAVITY`           = ``6.6738 \\times 10^{-8}``

`HUBBLE`            = ``3.2407789 \\times 10^{-18}``

And the derived parameters,

Hubble = `HUBBLE` * `UnitTime_in_s`                                              = ``100``

G      = `GRAVITY` * `UnitLength_in_cm`^-3 * `UnitMass_in_g` * `UnitTime_in_s`^2 = ``43.0187``

One gets,

`OverDensThresh` = 76.8495 [internal units of density]

And, for a cosmological simulation at readshift 0 (`cf_a3inv` = 1), this result in a physical density threshold of ``1.42857 \\times 10^{-5} \\, [\\mathrm{cm^{-3}}]``, or adding the proton mass a value of:

``\\log_{10} \\rho \\ [\\mathrm{M_\\odot \\, kpc^{-3}}] = 2.548``

"""
const COSMO_THRESHOLD_DENSITY = 353.059u"Msun*kpc^-3"

"""
Threshold density above which the gas cells/particles can turn into stars.

This value corresponds to `CritPhysDensity` ``= 0.318 \\, [\\mathrm{cm^{-3}}]`` in the `param.txt` file (used in cosmological and non-cosmological simulations). Which is converted to internal units within the code using `PhysDensThresh` = `CritPhysDensity` * `PROTONMASS` / `HYDROGEN_MASSFRAC` / `UnitDensity_in_cgs`. Then, to go to physical units again one has to do: `PhysDensThresh` * `UnitDensity_in_cgs` * `cf_a3inv` * `HubbleParam` * `HubbleParam`.

`PhysDensThresh` = ``1.03378 \\times 10^{6}`` [internal units of density]

For a cosmological simulation at readshift 0 (`cf_a3inv` = 1), this result in a physical density threshold of ``0.192 \\, [\\mathrm{cm^{-3}}]``, or adding the proton mass a value of:

``\\log_{10} \\rho \\, [\\mathrm{M_\\odot \\, kpc^{-3}}] = 6.677``
"""
const THRESHOLD_DENSITY = 4.749326e6u"Msun*kpc^-3"

@doc raw"""
Hubble constant in $\mathrm{Gyr^{-1}}$.

This value corresponds to $H_0 = 0.102201 \, \mathrm{Gyr}^{-1} = 100 \, \mathrm{km} \, \mathrm{s}^{-1} \, \mathrm{Mpc}^{-1}$.
"""
const HUBBLE_CONSTANT = 0.102201

"""
Subhalo numbers for the MW and M31 in Hestia simulations.
"""
const HESTIA_SUBHALOS = Dict(
    "Hestia17-11" => Dict(
        :subhalo_number_MW  => 1,
        :subhalo_number_M31 => 0,
    ),
    "Hestia09-18" => Dict(
        :subhalo_number_MW  => 3911,
        :subhalo_number_M31 => 2608,
    ),
    "Hestia37-11" => Dict(
        :subhalo_number_MW  => 920,
        :subhalo_number_M31 => 0,
    ),
)

@doc raw"""
Solar abundances.

They are defined as $12 + \log_{10}(N_\mathrm{X} / N_\mathrm{H})$, where $N_\mathrm{X}$ and $N_\mathrm{H}$ are the number densities of element $\mathrm{X}$ and hydrogen respectively.

# References

M. Asplund et al. (2009). *The Chemical Composition of the Sun*. Annual Review of Astronomy and Astrophysics, **47(1)**, 481–522. [doi:10.1146/annurev.astro.46.060407.145222](https://doi.org/10.1146/annurev.astro.46.060407.145222)
"""
const SOLAR_ABUNDANCE = Dict(
    :H  => 12,    # Hydrogen
    :He => 10.93, # Helium
    :C  => 8.43,  # Carbon
    :N  => 7.83,  # Nitrogen
    :O  => 8.69,  # Oxygen
    :Ne => 7.93,  # Neon
    :Mg => 7.60,  # Magnesium
    :Si => 7.51,  # Silicon
    :S  => 7.12,  # Sulfur
    :Ca => 6.34,  # Calcium
    :Fe => 7.50,  # Iron
)

"""
Standard atomic weights.

# References

T. Prohaska et al. (2022). *Standard atomic weights of the elements 2021 (IUPAC Technical Report)*. Pure and Applied Chemistry, **94(5)**, 573-600. [doi:10.1515/pac-2019-0603](https://doi.org/10.1515/pac-2019-0603)
"""
const ATOMIC_WEIGHTS = Dict(
    :H  => 1.0080, # Hydrogen
    :He => 4.0026, # Helium
    :C  => 12.011, # Carbon
    :N  => 14.007, # Nitrogen
    :O  => 15.999, # Oxygen
    :Ne => 20.180, # Neon
    :Mg => 24.305, # Magnesium
    :Si => 28.085, # Silicon
    :S  => 32.06,  # Sulfur
    :Ca => 40.078, # Calcium
    :Fe => 55.845, # Iron
)

@doc raw"""
Kennicutt-Schmidt law fits for molecular and neutral gas, from Bigiel et al. (2008) (Table 2, Average).

Power-law index, N, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $10 \, \mathrm{M_\odot \, pc^{-2}}$ are given.

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{HI, H_2, gas}}{10 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
const A_BIGIEL2008_MOLECULAR = −2.06 ± 0.17
const N_BIGIEL2008_MOLECULAR = 0.96 ± 0.07
const A_BIGIEL2008_NEUTRAL   = −2.39 ± 0.28
const N_BIGIEL2008_NEUTRAL   = 1.85 ± 0.70

@doc raw"""
Kennicutt-Schmidt law best-fit for molecular gas, from Bigiel et al. (2008) (Section 4.3, Equation 3).

Power-law index, N, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $10 \, \mathrm{M_\odot \, pc^{-2}}$ are given.

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{H_2}}{10 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
const A_BIGIEL2008_BF_MOLECULAR = −2.1 ± 0.2
const N_BIGIEL2008_BF_MOLECULAR = 1.0 ± 0.2

"""
Spatial resolution used in Bigiel et al. (2008).

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
const BIGIEL_PX_SIZE = 750.0u"pc"

"""
Spatial resolution used in Sun et al. (2023).

# References

J. Sun et al. (2023). *Star Formation Laws and Efficiencies across 80 Nearby Galaxies*. The Astrophysical Journal Letters, **945(2)**, L19. [doi:10.3847/2041-8213/acbd9c](https://doi.org/10.3847/2041-8213/acbd9c)
"""
const SUN_PX_SIZE = 1.5u"kpc"

@doc raw"""
Range of values for

```math
\Sigma_\mathrm{SFR} \, [\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}] \, ,
```
in the seven spiral in Table 1 of Bigiel et al. (2008), with associated molecular data.

The actual values for the SFR density are taken from Table 2 in Bigiel et al. (2010), using only the ones with associated molecular data.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)

F. Bigiel et al. (2010). *EXTREMELY INEFFICIENT STAR FORMATION IN THE OUTER DISKS OF NEARBY GALAXIES*. The Astrophysical Journal, **140(5)**, 1194. [doi:10.1088/0004-6256/140/5/1194](https://doi.org/10.1088/0004-6256/140/5/1194)
"""
const BIGIEL2008_SFR_RANGE = exp10.([-2.99, -0.33]) .* u"Msun * yr^-1 * kpc^-2"

"""
Path to the file with Table 2 from Bigiel et al. (2010).

# References

F. Bigiel et al. (2010). *EXTREMELY INEFFICIENT STAR FORMATION IN THE OUTER DISKS OF NEARBY GALAXIES*. The Astrophysical Journal, **140(5)**, 1194. [doi:10.1088/0004-6256/140/5/1194](https://doi.org/10.1088/0004-6256/140/5/1194)
"""
const BIGIEL2010_TABLE_2 = joinpath(
    @__DIR__,
    "../../experimental_data/bigiel_2010/aj351668t2_mrt.txt",
)

"""
Path to the file with Table 3 from Bigiel et al. (2010).

# References

F. Bigiel et al. (2010). *EXTREMELY INEFFICIENT STAR FORMATION IN THE OUTER DISKS OF NEARBY GALAXIES*. The Astrophysical Journal, **140(5)**, 1194. [doi:10.1088/0004-6256/140/5/1194](https://doi.org/10.1088/0004-6256/140/5/1194)
"""
const BIGIEL2010_TABLE_3 = joinpath(
    @__DIR__,
    "../../experimental_data/bigiel_2010/aj351668t3_mrt.txt",
)

"""
Path to the file with Table A1 from Sun et al. (2023).

# References

J. Sun et al. (2023). *Star Formation Laws and Efficiencies across 80 Nearby Galaxies*. The Astrophysical Journal Letters, **945(2)**, L19. [doi:10.3847/2041-8213/acbd9c](https://doi.org/10.3847/2041-8213/acbd9c)
"""
const SUN2023_TABLE = joinpath(
    @__DIR__,
    "../../experimental_data/sun_2023/apjlacbd9ct3_mrt.txt",
)

@doc raw"""
Range of values for

```math
\Sigma_\mathrm{SFR} \, [\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}] \, ,
```
from the combine data (Table 1 and 2) in Kennicutt (1998).

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
const KS98_SFR_RANGE = exp10.([-3.55, 2.98]) .* u"Msun * yr^-1 * kpc^-2"

@doc raw"""
Slope of the Kennicutt-Schmidt law, taken from Kennicutt (1998) (Section 4, Equation 4).

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{gas}}{1 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \mathrm{M_\odot \, yr^{-1} \, kpc^{-2}} \, ,
```

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
const N_KS98 = 1.4 ± 0.15

@doc raw"""
Intercept of the Kennicutt-Schmidt law, taken from Kennicutt (1998) (Section 4, Equation 4).

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{gas}}{1 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \mathrm{M_\odot \, yr^{-1} \, kpc^{-2}} \, ,
```

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
const a_KS98 = 2.5e-4 ± 0.7e-4

"""
Path to the file with the Milky Way profiles from Mollá et al. (2015).

# References

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function*. Monthly Notices of the Royal Astronomical Society **451(4)**, 3693–3708. [doi:10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)
"""
const MOLLA2015_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/mollá_2015.csv")

"""
Path to the file with the global galactic properties from Feldmann (2020).

# References

R. Feldmann (2020). *The link between star formation and gas in nearby galaxies*. Communications Physics **3(226)**. [doi:10.1038/s42005-020-00493-0](https://doi.org/10.1038/s42005-020-00493-0)
"""
const FELDMANN2020_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/feldmann_2020.csv")

"""
Path to the file with the data from Leroy et al. (2008).

# References

A. K. Leroy et al. (2008). *THE STAR FORMATION EFFICIENCY IN NEARBY GALAXIES: MEASURING WHERE GAS FORMS STARS EFFECTIVELY*. The Astronomical Journal **136(6)**, 2782–2845. [doi:10.1088/0004-6256/136/6/2782](https://doi.org/10.1088/0004-6256/136/6/2782)

"""
const LEROY2008_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/leroy_2008.jld2")

"""
Path to the file with the fits from McMillan (2011).

# References

P. J. McMillan (2011). *Mass models of the Milky Way*. Monthly Notices of the Royal Astronomical Society **414(3)**, 2446–2457. [doi:10.1111/j.1365-2966.2011.18564.x](https://doi.org/10.1111/j.1365-2966.2011.18564.x)
"""
const MCMILLAN2011_DATA_PATH = joinpath(@__DIR__, "../../experimental_data/mcmillan_2011.jld2")

"""
Reference pressure for the molecular fraction-pressure relation, taken from Blitz et al. (2006) (Table 2, "Mean" row, Third column).

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
const P0 = 3.5e4u"K*cm^-3" * Unitful.k

"""
Reference exponent for the molecular fraction-pressure relation, taken from Blitz et al. (2006) (Table 2, "Mean" row, Second column).

We use -α here.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
const ALPHA_BLITZ = -0.92

###############
# Type aliases
###############

"""
Color type.
"""
const ColorType = Union{ColorTypes.RGB,ColorTypes.RGBA,Symbol}

"""
Line style type.
"""
const LineStyleType = Union{Tuple{String,Symbol},Nothing,String,Symbol}

"""
Index type.
"""
const IndexType = Union{
    Colon,
    Integer,
    UnitRange{<:Integer},
    StepRange{<:Integer,<:Integer},
    Vector{<:Integer},
}

"""
Reduced index type.
"""
const ReducedIndexType = Union{
    Integer,
    UnitRange{<:Integer},
    StepRange{<:Integer,<:Integer},
    Vector{<:Integer},
}

# Dimensions of specific energy
@derived_dimension SpecificEnergy Unitful.𝐋^2 * Unitful.𝐓^-2 true

# Dimensions of surface density
@derived_dimension SurfaceDensity Unitful.𝐌 * Unitful.𝐋^-2 true

# Dimensions of mass flow surface density
@derived_dimension MassFlowDensity Unitful.𝐌 * Unitful.𝐓^-1 * Unitful.𝐋^-2 true

# Dimensions of angular momentum
@derived_dimension AngularMomentum Unitful.𝐌 * Unitful.𝐋^2 * Unitful.𝐓^-1 true

# Dimensions of number density
@derived_dimension NumberDensity Unitful.𝐋^-3 true

#########################
# Makie.jl default theme
#########################

"""
Default list of marker types.
"""
const MARKERS = [
    :circle,
    :rect,
    :dtriangle,
    :utriangle,
    :cross,
    :diamond,
    :ltriangle,
    :rtriangle,
    :pentagon,
    :xcross,
    :hexagon,
]

"""
Default list of line styles.
"""
const LINE_STYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]

"""
Default cycler.
"""
const CYCLE = Cycle([:color, :linestyle, :marker], covary=true)

"""
Default plot theme.

Regarding the graphic units used, we know that ``1 \\, \\mathrm{mm} = 2.83466 \\, \\mathrm{pt}`` and ``1 \\, \\mathrm{in} = 25.4 \\, \\mathrm{mm}``. Then, if we want ``1 \\, \\mathrm{[code\\,\\,]unit} = 0.1 \\, \\mathrm{mm}`` in vector graphics, we have to use `pt_per_unit` = 0.283466.

For pixel images, we control the ppi with `px_per_unit`. A resonable high ppi is 600. So, using `px_per_unit` = ``2.3622`` we get ``23.622 \\, \\mathrm{px/mm} \\sim 600 \\, \\mathrm{px/in}`` (remember that ``1 \\, \\mathrm{[code\\,\\,]unit} = 0.1 \\, \\mathrm{mm}``).
"""
const DEFAULT_THEME = Theme(
    #################################################################################################
    # Size of the figures in code units
    # For PDFs and SVGs, 880 [code ]unit = 8.8 cm
    # For PNGs, when printed to a size of 1 point = 0.1 mm, one will get a dpi of 600 (23.622 px/mm)
    #################################################################################################
    size=(880, 880),
    ######################################
    # 35 unit * 0.283466 pt/unit ~ 9.9 pt
    ######################################
    fontsize=35,
    #############################
    # (left, right, bottom, top)
    #############################
    figure_padding=(1, 15, 5, 15),
    palette=(color=Makie.wong_colors(), marker=MARKERS, linestyle=LINE_STYLES),
    CairoMakie=(px_per_unit=2.3622, pt_per_unit=0.283466),
    Axis=(
        xlabelpadding=15,
        xticklabelpad=10,
        xticksize=7,
        xgridvisible=false,
        spinewidth=3,
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(5),
        ylabelpadding=15,
        yticklabelpad=10,
        yticksize=7,
        ygridvisible=false,
        yminorticksvisible=true,
        yminorticks=IntervalsBetween(5),
        ######################################################################################
        # Aspect ratio of the figures. The options are:
        # nothing: The aspect ratio will be chosen by [Makie](https://docs.makie.org/stable/)
        # AxisAspect(n): The aspect ratio will be given by the number `n` = width / height
        # DataAspect(): The aspect ratio of the data will be used
        ######################################################################################
        aspect=AxisAspect(1),
    ),
    Legend=(
        tellheight=false,
        tellwidth=false,
        framevisible=false,
        colgap=20,
        halign=:right,
        valign=:bottom,
        nbanks=3,
        titlegap=-5,
        labelsize=30,
        linewidth=5,
        markersize=28,
        patchsize=(50, 50),
        linepoints=[Point2f(0.0, 0.5), Point2f(0.9, 0.5)],
        ###############################################
        # Vertices, relative to the default 1x1 square
        ###############################################
        polypoints=[
            Point2f(0.15, 0.15),
            Point2f(0.85, 0.15),
            Point2f(0.85, 0.85),
            Point2f(0.15, 0.85),
        ],
    ),
    Lines=(linewidth=5, cycle=CYCLE),
    VLines=(linewidth=3, cycle=CYCLE),
    HLines=(linewidth=3, cycle=CYCLE),
    ScatterLines=(
        linewidth=5,
        markersize=22,
        cycle=CYCLE,
    ),
    Scatter=(markersize=22, cycle=CYCLE),
    Band=(cycle=CYCLE, alpha=0.5),
    Errorbars=(whiskerwidth=10,),
    ########################################################################
    # Alternative colormaps:
    # colormap = :nipy_spectral - nan_color = ColorSchemes.nipy_spectral[1]
    # colormap = :cubehelix     - nan_color = ColorSchemes.cubehelix[1]
    ########################################################################
    Heatmap=(colormap=:CMRmap, nan_color=ColorSchemes.CMRmap[1]),
    Colorbar=(
        colormap=:CMRmap,
        size=25,
        ticklabelpad=10,
        minorticksvisible=true,
        ticksize=7,
        labelpadding=2,
    ),
    BarPlot=(
        color_over_background=:black,
        color_over_bar=:black,
        flip_labels_at=10,
        direction=:x,
        strokecolor=:black,
        strokewidth=1,
        bar_labels=:y,
        dodge_gap=0.04,
    ),
    Arrows=(lengthscale=0.02, arrowsize=7.0, linestyle=:solid, color=:white),
    Hist=(strokecolor=:black, strokewidth=1),
)

#############
# Structures
#############

"""
Dimensional information about a physical quantity.

# Fields

  - `hdf5_name::String`: HDF5 block name.
  - `dimensions::Unitful.Dimensions`: Physical dimensions of the quantity, e.g. `Unitful.𝐋 * Unitful.𝐓^-1`.
  - `unit::Union{Unitful.Units,Symbol}`: Units of the quantity within the simulation code. It can be a unit from [Unitful](https://github.com/PainterQubits/Unitful.jl) or [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl), or it can be the symbol `:internal` which denotes internal code units.
"""
struct Qty
    hdf5_name::String
    dimensions::Unitful.Dimensions
    unit::Union{Unitful.Units,Symbol}
end

"""
Data in the "Header" group of a HDF5 snapshot file.

# Fields

  - `box_size::Float64`: Total size of the simulation box.
  - `h0::Float64`: Hubble parameter.
  - `mass_table::Vector{Float64}`: Masses of particle types which have a constant mass.
  - `num_files::Int32`: Number of file chunks per snapshot.
  - `num_part::Vector{Int32}`: Number of particles (of each type) included in this file chunk.
  - `num_total::Vector{UInt32}`: Total number of particles (of each type) for this snapshot.
  - `omega_0::Float64`: The cosmological density parameter for matter.
  - `omega_l::Float64`: The cosmological density parameter for the cosmological constant.
  - `redshift::Float64`: The redshift.
  - `time::Float64`: The physical time/scale factor.
  - `l_unit::Unitful.Length`: Conversion factor from internal units of length to centimeters.
  - `m_unit::Unitful.Mass`: Conversion factor from internal units of mass to grams.
  - `v_unit::Unitful.Velocity`: Conversion factor from internal units of velocity to centimeters per second.
"""
@kwdef mutable struct SnapshotHeader
    box_size::Float64
    h0::Float64
    mass_table::Vector{Float64}
    num_files::Int32
    num_part::Vector{Int32}
    num_total::Vector{UInt32}
    omega_0::Float64
    omega_l::Float64
    redshift::Float64
    time::Float64
    l_unit::Unitful.Length
    m_unit::Unitful.Mass
    v_unit::Unitful.Velocity
end

"""
Data in the "Header" group of a HDF5 group catalog file.

Default values are for when there are no group catalog files.

# Fields

  - `box_size::Float64 = NaN`: Total size of the simulation box.
  - `h0::Float64 = NaN`: Hubble parameter.
  - `n_groups_part::Int32 = -1`: Number of halos (FoF groups) in this file chunk.
  - `n_groups_total::Int32 = -1`: Total number of halos (FoF groups) in this snapshot.
  - `n_subgroups_part::Int32 = -1`: Number of subhalos (subfind) in this file chunk.
  - `n_subgroups_total::Int32 = -1`: Total number of subhalos (subfind) in this snapshot.
  - `num_files::Int32 = -1`: Number of file chunks per snapshot.
  - `omega_0::Float64 = NaN`: The cosmological density parameter for matter.
  - `omega_l::Float64 = NaN`: The cosmological density parameter for the cosmological constant.
  - `redshift::Float64 = NaN`: The redshift.
  - `time::Float64 = NaN`: The physical time/scale factor.
"""
@kwdef mutable struct GroupCatHeader
    box_size::Float64 = NaN
    h0::Float64 = NaN
    n_groups_part::Int32 = -1
    n_groups_total::Int32 = -1
    n_subgroups_part::Int32 = -1
    n_subgroups_total::Int32 = -1
    num_files::Int32 = -1
    omega_0::Float64 = NaN
    omega_l::Float64 = NaN
    redshift::Float64 = NaN
    time::Float64 = NaN
end

"""
Metadata for a simulation.

# Fields

  - `path::String`: Full path to the simulation directory.
  - `index::Int`: An index associated with the simulation.
  - `slice::IndexType`: Slice of the simulation, i.e. which snapshots will be read. It can be an integer (a single snapshot), a vector of integers (several snapshots), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (all snapshots).
  - `cosmological::Bool`: If the simulation is cosmological,

      + `false` -> Newtonian simulation    (`ComovingIntegrationOn` = 0).
      + `true`  -> Cosmological simulation (`ComovingIntegrationOn` = 1).
  - `table::DataFrame`: A dataframe where each row is a snapshot, and the following 8 colums:

      + `:ids`            -> Dataframe index of each snapshot, i.e. if there are 10 snapshots in total it runs from 1 to 10.
      + `:numbers`        -> Number in the file name of each snapshot.
      + `:scale_factors`  -> Scale factor of each snapshot.
      + `:redshifts`      -> Redshift of each snapshot.
      + `:physical_times` -> Physical time since the Big Bang.
      + `:lookback_times` -> Physical time left to reach the last snapshot.
      + `:snapshot_paths` -> Full path to each snapshots.
      + `:groupcat_paths` -> Full path to each group catalog files.
"""
struct Simulation
    path::String
    index::Int
    slice::IndexType
    cosmological::Bool
    table::DataFrame
end

"""
Metadata for a snapshot.

# Fields

  - `path::String`: Full path to the snapshot.
  - `global_index::Int`: Index of the snapshot in the context of the whole simulation.
  - `slice_index::Int`: Index of the snapshot in the context of the slice.
  - `physical_time::Unitful.Time`: Physical time since the Big Bang.
  - `lookback_time::Unitful.Time`: Physical time left to reach the last snapshot.
  - `scale_factor::Float64`: Scale factor of the snapshot.
  - `redshift::Float64`: Redshift of the snapshot.
  - `header::SnapshotHeader`: Header of the snapshot.
"""
struct Snapshot
    path::String
    global_index::Int
    slice_index::Int
    physical_time::Unitful.Time
    lookback_time::Unitful.Time
    scale_factor::Float64
    redshift::Float64
    header::SnapshotHeader
end

"""
Metadata for a group catalog file.

# Fields

  - `path::Union{String,Missing}`: Full path to the group catalog file.
  - `header::GroupCatHeader`: Header of the group catalog.
"""
struct GroupCatalog
    path::Union{String,Missing}
    header::GroupCatHeader
end

"""
Unit conversion factors.

# Fields

  - `x_cgs::Unitful.Length`: Length, from internal units to ``\\mathrm{cm}``.
  - `x_cosmo::Unitful.Length`: Length, from internal units to ``\\mathrm{kpc}``.
  - `x_comoving::Unitful.Length`: Length, from internal units to ``\\mathrm{ckpc}``.
  - `v_cgs::Unitful.Velocity`: Velocity, from internal units to ``\\mathrm{cm \\, s^{-1}}``.
  - `v_cosmo::Unitful.Velocity`: Velocity, from internal units to ``\\mathrm{km \\, s^{-1}}``.
  - `m_cgs::Unitful.Mass`: Mass, from internal units to ``\\mathrm{g}``.
  - `m_cosmo::Unitful.Mass`: Mass, from internal units to ``\\mathrm{M_\\odot}``.
  - `t_cgs::Unitful.Time`: Time, from internal units to ``\\mathrm{s}``.
  - `t_cosmo::Unitful.Time`: Time, from internal units to ``\\mathrm{Myr}``.
  - `U_cgs::Unitful.Energy`: Specific energy, from internal units to ``\\mathrm{erg \\, g^{-1}}``.
  - `rho_cgs::Unitful.Density`: Density, from internal units to ``\\mathrm{g \\, cm^{-3}}``.
  - `P_Pa::Unitful.Pressure`: Pressure, from internal units to ``\\mathrm{Pa}``.
"""
struct InternalUnits

    x_cgs::Unitful.Length      # Length, from internal units to cm
    x_cosmo::Unitful.Length    # Length, from internal units to kpc
    x_comoving::Unitful.Length # Length, from internal units to ckpc

    v_cgs::Unitful.Velocity    # Velocity, from internal units to cm * s^-1
    v_cosmo::Unitful.Velocity  # Velocity, from internal units to km * s^-1

    m_cgs::Unitful.Mass        # Mass, from internal units to g
    m_cosmo::Unitful.Mass      # Mass, from internal units to M⊙

    t_cgs::Unitful.Time        # Time, from internal units to s
    t_cosmo::Unitful.Time      # Time, from internal units to Myr

    U_cgs::SpecificEnergy      # Specific energy, from internal units to erg * g^-1

    rho_cgs::Unitful.Density   # Density, from internal units to g * cm^-3

    P_Pa::Unitful.Pressure     # Pressure, from internal units to Pa

    """
        InternalUnits(; <keyword arguments>)

    Constructor for `InternalUnits`.

    For cosmological simulations, and when appropriate, the conversion factors turn comoving units into physical ones.

    # Arguments

      - `l_unit::Unitful.Length=ILLUSTRIS_L_UNIT`: Code parameter `UnitLength_in_cm`.
      - `m_unit::Unitful.Mass=ILLUSTRIS_M_UNIT`: Code parameter `UnitMass_in_g`.
      - `v_unit::Unitful.Velocity=ILLUSTRIS_V_UNIT`: Code parameter `UnitVelocity_in_cm_per_s`.
      - `a0::Float64=1.0`: Cosmological scale factor of the simulation.
      - `h0::Float64=1.0`: Hubble constant as "little h".
    """
    function InternalUnits(;
        l_unit::Unitful.Length=ILLUSTRIS_L_UNIT,
        m_unit::Unitful.Mass=ILLUSTRIS_M_UNIT,
        v_unit::Unitful.Velocity=ILLUSTRIS_V_UNIT,
        a0::Float64=1.0,
        h0::Float64=1.0,
    )

        #############
        # Base units
        #############

        # Length conversion factors
        x_cgs = l_unit * a0 / h0
        x_cosmo = x_cgs |> u"kpc"
        x_comoving = l_unit / h0 |> u"kpc"

        # Velocity conversion factors
        v_cgs = v_unit * sqrt(a0)
        v_cosmo = v_cgs |> u"km*s^-1"

        # Mass conversion factors
        m_cgs = m_unit / h0
        m_cosmo = m_cgs |> u"Msun"

        ################
        # Derived units
        ################

        # Temperature conversion factors
        t_cgs = x_cgs / v_cgs
        t_cosmo = t_cgs |> u"Myr"

        # Specific energy conversion factor
        U_cgs = v_unit^2 |> u"erg*g^-1"

        # Density conversion factor
        rho_cgs = m_cgs * x_cgs^-3

        # Thermal pressure conversion factor (it uses v_unit^2 instead of v_cgs^2,
        # which would add an extra factor of a0)
        P_Pa = v_unit^2 * m_cgs * x_cgs^-3 |> u"Pa"

        new(
            x_cgs,
            x_cosmo,
            x_comoving,
            v_cgs,
            v_cosmo,
            m_cgs,
            m_cosmo,
            t_cgs,
            t_cosmo,
            U_cgs,
            rho_cgs,
            P_Pa,
        )

    end

end

"""
Linear grid (1D).

# Fields

  - `grid::Vector{<:Number}`: Vector with the central value of each bin.
  - `ticks::Vector{<:Number}`: Vector with the edges of the bins.
  - `bin_widths::Vector{<:Number}`: Widths of the bins.
  - `log::Bool`: If the grid is logarithmic.
"""
struct LinearGrid
    grid::Vector{<:Number}
    ticks::Vector{<:Number}
    bin_widths::Vector{<:Number}
    log::Bool

    """
        LinearGrid(start::Number, stop::Number, n_bins::Int; <keyword arguments>)

    Constructor for `LinearGrid`.

    # Arguments

      - `start::Number`: Initial value of the grid.
      - `stop::Number`: Final value of the grid.
      - `n_bins::Int`: Number of bins.
      - `log::Bool=false`: If the bins will be logarithmic.
    """
    function LinearGrid(start::Number, stop::Number, n_bins::Int; log::Bool=false)

        # Example of a grid with 3 bins:
        #
        # +-----------------------------+-----------------------------+-----------------------------+
        # |   grid[1] -->|              |   grid[2] -->|              |   grid[3] -->|              |
        # +-----------------------------+-----------------------------+-----------------------------+
        #
        # +-----------------------------+-----------------------------+-----------------------------+
        # |<-- tick[1] (= `start`)      |<-- tick[2]       tick[3] -->|       tick[4] (= `stop`) -->|
        # +-----------------------------+-----------------------------+-----------------------------+

        (
            stop > start ||
            throw(ArgumentError("LinearGrid: `stop` must be larger than `start`, \
            but I got `stop` = $(stop) <= `start` = $(start)"))
        )

        if log

            (
                isPositive(start) ||
                throw(ArgumentError("LinearGrid: For a logarithmic grid you need a strictly \
                positive `start`, but I got `start` = $(start)"))
            )

            # Unit of length
            u_l = unit(start)

            log_start = log10(ustrip(start))
            log_stop  = log10(ustrip(u_l, stop))

            width      = (log_stop - log_start) / n_bins
            grid       = [exp10((i - 0.5) * width + log_start) * u_l for i in 1:n_bins]
            ticks      = [exp10(i * width + log_start) * u_l for i in 0:n_bins]
            bin_widths = [ticks[i + 1] - ticks[i] for i in 1:n_bins]

        else

            width      = (stop - start) / n_bins
            grid       = [(i - 0.5) * width + start for i in 1:n_bins]
            ticks      = [i * width + start for i in 0:n_bins]
            bin_widths = [width for _ in 1:n_bins]

        end

        new(grid, ticks, bin_widths, log)

    end
end

"""
Square grid (2D).

# Fields

  - `grid::Matrix{NTuple{2,<:Number}}`: Matrix with the physical coordinates of each pixel in the grid.
  - `x_ticks::Vector{<:Number}`: Full set of possible values for the x coordinate.
  - `y_ticks::Vector{<:Number}`: Full set of possible values for the y coordinate.
  - `physical_size::Number`: Side length of the square grid.
  - `n_bins::Int`: Number of bins per side of the grid.
  - `bin_width::Number`: Side length of each bin.
  - `bin_area::Number`: Area of each bin.
"""
struct SquareGrid
    grid::Matrix{NTuple{2,<:Number}}
    x_ticks::Vector{<:Number}
    y_ticks::Vector{<:Number}
    physical_size::Number
    n_bins::Int
    bin_width::Number
    bin_area::Number

    """
        SquareGrid(
            physical_size::Number,
            n_bins::Int;
            <keyword arguments>
        )

    Constructor for `SquareGrid`.

    # Arguments

      - `physical_size::Number`: Side length of the square grid.
      - `n_bins::Int`: Number of bins per dimesion of the grid.
      - `center::Vector{<:Number}=zeros(typeof(physical_size),3)`: 3D location of the center of the grid. The z axis is taken as the normal vector of the grid.
    """
    function SquareGrid(
        physical_size::Number,
        n_bins::Int;
        center::Vector{<:Number}=zeros(typeof(physical_size), 3),
    )

        # Example of a 3x3 grid:
        #
        # x_ticks       = [0, 1, 2]
        # y_ticks       = [0, 1, 2]
        # n_bins        = 3
        # physical_size = 3
        # center        = [1, 1, 1]
        #
        # Reference system:
        #
        #   y
        #   ↑
        #   |
        #   |
        #   |
        #    -------→ x
        #
        # +------------------+------------------+------------------+
        # |      i = 1       |      i = 4       |      i = 7       |
        # | grid[1] = (0, 2) | grid[4] = (1, 2) | grid[7] = (2, 2) |
        # +------------------+------------------+------------------+
        # |      i = 2       |      i = 5       |      i = 8       |
        # | grid[2] = (0, 1) | grid[5] = (1, 1) | grid[8] = (2, 1) |
        # +------------------+------------------+------------------+
        # |      i = 3       |      i = 6       |      i = 9       |
        # | grid[3] = (0, 0) | grid[6] = (1, 0) | grid[9] = (2, 0) |
        # +------------------+------------------+------------------+

        # Compute the bin dimensions
        bin_width = physical_size / n_bins
        bin_area  = bin_width * bin_width

        # Compute the x and y coordinates of each square bin
        shift   = 0.5 * (physical_size - bin_width)
        x_ticks = [(i - 1) * bin_width - shift + center[1] for i in 1:n_bins]
        y_ticks = [(i - 1) * bin_width - shift + center[2] for i in 1:n_bins]

        # Allocate memory
        grid = Matrix{NTuple{2,<:Number}}(undef, n_bins, n_bins)

        # Compute the position of each grid point
        for i in eachindex(grid)

            # The grid index `i` goes from top to bottom first, and then left to right,
            # starting at the top left of the grid
            i_x = ceil(Int, i / n_bins)
            i_y = mod1(i, n_bins) - 1

            # The coordinates are cartesian, so `y` goes from bottom to top and `x` goes from left to right,
            # starting at the bottom left of the grid
            grid[i] = (x_ticks[i_x], y_ticks[end - i_y])

        end

        new(grid, x_ticks, y_ticks, physical_size, n_bins, bin_width, bin_area)

    end
end

"""
Cubic grid (3D).

# Fields

  - `grid::Array{NTuple{3,<:Number},3}`: Matrix with the physical coordinates of each voxel in the grid.
  - `x_ticks::Vector{<:Number}`: Full set of possible values for the x coordinate.
  - `y_ticks::Vector{<:Number}`: Full set of possible values for the y coordinate.
  - `z_ticks::Vector{<:Number}`: Full set of possible values for the z coordinate.
  - `physical_size::Number`: Side length of the cubic grid.
  - `n_bins::Int`: Number of bins per side of the grid.
  - `bin_width::Number`: Side length of each bin.
  - `bin_area::Number`: Face area of each bin.
  - `bin_volume::Number`: Volume of each bin.
"""
struct CubicGrid
    grid::Array{NTuple{3,<:Number},3}
    x_ticks::Vector{<:Number}
    y_ticks::Vector{<:Number}
    z_ticks::Vector{<:Number}
    physical_size::Number
    n_bins::Int
    bin_width::Number
    bin_area::Number
    bin_volume::Number

    """
        CubicGrid(
            physical_size::Number,
            n_bins::Int;
            <keyword arguments>
        )

    Constructor for `CubicGrid`.

    # Arguments

      - `physical_size::Number`: Side length of the cubic grid.
      - `n_bins::Int`: Number of bins per dimesion of the grid.
      - `center::Vector{<:Number}=zeros(typeof(physical_size),3)`: 3D location of the center of the grid.
    """
    function CubicGrid(
        physical_size::Number,
        n_bins::Int;
        center::Vector{<:Number}=zeros(typeof(physical_size),3),
    )

        # Example of a 3x3x3 grid:
        #
        # x_ticks       = [0, 1, 2]
        # y_ticks       = [0, 1, 2]
        # z_ticks       = [0, 1, 2]
        # n_bins        = 3
        # physical_size = 3
        # center        = [1, 1, 1]
        #
        # Reference system:
        #
        #   y
        #   ↑
        #   |
        #   |
        #   | z
        #   +------→ x
        #
        # +----------------------+----------------------+----------------------+
        # |        i = 1         |        i = 4         |        i = 7         |
        # |  grid[1] = (0, 2, 0) |  grid[4] = (1, 2, 0) |  grid[7] = (2, 2, 0) |
        # +----------------------+----------------------+----------------------+
        # |        i = 2         |        i = 5         |        i = 8         |
        # |  grid[2] = (0, 1, 0) |  grid[5] = (1, 1, 0) |  grid[8] = (2, 1, 0) |
        # +----------------------+----------------------+----------------------+
        # |        i = 3         |        i = 6         |        i = 9         |
        # |  grid[3] = (0, 0, 0) |  grid[6] = (1, 0, 0) |  grid[9] = (2, 0, 0) |
        # +----------------------+----------------------+----------------------+
        #
        # +----------------------+----------------------+----------------------+
        # |        i = 10        |          i = 13      |          i = 16      |
        # | grid[10] = (0, 2, 1) | grid[13] = (1, 2, 1) | grid[16] = (2, 2, 1) |
        # +----------------------+----------------------+----------------------+
        # |        i = 11        |        i = 14        |        i = 17        |
        # | grid[11] = (0, 1, 1) | grid[14] = (1, 1, 1) | grid[17] = (2, 1, 1) |
        # +----------------------+----------------------+----------------------+
        # |        i = 12        |        i = 15        |        i = 18        |
        # | grid[12] = (0, 0, 1) | grid[15] = (1, 0, 1) | grid[18] = (2, 0, 1) |
        # +----------------------+----------------------+----------------------+
        #
        # +----------------------+----------------------+----------------------+
        # |        i = 19        |        i = 22        |        i = 25        |
        # | grid[19] = (0, 2, 2) | grid[22] = (1, 2, 2) | grid[25] = (2, 2, 2) |
        # +----------------------+----------------------+----------------------+
        # |        i = 22        |        i = 23        |        i = 26        |
        # | grid[20] = (0, 1, 2) | grid[23] = (1, 1, 2) | grid[26] = (2, 1, 2) |
        # +----------------------+----------------------+----------------------+
        # |        i = 21        |        i = 24        |        i = 27        |
        # | grid[21] = (0, 0, 2) | grid[24] = (1, 0, 2) | grid[27] = (2, 0, 2) |
        # +----------------------+----------------------+----------------------+

        # Compute the bin dimensions
        bin_width  = physical_size / n_bins
        bin_area   = bin_width * bin_width
        bin_volume = bin_area * bin_width

        # Compute the x, y, and z coordinates of each cubic bin
        shift   = 0.5 * (physical_size - bin_width)
        x_ticks = [(i - 1) * bin_width - shift + center[1] for i in 1:n_bins]
        y_ticks = [(i - 1) * bin_width - shift + center[2] for i in 1:n_bins]
        z_ticks = [(i - 1) * bin_width - shift + center[3] for i in 1:n_bins]

        # Allocate memory
        grid = Array{NTuple{3,<:Number},3}(undef, n_bins, n_bins, n_bins)

        # Compute the position of each grid point
        for i in eachindex(grid)

            flat_idx = mod1(i, n_bins * n_bins)

            i_x = ceil(Int, flat_idx / n_bins)
            i_y = mod1(flat_idx, n_bins) - 1
            i_z = ceil(Int, i / (n_bins * n_bins))

            grid[i] = (x_ticks[i_x], y_ticks[end - i_y], z_ticks[i_z])

        end

        new(grid, x_ticks, y_ticks, z_ticks, physical_size, n_bins, bin_width, bin_area, bin_volume)

    end
end

"""
Circular grid (2D or 3D).

Series of concentric rings or spherical shells.

# Fields

  - `grid::Vector{<:Number}`: Vector with the distance of each bin to the center of the grid.
  - `ticks::Vector{<:Number}`: Vector with the edges of the bins.
  - `center::Vector{<:Number}`: 3D location of the center of the grid. In the 2D case the grid is assumed to be in the xy plane.
  - `bin_area::Vector{<:Number}`: Area of each ring.
  - `bin_volumes::Vector{<:Number}`: Volume of each spherical shell.
  - `log::Bool`: If the grid is logarithmic.
"""
struct CircularGrid
    grid::Vector{<:Number}
    ticks::Vector{<:Number}
    center::Vector{<:Number}
    bin_areas::Vector{<:Number}
    bin_volumes::Vector{<:Number}
    log::Bool

    """
        CircularGrid(
            radius::Number,
            n_bins::Int;
            <keyword arguments>
        )

    Constructor for `CircularGrid`.

    # Arguments

      - `radius::Number`: Radius of the grid (equal to the last bin tick).
      - `n_bins::Int`: Number of bins.
      - `center::Vector{<:Number}=zeros(typeof(radius), 3)`: 3D location of the center of the grid. In the 2D case the grid is assumed to be in the xy plane.
      - `log::Bool=false`: If the bins will be logarithmic.
      - `shift::Number=zero(radius)`: Distance of the first bin tick to the center.
    """
    function CircularGrid(
        radius::Number,
        n_bins::Int;
        center::Vector{<:Number}=zeros(typeof(radius), 3),
        log::Bool=false,
        shift::Number=zero(radius),
    )

        # Example of a grid with 3 bins:
        #
        # |<-----------------------------------------radius---------------------------------------->|
        #
        # +-----------+-------------------------+-------------------------+-------------------------+
        # |<--shift-->| grid[1] -->|            | grid[2] -->|            | grid[3] -->|            |
        # +-----------+-------------------------+-------------------------+-------------------------+
        #
        # +-----------+-------------------------+-------------------------+-------------------------+
        # |<--shift-->|<-- tick[1]   tick[2] -->|              tick[3] -->|              tick[4] -->|
        # +-----------+-------------------------+-------------------------+-------------------------+

        (
            isPositive(radius) ||
            throw(ArgumentError("CircularGrid: `radius` must be strictly positive, \
            but I got `radius` = $(radius)"))
        )

        if log

            (
                isPositive(shift) ||
                throw(ArgumentError("CircularGrid: For a logarithmic grid you need a \
                strictly positive `shift`, but I got `shift` = $(shift)"))
            )

            # Length unit
            u_l = unit(radius)

            log_shift  = log10(ustrip(u_l, shift))
            log_radius = log10(ustrip(radius))

            width = (log_radius - log_shift) / n_bins
            grid  = [exp10((i - 0.5) * width + log_shift) * u_l for i in 1:n_bins]
            ticks = [exp10(i * width + log_shift) * u_l for i in 0:n_bins]

        else

            width = (radius - shift) / n_bins
            grid  = [(i - 0.5) * width + shift for i in 1:n_bins]
            ticks = [i * width + shift for i in 0:n_bins]

        end

        bin_areas   = [area(ticks[i + 1]) - area(ticks[i]) for i in 1:n_bins]
        bin_volumes = [volume(ticks[i + 1]) - volume(ticks[i]) for i in 1:n_bins]

        new(grid, ticks, center, bin_areas, bin_volumes, log)

    end
end

"""
Plotting parameters for a quantity.

# Fields

  - `request::Dict{Symbol,Vector{String}} = Dict{Symbol,Vector{String}}()`: Data request for [`readSnapshot`](@ref). It must have the shape `cell/particle type` -> [`block`, `block`, `block`, ...].
  - `var_name::AbstractString = ""`: Name of the quantity for the plot axis. It should not include units or scaling factors.
  - `exp_factor::Int = 0`: Numerical exponent to scale down the axis, e.g. if `x_exp_factor` = 10 the values will be divided by ``10^{10}``. The default is no scaling.
  - `unit::Unitful.Units = Unitful.NoUnits`: Target unit for the axis.
  - `axis_label::AbstractString = "auto_label"`: Label for the axis. It can contain the string `auto_label`, which will be replaced by the default label: `var_name` / 10^`exp_factor` `unit`.
"""
@kwdef mutable struct PlotParams
    request::Dict{Symbol,Vector{String}} = Dict{Symbol,Vector{String}}()
    var_name::AbstractString = ""
    exp_factor::Int = 0
    unit::Unitful.Units = Unitful.NoUnits
    axis_label::AbstractString = "auto_label"
end

##########################
# Code specific constants
##########################

include("./arepo.jl")

######################
# Cell/particle types
######################

"""
Type of cell/particle corresponding to each code index.
"""
const INDEX_PARTICLE = Dict(n => symbol for (symbol, n) in PARTICLE_INDEX)

"""
Type of cell/particle corresponding to each internal code name (data group in the HDF5 output).
"""
const PARTICLE_TYPE = Dict("PartType$n" => symbol for (symbol, n) in PARTICLE_INDEX)

"""
Internal code name (data group in the HDF5 output) corresponding to each type of cell/particle.
"""
const PARTICLE_CODE_NAME = Dict(symbol => "PartType$n" for (symbol, n) in PARTICLE_INDEX)

"""
Default filter dictionary that does not exclude any cells/particles.
"""
const PASS_ALL = Dict(key => (:) for key in keys(PARTICLE_INDEX))

"""
Filter that excludes every cell/particle.
"""
const PASS_NONE = Dict(key => Int[] for key in keys(PARTICLE_INDEX))

###################
# Tracked elements
###################

"""
Symbol list for the gas abundance quantities.
"""
const GAS_ABUNDANCE = [Symbol(element, "_gas_abundance") for element in keys(ELEMENT_INDEX)]

"""
Symbol list for the stellar abundance quantities.
"""
const STELLAR_ABUNDANCE = [Symbol(element, "_stellar_abundance") for element in keys(ELEMENT_INDEX)]
