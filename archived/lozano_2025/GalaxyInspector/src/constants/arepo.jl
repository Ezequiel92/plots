####################################################################################################
# Constants and data structures
####################################################################################################

#################
# Code constants
#################

"""
Base name of the snapshot files, set in the code variable `SnapshotFileBase`.
"""
const SNAP_BASENAME = "snap"

"""
Base name of the group catalog files.
"""
const GC_BASENAME = "fof_subhalo_tab"

"""
Mass fraction of hydrogen.
"""
const HYDROGEN_MASSFRAC = 0.76

"""
Mass of the tracers in internal code units. Its value comes from `All.TargetGasMass = All.TargetGasMassFactor * All.ReferenceGasPartMass` in the code.
"""
const TRACER_MASS = 3.65456e-06

"""
Solar metallicity, as used in Arepo.

# References

M. Asplund et al. (2006). *The new solar abundances - Part I: the observations*. Communications in Asteroseismology, **147**. [doi:10.1553/cia147s76](https://doi.org/10.1553/cia147s76)
"""
const SOLAR_METALLICITY = 0.0127

######################
# Cell/particle types
######################

"""
Current cell/particle index in use.
"""
const PARTICLE_INDEX = LONG_PARTICLE_INDEX

"""
Current human readable name of each cell/particle type in use.
"""
const PARTICLE_NAMES = LONG_PARTICLE_NAMES

###################
# Tracked elements
###################

"""
Code index for each tracked element.
"""
const ELEMENT_INDEX = Dict(
    :H     => 1,  # Hydrogen
    :He    => 2,  # Helium
    :C     => 3,  # Carbon
    :N     => 4,  # Nitrogen
    :O     => 5,  # Oxygen
    :Ne    => 6,  # Neon
    :Mg    => 7,  # Magnesium
    :Si    => 8,  # Silicon
    :Fe    => 9,  # Iron
    :Other => 10, # All other
)

"""
List of element indices above helium.
"""
const METAL_LIST = [3, 4, 5, 6, 7, 8, 9, 10]

#############################################################
# Quantities that can be in a snapshot or group catalog file
#############################################################

"""
Dictionary of dimensional properties for the quantities in the code.
"""
const QUANTITIES = Dict(

    ######################
    # Snapshot quantities
    ######################

    "ACCE" => Qty("Acceleration", Unitful.𝐋 * Unitful.𝐓^-2, :internal),
    "AGE " => Qty("StellarFormationTime", Unitful.𝐓, :internal),
    "AREA" => Qty("SurfaceArea", Unitful.𝐋^2, :internal),
    "CLKT" => Qty("", Unitful.𝐓, :internal),
    "CMCE" => Qty("CenterOfMass", Unitful.𝐋, :internal),
    "CSND" => Qty("SoundSpeed", Unitful.𝐋 * Unitful.𝐓^-1, :internal),
    "GAGE" => Qty("GFM_StellarFormationTime", Unitful.𝐓, :internal),
    "GIMA" => Qty("GFM_InitialMass", Unitful.𝐌, :internal),
    "GME2" => Qty("GFM_Metals", Unitful.NoDims, Unitful.NoUnits),
    "GMET" => Qty("GFM_Metals", Unitful.NoDims, Unitful.NoUnits),
    "GZ  " => Qty("GFM_Metallicity", Unitful.NoDims, Unitful.NoUnits),
    "GZ2 " => Qty("GFM_Metallicity", Unitful.NoDims, Unitful.NoUnits),
    "Hepp" => Qty("DoublyIonizedHeliumAbundance", Unitful.NoDims, Unitful.NoUnits),
    "ID  " => Qty("ParticleIDs", Unitful.NoDims, Unitful.NoUnits),
    "PAID" => Qty("ParentID", Unitful.NoDims, Unitful.NoUnits),
    "TRID" => Qty("TracerID", Unitful.NoDims, Unitful.NoUnits),
    "MACH" => Qty("MachNumber", Unitful.NoDims, Unitful.NoUnits),
    "MASS" => Qty("Masses", Unitful.𝐌, :internal),
    "NE  " => Qty("ElectronAbundance", Unitful.NoDims, Unitful.NoUnits),
    "NFAC" => Qty("NumFacesCell", Unitful.NoDims, Unitful.NoUnits),
    "NH  " => Qty("NeutralHydrogenAbundance", Unitful.NoDims, Unitful.NoUnits),
    "NHP " => Qty("IonizedHydrogenAbundance", Unitful.NoDims, Unitful.NoUnits),
    "NHeO" => Qty("NeutralHeliumAbundance", Unitful.NoDims, Unitful.NoUnits),
    "NHep" => Qty("SinglyIonizedHeliumAbundance", Unitful.NoDims, Unitful.NoUnits),
    "POS " => Qty("Coordinates", Unitful.𝐋, :internal),
    "PRES" => Qty("Pressure", Unitful.𝐌 * Unitful.𝐋^-1 * Unitful.𝐓^-2, :internal),
    "RHO " => Qty("Density", Unitful.𝐌 * Unitful.𝐋^-3, :internal),
    "SFR " => Qty("StarFormationRate", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "SOFT" => Qty("Softenings", Unitful.𝐋, :internal),
    "TBH"  => Qty("TimebinHydro", Unitful.NoDims, Unitful.NoUnits),
    "TEMP" => Qty("Temperature", Unitful.𝚯, u"K"),
    "TSTP" => Qty("TimeStep", Unitful.𝐓, :internal),
    "U   " => Qty("InternalEnergy", Unitful.𝐋^2 * Unitful.𝐓^-2, :internal),
    "VEL " => Qty("Velocities", Unitful.𝐋 * Unitful.𝐓^-1, :internal),
    "VEVE" => Qty("VertexVelocity", Unitful.𝐋 * Unitful.𝐓^-1, :internal),
    "VOL " => Qty("Volume", Unitful.𝐋^3, :internal),
    "Z   " => Qty("Metallicity", Unitful.NoDims, Unitful.NoUnits),

    #####################
    # sfr.txt quantities
    #####################

    # Time or scale factor
    "SFC1" => Qty("", Unitful.𝐓, :internal),
    # Total stellar mass to be formed prior to stochastic sampling
    "SFC2" => Qty("", Unitful.𝐌, :internal),
    # Instantaneous star formation rate of all cells
    "SFC3" => Qty("", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    # Instantaneous star formation rate of active cells
    "SFC4" => Qty("", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    # Total mass in stars formed after stochastic sampling
    "SFC5" => Qty("", Unitful.𝐌, :internal),
    # Cumulative stellar mass formed
    "SFC6" => Qty("", Unitful.𝐌, :internal),

    ####################
    # EL_SFR quantities
    ####################

    # Integration time, in gas cells and stellar particles
    "ODIT" => Qty("ODE_IntegrationTime", Unitful.𝐓, u"Myr"),
    # Scale factor, in gas cells and stellar particles
    "PARA" => Qty("ODE_ParameterA", Unitful.NoDims, Unitful.NoUnits),
    # Star formation time parameter, in gas cells and stellar particles
    "TAUS" => Qty("ODE_TauS", Unitful.𝐓, u"Myr"),
    # Gas density, for gas cells and stellar particles
    "RHOC" => Qty("ODE_ParameterCellDensity", Unitful.𝐋^-3, u"cm^-3"),
    # Gas metallicity, for gas cells and stellar particles
    "PARZ" => Qty("ODE_ParameterMetallicity", Unitful.NoDims, Unitful.NoUnits),
    # Photodissociation parameter, for gas cells and stellar particles
    "ETAD" => Qty("ODE_ParameterEtaD", Unitful.NoDims, Unitful.NoUnits),
    # Photoionization parameter, for gas cells and stellar particles
    "ETAI" => Qty("ODE_ParameterEtaI", Unitful.NoDims, Unitful.NoUnits),
    # Mass recycling parameter, for gas cells and stellar particles
    "PARR" => Qty("ODE_ParameterR", Unitful.NoDims, Unitful.NoUnits),
    # Gas fractions, for gas cells and stellar particles
    "FRAC" => Qty("ODE_Fractions", Unitful.NoDims, Unitful.NoUnits),
    # Star formation flag, for gas cells
    "SFFL" => Qty("ODE_SfFlag", Unitful.NoDims, Unitful.NoUnits),
    # Cold gas fraction, for gas cells
    "COLF" => Qty("ODE_ColdMassFrac", Unitful.NoDims, Unitful.NoUnits),
    # Parent gas mass (at the moment of star formation), for stellar particles
    "GMAS" => Qty("ODE_GasMass", Unitful.𝐌, :internal),
    # Parent SFR (at the moment of star formation), for stellar particles
    "GSFR" => Qty("ODE_GasSFR", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    # Parent gas pressure (at the moment of star formation), for stellar particles
    "GPRE" => Qty("ODE_GasPressure", Unitful.𝐌 * Unitful.𝐋^-1 * Unitful.𝐓^-2, :internal),
    # Parent position (at the moment of star formation), for stellar particles
    "GPOS" => Qty("ODE_GasPosition", Unitful.𝐋, :internal),
    # Parent velocity (at the moment of star formation), for stellar particles
    "GVEL" => Qty("ODE_GasVelocity", Unitful.𝐋 * Unitful.𝐓^-1, :internal),

    # Legacy quantities
    "ACIT" => Qty("ODE_AccumulatedIntegrationTime", Unitful.𝐓, u"Myr"),
    "CTIM" => Qty("ODE_CurrentTime", Unitful.𝐓, u"Myr"),
    "DTIM" => Qty("ODE_DeltaTime", Unitful.𝐓, u"Myr"),
    "FSPR" => Qty("ODE_f_star_prev", Unitful.NoDims, Unitful.NoUnits),

    ##############################
    # Halo (FoF group) quantities
    ##############################

    "G_BHMass"             => Qty("GroupBHMass", Unitful.𝐌, :internal),
    "G_CM"                 => Qty("GroupCM", Unitful.𝐋, :internal),
    "G_FirstSub"           => Qty("GroupFirstSub", Unitful.NoDims, Unitful.NoUnits),
    "G_GasMetalFractions"  => Qty("GroupGasMetalFractions", Unitful.NoDims, Unitful.NoUnits),
    "G_GasMetallicity"     => Qty("GroupGasMetallicity", Unitful.NoDims, Unitful.NoUnits),
    "G_Len"                => Qty("GroupLen", Unitful.NoDims, Unitful.NoUnits),
    "G_LenType"            => Qty("GroupLenType", Unitful.NoDims, Unitful.NoUnits),
    "G_Mass"               => Qty("GroupMass", Unitful.𝐌, :internal),
    "G_MassType"           => Qty("GroupMassType", Unitful.𝐌, :internal),
    "G_Nsubs"              => Qty("GroupNsubs", Unitful.NoDims, Unitful.NoUnits),
    "G_Pos"                => Qty("GroupPos", Unitful.𝐋, :internal),
    "G_SFR"                => Qty("GroupSFR", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "G_StarMetalFractions" => Qty("GroupStarMetalFractions", Unitful.NoDims, Unitful.NoUnits),
    "G_StarMetallicity"    => Qty("GroupStarMetallicity", Unitful.NoDims, Unitful.NoUnits),
    "G_WindMass"           => Qty("GroupWindMass", Unitful.𝐌, :internal),
    "G_M_Crit200"          => Qty("Group_M_Crit200", Unitful.𝐌, :internal),
    "G_M_Crit500"          => Qty("Group_M_Crit500", Unitful.𝐌, :internal),
    "G_M_Mean200"          => Qty("Group_M_Mean200", Unitful.𝐌, :internal),
    "G_M_TopHat200"        => Qty("Group_M_TopHat200", Unitful.𝐌, :internal),
    "G_R_Crit200"          => Qty("Group_R_Crit200", Unitful.𝐋, :internal),
    "G_R_Crit500"          => Qty("Group_R_Crit500", Unitful.𝐋, :internal),
    "G_R_Mean200"          => Qty("Group_R_Mean200", Unitful.𝐋, :internal),
    "G_R_TopHat200"        => Qty("Group_R_TopHat200", Unitful.𝐋, :internal),
    "G_Vel"                => Qty("GroupVel", Unitful.𝐋 * Unitful.𝐓^-1, :gvel),

    ###############################
    # Subhalo (subfind) quantities
    ###############################

    "S_BHMass"                       => Qty("SubhaloBHMass", Unitful.𝐌, :internal),
    "S_CM"                           => Qty("SubhaloCM", Unitful.𝐋, :internal),
    "S_GasMetalFractions"            => Qty("SubhaloGasMetalFractions", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetalFractionsHalfRad"     => Qty("SubhaloGasMetalFractionsHalfRad", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetalFractionsMaxRad"      => Qty("SubhaloGasMetalFractionsMaxRad", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetalFractionsSfr"         => Qty("SubhaloGasMetalFractionsSfr", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetalFractionsSfrWeighted" => Qty("SubhaloGasMetalFractionsSfrWeighted", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetallicity"               => Qty("SubhaloGasMetallicity", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetallicityHalfRad"        => Qty("SubhaloGasMetallicityHalfRad", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetallicityMaxRad"         => Qty("SubhaloGasMetallicityMaxRad", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetallicitySfr"            => Qty("SubhaloGasMetallicitySfr", Unitful.NoDims, Unitful.NoUnits),
    "S_GasMetallicitySfrWeighted"    => Qty("SubhaloGasMetallicitySfrWeighted", Unitful.NoDims, Unitful.NoUnits),
    "S_GrNr"                         => Qty("SubhaloGrNr", Unitful.NoDims, Unitful.NoUnits),
    "S_HalfmassRad"                  => Qty("SubhaloHalfmassRad", Unitful.𝐋, :internal),
    "S_HalfmassRadType"              => Qty("SubhaloHalfmassRadType", Unitful.𝐋, :internal),
    "S_IDMostbound"                  => Qty("SubhaloIDMostbound", Unitful.NoDims, Unitful.NoUnits),
    "S_Len"                          => Qty("SubhaloLen", Unitful.NoDims, Unitful.NoUnits),
    "S_LenType"                      => Qty("SubhaloLenType", Unitful.NoDims, Unitful.NoUnits),
    "S_Mass"                         => Qty("SubhaloMass", Unitful.𝐌, :internal),
    "S_MassInHalfRad"                => Qty("SubhaloMassInHalfRad", Unitful.𝐌, :internal),
    "S_MassInHalfRadType"            => Qty("SubhaloMassInHalfRadType", Unitful.𝐌, :internal),
    "S_MassInMaxRad"                 => Qty("SubhaloMassInMaxRad", Unitful.𝐌, :internal),
    "S_MassInMaxRadType"             => Qty("SubhaloMassInMaxRadType", Unitful.𝐌, :internal),
    "S_MassInRad"                    => Qty("SubhaloMassInRad", Unitful.𝐌, :internal),
    "S_MassInRadType"                => Qty("SubhaloMassInRadType", Unitful.𝐌, :internal),
    "S_MassType"                     => Qty("SubhaloMassType", Unitful.𝐌, :internal),
    "S_Parent"                       => Qty("SubhaloParent", Unitful.NoDims, Unitful.NoUnits),
    "S_Pos"                          => Qty("SubhaloPos", Unitful.𝐋, :internal),
    "S_SFR"                          => Qty("SubhaloSFR", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "S_SFRinHalfRad"                 => Qty("SubhaloSFRinHalfRad", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "S_SFRinMaxRad"                  => Qty("SubhaloSFRinMaxRad", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "S_SFRinRad"                     => Qty("SubhaloSFRinRad", Unitful.𝐌 * Unitful.𝐓^-1, u"Msun*yr^-1"),
    "S_StarMetalFractions"           => Qty("SubhaloStarMetalFractions", Unitful.NoDims, Unitful.NoUnits),
    "S_StarMetalFractionsHalfRad"    => Qty("SubhaloStarMetalFractionsHalfRad", Unitful.NoDims, Unitful.NoUnits),
    "S_StarMetalFractionsMaxRad"     => Qty("SubhaloStarMetalFractionsMaxRad", Unitful.NoDims, Unitful.NoUnits),
    "S_StarMetallicity"              => Qty("SubhaloStarMetallicity", Unitful.NoDims, Unitful.NoUnits),
    "S_StarMetallicityHalfRad"       => Qty("SubhaloStarMetallicityHalfRad", Unitful.NoDims, Unitful.NoUnits),
    "S_StarMetallicityMaxRad"        => Qty("SubhaloStarMetallicityMaxRad", Unitful.NoDims, Unitful.NoUnits),
    "S_StellarPhotometrics"          => Qty("SubhaloStellarPhotometrics", Unitful.NoDims, Unitful.NoUnits),
    "S_StellarPhotometricsMassInRad" => Qty("SubhaloStellarPhotometricsMassInRad", Unitful.𝐌, :internal),
    "S_StellarPhotometricsRad"       => Qty("SubhaloStellarPhotometricsRad", Unitful.𝐋, :internal),
    "S_Vel"                          => Qty("SubhaloVel", Unitful.𝐋 * Unitful.𝐓^-1, u"km*s^-1"),
    "S_VelDisp"                      => Qty("SubhaloVelDisp", Unitful.𝐋 * Unitful.𝐓^-1, u"km*s^-1"),
    "S_Vmax"                         => Qty("SubhaloVmax", Unitful.𝐋 * Unitful.𝐓^-1, u"km*s^-1"),
    "S_VmaxRad"                      => Qty("SubhaloVmaxRad", Unitful.𝐋, :internal),
    "S_WindMass"                     => Qty("SubhaloWindMass", Unitful.𝐌, :internal),
)
