####################################################################################################
#
#    ______        __                     ____                                 __                  _  __
#   / ____/____ _ / /____ _ _  __ __  __ /  _/____   _____ ____   ___   _____ / /_ ____   _____   (_)/ /
#  / / __ / __ `// // __ `/| |/_// / / / / / / __ \ / ___// __ \ / _ \ / ___// __// __ \ / ___/  / // /
# / /_/ // /_/ // // /_/ /_>  < / /_/ /_/ / / / / /(__  )/ /_/ //  __// /__ / /_ / /_/ // /_    / // /
# \____/ \__,_//_/ \__,_//_/|_| \__, //___//_/ /_//____// .___/ \___/ \___/ \__/ \____//_/(_)__/ //_/
#                              /____/                  /_/                                  /___/
####################################################################################################

####################################################################################################
# A Julia module for the analysis of galaxy simulations
####################################################################################################

module GalaxyInspector

using CSV,
    CairoMakie,
    ColorSchemes,
    Colors,
    DataFrames,
    DelimitedFiles,
    FileIO,
    GLM,
    Glob,
    HDF5,
    Images,
    JLD2,
    LaTeXStrings,
    LinearAlgebra,
    Logging,
    Measurements,
    NearestNeighbors,
    ProgressMeter,
    QuadGK,
    Rotations,
    Statistics,
    StatsBase,
    Unitful,
    UnitfulAstro

####################################################################################################
# Optimization
####################################################################################################

@eval Base.Experimental.@optlevel 3

####################################################################################################
# Submodules
####################################################################################################

include("constants/globals.jl")

include("auxiliary_functions.jl")

include("analysis/data_acquisition.jl")
include("analysis/compute_quantities/positions.jl")
include("analysis/compute_quantities/velocities.jl")
include("analysis/compute_quantities/masses.jl")
include("analysis/compute_quantities/other.jl")
include("analysis/filters.jl")
include("analysis/tracers.jl")
include("analysis/transformations.jl")
include("analysis/data_analysis.jl")

include("plotting/post_processing.jl")
include("plotting/pipelines.jl")
include("plotting/convenience.jl")

####################################################################################################
# Public functions
####################################################################################################

# From `analysis/data_acquisition.jl`
export readGroupCatalog
export readSnapshot
export getBlock
export makeDataDict

# From `plotting/pipelines.jl`
export plotSnapshot
export plotTimeSeries

# From `plotting/convenience.jl`
export snapshotReport
export simulationReport
export sfrTXT
export cpuTXT
export stellarBirthHalos
export densityMap
export gasSFRMap
export densityMapVelField
export metallicityMap
export temperatureMap
export scatterPlot
export scatterDensityMap
export atomicMolecularTransition
export gasBarPlot
export timeSeries
export gasEvolution
export virialAccretionEvolution
export discAccretionEvolution
export rotationCurve
export densityProfile
export massProfile
export velocityProfile
export stellarHistory
export lineHistogram
export compareFeldmann2020
export compareMolla2015
export kennicuttSchmidtLaw
export fitVSFLaw
export massMetallicityRelation
export gasVelocityCubes
export stellarVelocityCubes
export clumpingFactor

end
