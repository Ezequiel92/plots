push!(LOAD_PATH, "./src/")
using Documenter, GalaxyInspector

# True if it is being deployed, false if it is being compile locally
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

# Compile documentation
makedocs(
    sitename="GalaxyInspector.jl",
    authors="Ezequiel Lozano",
    format   = Documenter.HTML(
        prettyurls=CI,
        assets=["./assets/GI-docs.css"],
        warn_outdated=true,
        collapselevel=1,
        edit_link="main",
        size_threshold=500000,
        size_threshold_warn=400000,
    ),
    modules=[GalaxyInspector],
    pages=[
        "Introduction" => "intro.md",
        "API" => Any[
            "api/plotting/pipelines.md",
            "api/plotting/convenience.md",
            "api/plotting/post_processing.md",
            "api/constants/globals.md",
            "api/constants/arepo.md",
            "api/analysis/data_analysis.md",
            "api/analysis/data_acquisition.md",
            "api/analysis/compute_quantities/positions.md",
            "api/analysis/compute_quantities/velocities.md",
            "api/analysis/compute_quantities/masses.md",
            "api/analysis/compute_quantities/other.md",
            "api/auxiliary_functions.md",
            "api/analysis/filters.md",
            "api/analysis/tracers.md",
            "api/analysis/transformations.md",
        ],
        "Index" => "index.md",
    ],
    warnonly=[:missing_docs],
)

# Deploy the documentation to GitHub pages
if CI
    deploydocs(
        repo="github.com/Ezequiel92/GalaxyInspector.git",
        devbranch="main",
        versions=["dev" => "dev"],
    )
end
