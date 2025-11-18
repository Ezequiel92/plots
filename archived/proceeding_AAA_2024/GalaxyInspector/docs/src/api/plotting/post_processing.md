# Post processing functions

These functions are used internally and none are exported. 

A post-processing function must take a [Makie](https://docs.makie.org/stable/) figure, add something to it, and return how to label the additions (or `nothing` when no new labels should be drawn).

### Expected signature:

```julia
  post_processing(figure, args...; kwargs...) -> ([marker, ...], [label, ...])
```

where:

  - `figure::Makie.Figure`
  - `marker::LegendElement`
  - `label::String`

---

```@autodocs
Modules = [GalaxyInspector]
Pages   = ["plotting/post_processing.jl"]
```