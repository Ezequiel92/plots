####################################################################################################
# Auxiliary functions
####################################################################################################

"""
    metaFormatter(level::LogLevel, _module, group, id, file, line)

Formatter for loggers.

See the documentation for [ConsoleLogger](https://docs.julialang.org/en/v1/stdlib/Logging/#Base.CoreLogging.ConsoleLogger)
"""
function metaFormatter(level::LogLevel, _module, group, id, file, line)

    @nospecialize

    color = Logging.default_logcolor(level)
    prefix = string(level == Warn ? "Warning" : string(level), " |")
    suffix::String = ""

    if file !== nothing
        suffix *= contractuser(file)::String
        if line !== nothing
            suffix *= ":$(isa(line, UnitRange) ? "$(first(line))-$(last(line))" : line)"
        end
    end

    if !isempty(suffix)
        suffix = "@ " * suffix * "\n"
    else
        suffix *= "\n"
    end

    return color, prefix, suffix

end

"""
    setLogging!(verb::Bool; <keyword arguments>)::Nothing

Set if logging messages will printed out. By default no logs are printed.

# Arguments

  - `verb::Bool`: If logs will be printed out using the default logger.
  - `stream::IO=stdout`: Sets where to print the logs. It can be a file.
"""
function setLogging!(verb::Bool; stream::IO=stdout)::Nothing

    logging[] = verb

    verb && global_logger(ConsoleLogger(stream; meta_formatter=metaFormatter))

    return nothing

end

"""
    ring(vec::Vector, index::Integer)::Vector

Make the indexing operation `vec[index]` work using modular arithmetic for the indices.

# Arguments

  - `vec::Vector`: Vector.
  - `index::Integer`: Index.

# Returns

  - `vec[mod1(index, length(vec))]`

# Examples

```julia-repl
julia> ring([1, 2, 3], 11)
2

julia> ring([1, 2, 3], 3)
3

julia> ring([1, 2, 3], -5)
1
```
"""
ring(vec::Vector, index::Integer) = vec[mod1(index, length(vec))]

"""
    parserWS(data::AbstractString)::Union{Float64,Missing}

Parse a string as a Float64, ignoring white spaces. If the string is empty return missing.

# Arguments

  - `data::AbstractString`: String to be parsed.

# Returns

  - Number in the string as a Float64.
"""
function parserWS(data::AbstractString)::Union{Float64,Missing}

	clean_data = strip(data)

	!isempty(clean_data) || return missing

	return parse(Float64, clean_data)

end

"""
    safeSelect(vec::Vector, index::IndexType)

Make the indexing operation `vec[index]` ignoring indices that are out of bounds.

# Arguments

  - `vec::Vector`: Vector.
  - `index::IndexType`: Indices. It can be an integer (a single element), a vector of integers (several elements), an `UnitRange` (e.g. 5:13), an `StepRange` (e.g. 5:2:13) or (:) (every element).

# Returns

  - `vec[index (minus out of bounds indices)]`

# Examples

```julia-repl
julia> safeSelect([1, 2, 3], 11)
Int[]

julia> safeSelect([1, 2, 3], 1:5)
3-element Vector{Int}:
 1
 2
 3

julia> safeSelect([1, 2, 3], 1:3:10)
1

julia> safeSelect([1, 2, 3], [1, 2, 5, 9])
2-element Vector{Int}:
 1
 2

julia> safeSelect([1, 2, 3], (:))
3-element Vector{Int}:
 1
 2
 3
```
"""
function safeSelect(vec::Vector, index::IndexType)

    index != (:) || return vec

    index_list = [index...]

    filter!(x -> x <= length(vec), index_list)

    (
        length(index_list) == length([index...]) || !logging[] ||
        @info("safeSelect: There are out of bounds indices")
    )

    if length(index_list) == 1
        return vec[index_list...]
    else
        return vec[index_list]
    end

end

"""
Create a copy of `list` with every negative value set to 0.
"""
setPositive(list::VecOrMat{<:Number}) = replace(x -> x >= zero(x) ? x : zero(x), list)
setPositive(x::Number) = x >= zero(x) ? x : zero(x)

"""
Test for strict positivity.
"""
isPositive(x::Number)::Bool = x > zero(x)
isPositive(x::AbstractArray)::Bool = all(isPositive, x)
isPositive(x...)::Bool = all(isPositive, x)

"""
New method for `Base.iszero` to compare [`IndexType`](@ref) with 0 as an integer.
"""
Base.iszero(x::IndexType)::Bool = x == 0

"""
New method for `Base.isempty` to check for empty [LaTeXStrings](https://github.com/JuliaStrings/LaTeXStrings.jl).
"""
Base.isempty(l_str::LaTeXString)::Bool = l_str == L""

"""
New methods for `Base.intersect` to use with the `Colon` type.
"""
Base.intersect(a1::Colon, a2::IndexType)::IndexType = a2
Base.intersect(a1::IndexType, a2::Colon)::IndexType = a1
Base.intersect(a1::Colon, a2::Colon)::Colon = (:)

"""
New methods for `Base.intersect` to use with the `Vector{Bool}` type.
"""
Base.intersect(a1::Vector{Bool}, a2::ReducedIndexType)::Vector{Int} = findall(a1) ∩ a2
Base.intersect(a1::ReducedIndexType, a2::Vector{Bool})::Vector{Int} = a1 ∩ findall(a2)
Base.intersect(a1::Vector{Bool}, a2::Vector{Bool})::Vector{Bool} = Vector{Bool}(a1 .&& a2)

"""
New methods for `Base.union` to use with the `Vector{Bool}` type.
"""
Base.union(a1::Vector{Bool}, a2::ReducedIndexType)::Vector{Int} = findall(a1) ∪ a2
Base.union(a1::ReducedIndexType, a2::Vector{Bool})::Vector{Int} = a1 ∪ findall(a2)
Base.union(a1::Vector{Bool}, a2::Vector{Bool})::Vector{Bool} = Vector{Bool}(a1 .|| a2)

"""
Area of a circle with radius `r`.
"""
function area(r::Number)::Number

    x = setPositive(r)

    return π * x * x

end

"""
Volume of a sphere with radius `r`.
"""
function volume(r::Number)::Number

    x = setPositive(r)

    return π * x * x * x * 1.333

end

"""
    evaluateNormal(data::Vector{<:Number})::Vector{<:Number}

Evaluate a normal distribution at the values in `data`.

The mean and standard deviation of the distribution are the ones from `data` itself.

# Arguments

  - `data::Vector{<:Number}`: Data vector used to compute the mean and standard deviation of the normal distribution.

# Returns

  - The normal distribution evaluated at the values in `data`.
"""
function evaluateNormal(data::Vector{<:Number})::Vector{<:Number}

    μ  = mean(data)
    σ2 = 2.0 * std(data; mean=μ)^2
    d2 = @. (data - μ)^2

    return @. (1.0 / sqrt(σ2 * π)) * exp(-d2 / σ2)

end

"""
Always returns `nothing`, for any type and number of arguments.
"""
getNothing(x...; y...)::Nothing = nothing

"""
Always returns an empty vector, for any type and number of arguments.
"""
getEmpty(x...; y...)::Vector = []

"""
    hvcatImages(
        blocks_per_row::Int,
        paths::Vector{String};
        <keyword arguments>
    )::Nothing

Join several images vertically and horizontally.

The elements will fill the rows and columns starting at the top left, going from left to right and from top to bottom (row-major order).

# Arguments

  - `blocks_per_row::Int`: Number of columns.
  - `paths::Vector{String}`: Paths to the images.
  - `output_path::String="./joined_images.png"`: Path to the output image.
"""
function hvcatImages(
    blocks_per_row::Int,
    paths::Vector{String};
    output_path::String="./joined_images.png",
)::Nothing

    isempty(paths) && throw(ArgumentError("hvcatImages: `paths` is empty"))

    new_image = hvcat(blocks_per_row, [load(path) for path in paths]...)

    save(output_path, new_image)

    return nothing

end

"""
    rangeCut!(
        raw_values::Vector{<:Number},
        range::Tuple{<:Number,<:Number};
        <keyword arguments>
    )::Bool

Delete every element in `raw_values` that is outside the given `range`.

# Arguments

  - `raw_values::Vector{<:Number}`: Dataset that will be pruned.
  - `range::Tuple{<:Number,<:Number}`: The range in question.
  - `keep_edges::Bool=true`: If the edges of the range will be kept.
  - `min_left::Int=1`: Minimum number of values that need to be left after pruning to proceed with the transformation.

# Returns

  - If a transformation was performed.
"""
function rangeCut!(
    raw_values::Vector{<:Number},
    range::Tuple{<:Number,<:Number};
    keep_edges::Bool=true,
    min_left::Int=1,
)::Bool

    # Shortcut computation for special cases
    !(isempty(raw_values) || all(isinf.(range))) || return false

    if keep_edges

        # Check that after the transformation at least `min_left` elements will be left
        count(x -> range[1] <= x <= range[2], raw_values) >= min_left || return false

        # Delete element outside of the provided range
        filter!(x -> range[1] <= x <= range[2], raw_values)

    else

        # Check that after the transformation at least `min_left` elements will be left
        count(x -> range[1] < x < range[2], raw_values) >= min_left || return false

        # Delete element outside of the provided range
        filter!(x -> range[1] < x < range[2], raw_values)

    end

    return true

end

"""
    rangeCut!(
        m_data::Vector{<:Number},
        s_data::Vector,
        range::Tuple{<:Number,<:Number};
        <keyword arguments>
    )::Bool

Delete every element in `m_data` that is outside the given `range`.

Every corresponding element in `s_data` (i.e. with the same index) will be deleted too.

# Arguments

  - `m_data::Vector{<:Number}`: Master dataset that will be pruned.
  - `s_data::Vector`: Slave dataset that will be pruned according to which values of `m_data` are outside `range`.
  - `range::Tuple{<:Number,<:Number}`: The range in question.
  - `keep_edges::Bool=true`: If the edges of the range will be kept.
  - `min_left::Int=1`: Minimum number of values that need to be left in the master dataset after pruning to proceed with the transformation.

# Returns

  - If a transformation was performed.
"""
function rangeCut!(
    m_data::Vector{<:Number},
    s_data::Vector,
    range::Tuple{<:Number,<:Number};
    keep_edges::Bool=true,
    min_left::Int=1,
)::Bool

    # Shortcut computation for special cases
    !(isempty(m_data) || all(isinf.(range))) || return false

    (
        length(s_data) >= length(m_data) ||
        throw(ArgumentError("rangeCut!: `s_data` must have at least as many elements as `m_data`, \
        but I got length(`s_data`) = $(length(s_data)) < length(`m_data`) = $(length(m_data))"))
    )

    if keep_edges

        # Find the elements outside of the provided range
        idxs = map(x -> x < range[1] || x > range[2], m_data)

        # Check that after the transformation at least `min_left` elements will be left
        count(.!idxs) >= min_left || return false

        # Delete element outside of the provided range
        deleteat!(m_data, idxs)
        deleteat!(s_data, idxs)

    else

        # Find the elements outside of the provided range
        idxs = map(x -> x <= range[1] || x >= range[2], m_data)

        # Check that after the transformation at least `min_left` elements will be left
        count(.!idxs) >= min_left || return false

        # Delete element outside of the provided range
        deleteat!(m_data, idxs)
        deleteat!(s_data, idxs)

    end

    return true

end

"""
    sanitizeData!(
        raw_values::Vector{<:Number};
        <keyword arguments>
    )::NTuple{2,Bool}

Do the following transformations over `raw_values`, in order:

  - Trim it to fit within the domain of the function `func_domain`.
  - Trim it to fit within `range`.
  - Scale it down by a factor of 10^`exp_factor`.

By default, no transformation is done.

# Arguments

  - `raw_values::Vector{<:Number}`: Dataset to be sanitized.
  - `func_domain::Function=identity`: `raw_values` will be trimmed to fit within the domain of the function `func_domain`. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `range::Tuple{<:Number,<:Number}=(-Inf, Inf)`: Every element in `raw_values` that falls outside of `range` will be deleted.
  - `keep_edges::Bool=true`: If the edges of `range` will be kept.
  - `min_left::Int=1`: Minimum number of values that need to be left after each transformation to proceed with it.
  - `exp_factor::Int=0`: Every element in `raw_values` will be divided by 10^`exp_factor`.

# Returns

  - A tuple with two flags:

      + If `raw_values` was mutated to fit within the domain of `func_domain`.
      + If `raw_values` was mutated to fit within `range`.
"""
function sanitizeData!(
    raw_values::Vector{<:Number};
    func_domain::Function=identity,
    range::Tuple{<:Number,<:Number}=(-Inf, Inf),
    keep_edges::Bool=true,
    min_left::Int=1,
    exp_factor::Int=0,
)::NTuple{2,Bool}

    !isempty(raw_values) || return false, false

    d_unit = unit(first(raw_values))

    # Trim `raw_values` to fit within the domain of `func_domain`
    if func_domain ∈ [identity, Makie.pseudolog10, Makie.Symlog10]

        domain_flag = false

    elseif func_domain == sqrt

        domain_flag = rangeCut!(raw_values, (0.0, Inf) .* d_unit; keep_edges=true, min_left)

    elseif func_domain == Makie.logit

        domain_flag = rangeCut!(raw_values, (0.0, 1.0) .* d_unit; keep_edges=false, min_left)

    elseif func_domain ∈ [log, log2, log10]

        domain_flag = rangeCut!(raw_values, (0.0, Inf) .* d_unit; keep_edges=false, min_left)

    else

        throw(ArgumentError("sanitizeData!: The function $(func_domain) is not supported. See \
        the list of supported scaling functions in the [Makie](https://docs.makie.org/stable/) \
        documentation"))

    end

    # Trim `raw_values` to fit within `range`
    range_flag = rangeCut!(raw_values, range; keep_edges, min_left)

    (
        !(isa(raw_values, Vector{<:Integer}) && !iszero(exp_factor) && logging[]) ||
        @warn("sanitizeData!: Elements of `raw_values` are of type `Integer`, this may result \
        in errors or unwanted truncation when using `exp_factor` != 0")
    )

    # Scale `raw_values` down by a factor of 10^`exp_factor`
    iszero(exp_factor) || (raw_values ./= exp10(exp_factor))

    return domain_flag, range_flag

end

"""
    sanitizeData!(
        x_data::Vector{<:Number},
        y_data::Vector{<:Number};
        <keyword arguments>
    )::NTuple{4,Bool}

Do the following transformations over `x_data` and `y_data`, in order:

  - Trim them to fit within the domain of the functions `func_domain[1]` and `func_domain[2]`, respectively.
  - Trim them to fit within `range[1]` and `range[2]`, respectively.
  - Scale them down by a factor 10^`exp_factor[1]` and 10^`exp_factor[2]`, respectively.

By default, no transformation is done.

!!! note

    The datasets must have the same length, and any operation that deletes an element, will delete the corresponding element (i.e. with the same index) in the other dataset, so that the dataset will remain of equal length.

# Arguments

  - `x_data::Vector{<:Number}`: First dataset to be sanitized.
  - `y_data::Vector{<:Number}`: Second dataset to be sanitized.
  - `func_domain::NTuple{2,Function}=(identity, identity)`: `x_data` will be trimmed to fit within the domain of the function `func_domain[1]`, and `y_data` will be trimmed to fit within the domain of the function `func_domain[2]`. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `range::Tuple{Tuple{<:Number,<:Number},Tuple{<:Number,<:Number}}=((-Inf, Inf), (-Inf, Inf))`: Every element in `x_data` that falls outside of `range[1]` will be deleted, and every element in `y_data` that falls outside of `range[2]` will be deleted.
  - `keep_edges::NTuple{2,Bool}=(true, true)`: If the edges of each corresponding `range` will be kept.
  - `min_left::Int=1`: Minimum number of values that need to be left in each dataset after any of the transformations to proceed with them.
  - `exp_factor::NTuple{2,Int}=(0, 0)`: Every element in `x_data` will be divided by 10^`exp_factor[1]`, and every element in `y_data` will be divided by 10^`exp_factor[2]`.

# Returns

  - A tuple with four flags:

      + If `x_data` was successfully modified to fit within the domain of `func_domain[1]`.
      + If `y_data` was successfully modified to fit within the domain of `func_domain[2]`.
      + If `x_data` was successfully modified to fit within `range[1]`.
      + If `y_data` was successfully modified to fit within `range[2]`.
"""
function sanitizeData!(
    x_data::Vector{<:Number},
    y_data::Vector{<:Number};
    func_domain::NTuple{2,Function}=(identity, identity),
    range::Tuple{Tuple{<:Number,<:Number},Tuple{<:Number,<:Number}}=((-Inf, Inf), (-Inf, Inf)),
    keep_edges::NTuple{2,Bool}=(true, true),
    min_left::Int=1,
    exp_factor::NTuple{2,Int}=(0, 0),
)::NTuple{4,Bool}

    (
        length(x_data) == length(y_data) ||
        throw(ArgumentError("sanitizeData!: `x_data` and `y_data` must have the same length, \
        but I got length(x_data) = $(length(x_data)) != length(y_data) = $(length(y_data))"))
    )

    x_unit = isempty(x_data) ? Unitful.NoUnits : unit(first(x_data))

    # Trim the data to fit within the domain of `func_domain[1]`
    if func_domain[1] ∈ [identity, Makie.pseudolog10, Makie.Symlog10]

        x_domain_flag = false

    elseif func_domain[1] == sqrt

        x_domain_flag = rangeCut!(x_data, y_data, (0.0, Inf) .* x_unit; keep_edges=true, min_left)

    elseif func_domain[1] == Makie.logit

        x_domain_flag = rangeCut!(x_data, y_data, (0.0, 1.0) .* x_unit; keep_edges=false, min_left)

    elseif func_domain[1] ∈ [log, log2, log10]

        x_domain_flag = rangeCut!(x_data, y_data, (0.0, Inf) .* x_unit; keep_edges=false, min_left)

    else

        throw(ArgumentError("sanitizeData!: The function $(func_domain[1]) is not supported. See \
        the list of supported scaling functions in the [Makie](https://docs.makie.org/stable/) \
        documentation"))

    end

    y_unit = isempty(y_data) ? Unitful.NoUnits : unit(first(y_data))

    # Trim the data to fit within the domain of `func_domain[2]`
    if func_domain[2] ∈ [identity, Makie.pseudolog10, Makie.Symlog10]

        y_domain_flag = false

    elseif func_domain[2] == sqrt

        y_domain_flag = rangeCut!(y_data, x_data, (0.0, Inf) .* y_unit; keep_edges=true, min_left)

    elseif func_domain[2] == Makie.logit

        y_domain_flag = rangeCut!(y_data, x_data, (0.0, 1.0) .* y_unit; keep_edges=false, min_left)

    elseif func_domain[2] ∈ [log, log2, log10]

        y_domain_flag = rangeCut!(y_data, x_data, (0.0, Inf) .* y_unit; keep_edges=false, min_left)

    else

        throw(ArgumentError("sanitizeData!: The function $(func_domain[2]) is not supported. See \
        the list of supported scaling functions in the [Makie](https://docs.makie.org/stable/) \
        documentation"))

    end

    # Trim data to fit within `range[1]`
    x_range_flag = rangeCut!(x_data, y_data, range[1]; keep_edges=keep_edges[1], min_left)

    # Trim data to fit within `range[2]`
    y_range_flag = rangeCut!(y_data, x_data, range[2]; keep_edges=keep_edges[2], min_left)

    (
        !(isa(x_data, Vector{<:Integer}) && !iszero(exp_factor[1]) && logging[]) ||
        @warn("sanitizeData!: Elements of `x_data` are of type Integer, this may result \
        in errors or unwanted truncation when using `exp_factor[1]` != 0")
    )

    (
        !(isa(y_data, Vector{<:Integer}) && !iszero(exp_factor[2]) && logging[]) ||
        @warn("sanitizeData!: Elements of `y_data` are of type Integer, this may result \
        in errors or unwanted truncation when using `exp_factor[2]` != 0")
    )

    # Scale the data down by the factors `exp_factor`
    iszero(exp_factor[1]) || (x_data ./= exp10(exp_factor[1]))
    iszero(exp_factor[2]) || (y_data ./= exp10(exp_factor[2]))

    return x_domain_flag, y_domain_flag, x_range_flag, y_range_flag

end

"""
    scaledBins(
        values::Vector{<:Number},
        n_bins::Int;
        <keyword arguments>
    )::Vector{Float64}

Compute a set of bin edges, to encompass a given list of values.

# Arguments

  - `values::Vector{<:Number}`: Values to be binned.
  - `n_bins::Int`: Number of bins.
  - `scaling::Function=identity`: Scaling function. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity.
  - `limits::Tuple{<:Number,<:Number}=(-Inf, Inf)`: Set it to a value different than `(-Inf, Inf)` if you want to fix the limits of the bins.

# Returns

  - A sorted list of bin edges.
"""
function scaledBins(
    values::Vector{<:Number},
    n_bins::Int;
    scaling::Function=identity,
    limits::Tuple{<:Number,<:Number}=(-Inf, Inf),
)::Vector{<:Number}

    (
        !(limits[1] > limits[2]) ||
        throw(ArgumentError("scaledBins: `limits` must be (min, max), but I got `limits`[1] = \
        $(limits[1]) > `limits`[2] = $(limits[2])"))
    )

    # Compute the limits of the bins
    min = isinf(limits[1]) ? scaling(ustrip(minimum(values))) : scaling(ustrip(limits[1]))
    max = isinf(limits[2]) ? scaling(ustrip(maximum(values))) : scaling(ustrip(limits[2]))

    # For a small range, increase it by 0.2 * abs(max)
    if ((range = max - min) <= 1e-4 * abs(max)) && isinf(limits[1]) && isinf(limits[2])
        range += 0.2 * abs(max)
    end

    # Compute the width of the bins
    width = range / n_bins

    # Get the inverse function of `scaling`
    inverse = Makie.inverse_transform(scaling)

    # Get the unit of the values
    v_unit = unit(first(values))

    return [inverse(min + width * i) * v_unit for i in 0:n_bins]

end

"""
    listHistogram1D(
        positions::Vector{<:Number},
        values::Vector{<:Number},
        grid::Union{LinearGrid,CircularGrid},
    )::Vector{Vector{<:Number}}

Compute a 1D histogram of `values`, returning the full list of `values` within each bin.

# Arguments

  - `positions::Vector{<:Number}`: Positions of the `values` within a 1D axis. This determines to which bin each value will be added.
  - `values::Vector{<:Number}`: The values that will be added up in each bin, according to their `positions`.
  - `grid::Union{LinearGrid,CircularGrid}`: A linear or circular grid.

# Returns

    - A vector with the lists of `values` within each bin.
"""
function listHistogram1D(
    positions::Vector{<:Number},
    values::Vector{<:Number},
    grid::Union{LinearGrid,CircularGrid},
)::Vector{Vector{<:Number}}

    (
        length(values) == length(positions) ||
        throw(ArgumentError("listHistogram1D: `values` must have as many elements as \
        `positions`, but I got length(`values`) = $(length(values)) and \
        length(`positions`) = $(length(positions))"))
    )

    n_bins = length(grid.grid)

    if grid.log
        l_u = unit(first(grid.ticks))
        positions = log10.(ustrip.(l_u, positions))
        p_min = log10(ustrip(grid.ticks[1]))
        p_max = log10(ustrip(l_u, grid.ticks[end]))
    else
        p_min = grid.ticks[1]
        p_max = grid.ticks[end]
    end

    # Compute the bin width
    width = (p_max - p_min) / n_bins

    # Allocate memory
    histogram = [eltype(values)[] for _ in 1:n_bins]

    # Compute the histogram; ignoring NaNs and positions outside the grid range
    for (position, value) in zip(positions, values)

        if isnan(position) || isnan(value)
            continue
        elseif position < p_min || p_max < position
            continue
        elseif position == p_min
            idx = 1
        elseif position == p_max
            idx = n_bins
        else
            idx = ceil(Int, (position - p_min) / width)
        end

        push!(histogram[idx], value)

    end

    return histogram

end

"""
    listHistogram1D(
        positions::Vector{<:Number},
        values::Vector{<:Number},
        edges::Vector{<:Number},
    )::Vector{Vector{<:Number}}

Compute a 1D histogram of `values`, returning the full list of `values` within each bin.

# Arguments

  - `positions::Vector{<:Number}`: Positions of the `values` within a 1D axis. This determines to which bin each value will be added.
  - `values::Vector{<:Number}`: The values that will be added up in each bin, according to their `positions`.
  - `edges::Vector{<:Number}`: A sorted list of bin edges.

# Returns

  - A vector with the lists of `values` within each bin.
"""
function listHistogram1D(
    positions::Vector{<:Number},
    values::Vector{<:Number},
    edges::Vector{<:Number},
)::Vector{Vector{<:Number}}

    (
        length(values) == length(positions) ||
        throw(ArgumentError("listHistogram1D: `values` must have as many elements as \
        `positions`, but I got length(`values`) = $(length(values)) and \
        length(`positions`) = $(length(positions))"))
    )

    issorted(edges) || sort!(edges)

    n_bins = length(edges) - 1

    p_min = first(edges)
    p_max = last(edges)

    # Allocate memory
    histogram = [eltype(values)[] for _ in 1:n_bins]

    # Compute the histogram; ignoring NaNs and positions outside of range
    for (position, value) in zip(positions, values)

        if isnan(position) || isnan(value)
            continue
        elseif position < p_min || p_max < position
            continue
        elseif position == p_min
            idx = 1
        else
            idx = searchsortedfirst(edges, position) - 1
        end

        push!(histogram[idx], value)

    end

    return histogram

end

"""
    listHistogram3D(positions::Matrix{<:Number}, grid::CubicGrid)::Array{Vector{Int},3}

Compute a 3D histogram of `positions`, returning the full list of indices within each bin.

# Arguments

  - `positions::Matrix{<:Number}`: Positions of the points in the grid. Each column correspond to a point and each row is a dimension. This determines to which bin the index of each point will be added.
  - `grid::CubicGrid`: A cubic grid.

# Returns

    - A tensor with the indices of the points within each bin.
"""
function listHistogram3D(positions::Matrix{<:Number}, grid::CubicGrid)::Array{Vector{Int},3}

    # Half bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    x_borders = (grid.x_ticks[1] - h_bin_width, grid.x_ticks[end] + h_bin_width)
    y_borders = (grid.y_ticks[1] - h_bin_width, grid.y_ticks[end] + h_bin_width)
    z_borders = (grid.z_ticks[1] - h_bin_width, grid.z_ticks[end] + h_bin_width)

    # Allocate memory
    histogram = Array{Vector{Int}}(undef, size(grid.grid))
    for i in eachindex(histogram)
        histogram[i] = Int[]
    end

    for (idx, point) in pairs(eachcol(positions))

        !any(isnan, point) || continue

        x = point[1]
        y = point[2]
        z = point[3]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue
        !(z > z_borders[2] || z < z_borders[1]) || continue

        if x == x_borders[1]
            i_x = 1
        elseif x == x_borders[2]
            i_x = grid.n_bins
        else
            i_x = ceil(Int, (x - x_borders[1]) / grid.bin_width)
        end

        if y == y_borders[1]
            i_y = 1
        elseif y == y_borders[2]
            i_y = grid.n_bins
        else
            i_y = ceil(Int, (y - y_borders[1]) / grid.bin_width)
        end

        if z == z_borders[1]
            i_z = 1
        elseif z == z_borders[2]
            i_z = grid.n_bins
        else
            i_z = ceil(Int, (z - z_borders[1]) / grid.bin_width)
        end

        push!(histogram[grid.n_bins - i_y + 1, i_x, i_z], idx)

    end

    return histogram

end

"""
    histogram1D(
        positions::Vector{<:Number},
        values::Vector{<:Number},
        grid::Union{LinearGrid,CircularGrid};
        <keyword arguments>
    )::Vector{<:Number}

Compute a 1D histogram of `values`.

# Arguments

  - `positions::Vector{<:Number}`: Positions of the `values` within a 1D axis. This determines to which bin each value will be added.
  - `values::Vector{<:Number}`: The values that will be added up in each bin, according to their `positions`.
  - `grid::Union{LinearGrid,CircularGrid}`: A linear or circular grid.
  - `total::Bool=true`: If the sum (default) or the mean of `values` will be computed for each bin.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A vector with the histogram values.
"""
function histogram1D(
    positions::Vector{<:Number},
    values::Vector{<:Number},
    grid::Union{LinearGrid,CircularGrid};
    total::Bool=true,
    empty_nan::Bool=true,
)::Vector{<:Number}

    (
        length(values) == length(positions) ||
        throw(ArgumentError("histogram1D: `values` must have as many elements as \
        `positions`, but I got length(`values`) = $(length(values)) and \
        length(`positions`) = $(length(positions))"))
    )

    n_bins = length(grid.grid)

    if grid.log
        l_u = unit(first(grid.ticks))
        positions = log10.(ustrip.(l_u, positions))
        p_min = log10(ustrip(grid.ticks[1]))
        p_max = log10(ustrip(l_u, grid.ticks[end]))
    else
        p_min = grid.ticks[1]
        p_max = grid.ticks[end]
    end

    # Compute the bin width
    width = (p_max - p_min) / n_bins

    # Allocate memory
    histogram = zeros(eltype(values), n_bins)
    counts = zeros(Int, n_bins)

    # Compute the histogram; ignoring NaNs and positions outside the grid range
    for (position, value) in zip(positions, values)

        if isnan(position) || isnan(value)
            continue
        elseif position < p_min || p_max < position
            continue
        elseif position == p_min
            idx = 1
        elseif position == p_max
            idx = n_bins
        else
            idx = ceil(Int, (position - p_min) / width)
        end

        histogram[idx] += value
        counts[idx] += 1

    end

    if empty_nan
        # Set empty bins to NaN
        nan = NaN * unit(first(values))

        for (i, count) in pairs(counts)
            if iszero(count)
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        for (i, count) in pairs(counts)
            if !iszero(count)
                histogram[i] /= count
            end
        end
    end

    return histogram

end

"""
    histogram1D(
        positions::Vector{<:Number},
        values::Vector{<:Number},
        edges::Vector{<:Number};
        <keyword arguments>
    )::Vector{<:Number}

Compute a 1D histogram of `values`.

# Arguments

  - `positions::Vector{<:Number}`: Positions of the `values` within a 1D axis. This determines to which bin each value will be added.
  - `values::Vector{<:Number}`: The values that will be added up in each bin, according to their `positions`.
  - `edges::Vector{<:Number}`: A sorted list of bin edges.
  - `total::Bool=true`: If the sum (default) or the mean of `values` will be computed for each bin.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A vector with the histogram values.
"""
function histogram1D(
    positions::Vector{<:Number},
    values::Vector{<:Number},
    edges::Vector{<:Number};
    total::Bool=true,
    empty_nan::Bool=true,
)::Vector{<:Number}

    (
        length(values) == length(positions) ||
        throw(ArgumentError("histogram1D: `values` must have as many elements as \
        `positions`, but I got length(`values`) = $(length(values)) and \
        length(`positions`) = $(length(positions))"))
    )

    issorted(edges) || sort!(edges)

    n_bins = length(edges) - 1

    p_min = first(edges)
    p_max = last(edges)

    # Allocate memory
    histogram = zeros(eltype(values), n_bins)
    counts = zeros(Int, n_bins)

    # Compute the histogram; ignoring NaNs and positions outside of range
    for (position, value) in zip(positions, values)

        if isnan(position) || isnan(value)
            continue
        elseif position < p_min || p_max < position
            continue
        elseif position == p_min
            idx = 1
        else
            idx = searchsortedfirst(edges, position) - 1
        end

        histogram[idx] += value
        counts[idx] += 1

    end

    if empty_nan
        # Set empty bins to NaN
        nan = NaN * unit(first(values))

        for (i, count) in pairs(counts)
            if iszero(count)
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        for (i, count) in pairs(counts)
            if !iszero(count)
                histogram[i] /= count
            end
        end
    end

    return histogram

end

"""
    histogram1D(
        positions::Vector{<:Number},
        grid::Union{LinearGrid,CircularGrid},
    )::Vector{Int}

Compute a 1D histogram of `positions`.

# Arguments

  - `positions::Vector{<:Number}`: Values for which the histogram will be constructed.
  - `grid::Union{LinearGrid,CircularGrid}`: A linear or circular grid.

# Returns

  - A vector with the counts.
"""
function histogram1D(
    positions::Vector{<:Number},
    grid::Union{LinearGrid,CircularGrid},
)::Vector{Int}

    n_bins = length(grid.grid)

    if grid.log
        l_u = unit(first(grid.ticks))
        positions = log10.(ustrip.(l_u, positions))
        p_min = log10(ustrip(grid.ticks[1]))
        p_max = log10(ustrip(l_u, grid.ticks[end]))
    else
        p_min = grid.ticks[1]
        p_max = grid.ticks[end]
    end

    # Compute the bin width
    width = (p_max - p_min) / n_bins

    # Allocate memory
    histogram = zeros(Int, n_bins)

    # Compute the histogram; ignoring NaNs and positions outside the grid range
    for position in positions

        if isnan(position)
            continue
        elseif position < p_min || p_max < position
            continue
        elseif position == p_min
            idx = 1
        elseif position == p_max
            idx = n_bins
        else
            idx = ceil(Int, (position - p_min) / width)
        end

        histogram[idx] += 1

    end

    return histogram

end

"""
    histogram1D(positions::Vector{<:Number}, edges::Vector{<:Number})::Vector{Int}

Compute a 1D histogram of `positions`.

# Arguments

  - `positions::Vector{<:Number}`: Values for which the histogram will be constructed.
  - `edges::Vector{<:Number}`: A sorted list of bin edges.

# Returns

  - A vector with the counts.
"""
function histogram1D(positions::Vector{<:Number}, edges::Vector{<:Number})::Vector{Int}

    issorted(edges) || sort!(edges)

    n_bins = length(edges) - 1

    p_min = first(edges)
    p_max = last(edges)

    # Allocate memory
    histogram = zeros(Int, n_bins)

    # Compute the histogram; ignoring NaNs and positions outside of range
    for position in positions

        if isnan(position)
            continue
        elseif position < p_min || p_max < position
            continue
        elseif position == p_min
            idx = 1
        else
            idx = searchsortedfirst(edges, position) - 1
        end

        histogram[idx] += 1

    end

    return histogram

end

"""
    histogram2D(
        positions::Matrix{<:Number},
        values::Vector{<:Number},
        grid::SquareGrid;
        <keyword arguments>
    )::Matrix{<:Number}

Compute a 2D histogram of `values`.

# Arguments

  - `positions::Matrix{<:Number}`: Positions of the values in the grid. Each column correspond to a value and each row is a dimension. This determines to which bin each value will be added.
  - `values::Vector{<:Number}`: The values that will be added up in each square bin, according to their `positions`.
  - `grid::SquareGrid`: A square grid.
  - `total::Bool=true`: If the sum (default) or the mean of `values` will be computed in each bin.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A matrix with the histogram values.
"""
function histogram2D(
    positions::Matrix{<:Number},
    values::Vector{<:Number},
    grid::SquareGrid;
    total::Bool=true,
    empty_nan::Bool=true,
)::Matrix{<:Number}

    (
        length(values) == size(positions, 2) ||
        throw(ArgumentError("histogram2D: `values` must have as many elements as `positions` \
        has columns, but I got length(`values`) = $(length(values)) and size(`positions`, 2) = \
        $(size(positions, 2))"))
    )

    # Half bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    x_borders = (grid.x_ticks[1] - h_bin_width, grid.x_ticks[end] + h_bin_width)
    y_borders = (grid.y_ticks[1] - h_bin_width, grid.y_ticks[end] + h_bin_width)

    # Allocate memory
    histogram = zeros(eltype(values), size(grid.grid))
    counts = zeros(Int, size(grid.grid))

    for (i, point) in pairs(eachcol(positions))

        !isnan(values[i])  || continue
        !any(isnan, point) || continue

        x = point[1]
        y = point[2]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue

        if x == x_borders[1]
            i_x = 1
        elseif x == x_borders[2]
            i_x = grid.n_bins
        else
            i_x = ceil(Int, (x - x_borders[1]) / grid.bin_width)
        end

        if y == y_borders[1]
            i_y = 1
        elseif y == y_borders[2]
            i_y = grid.n_bins
        else
            i_y = ceil(Int, (y - y_borders[1]) / grid.bin_width)
        end

        histogram[grid.n_bins - i_y + 1, i_x] += values[i]
        counts[grid.n_bins - i_y + 1, i_x] += 1

    end

    if empty_nan
        # Set empty bins to NaN
        nan = NaN * unit(first(values))

        for i in eachindex(histogram)
            if iszero(counts[i])
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        for i in eachindex(histogram)
            if !iszero(counts[i])
                histogram[i] /= counts[i]
            end
        end
    end

    return histogram

end

"""
    histogram2D(
        positions::Matrix{<:Number},
        values::Vector{<:Number},
        x_edges::Vector{<:Number},
        y_edges::Vector{<:Number};
        <keyword arguments>
    )::Matrix{<:Number}

Compute a 2D histogram of `values`.

# Arguments

  - `positions::Matrix{<:Number}`: Positions of the values in the grid. Each column correspond to a value and each row is a dimension. This determines to which bin each value will be added.
  - `values::Vector{<:Number}`: The values that will be added up in each square bin, according to their `positions`.
  - `x_edges::Vector{<:Number}`: A sorted list of bin edges for the x axis.
  - `y_edges::Vector{<:Number}`: A sorted list of bin edges for the y axis.
  - `total::Bool=true`: If the sum (default) or the mean of `values` will be computed in each bin.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A matrix with the histogram values.
"""
function histogram2D(
    positions::Matrix{<:Number},
    values::Vector{<:Number},
    x_edges::Vector{<:Number},
    y_edges::Vector{<:Number};
    total::Bool=true,
    empty_nan::Bool=true,
)::Matrix{<:Number}

    (
        length(values) == size(positions, 2) ||
        throw(ArgumentError("histogram2D: `values` must have as many elements as `positions` \
        has columns, but I got length(`values`) = $(length(values)) and size(`positions`, 2) = \
        $(size(positions, 2))"))
    )

    issorted(x_edges) || sort!(x_edges)
    issorted(y_edges) || sort!(y_edges)

    n_x_bins = length(x_edges) - 1
    n_y_bins = length(y_edges) - 1

    x_borders = (first(x_edges), last(x_edges))
    y_borders = (first(y_edges), last(y_edges))

    # Allocate memory
    histogram = zeros(eltype(values), (n_x_bins, n_y_bins))
    counts = zeros(Int, (n_x_bins, n_y_bins))

    for (i, point) in pairs(eachcol(positions))

        !isnan(values[i])  || continue
        !any(isnan, point) || continue

        x = point[1]
        y = point[2]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue

        if x == x_borders[1]
            i_x = 1
        else
            i_x = searchsortedfirst(x_edges, x) - 1
        end

        if y == y_borders[1]
            i_y = 1
        else
            i_y = searchsortedfirst(y_edges, y) - 1
        end

        histogram[n_y_bins - i_y + 1, i_x] += values[i]
        counts[n_y_bins - i_y + 1, i_x] += 1

    end

    if empty_nan
        # Set empty bins to NaN
        nan = NaN * unit(first(values))

        for i in eachindex(histogram)
            if iszero(counts[i])
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        for i in eachindex(histogram)
            if !iszero(counts[i])
                histogram[i] /= counts[i]
            end
        end
    end

    return histogram

end

"""
    histogram2D(positions::Matrix{<:Number}, grid::SquareGrid)::Matrix{Int}

Compute a 2D histogram of `positions`.

# Arguments

  - `positions::Matrix{<:Number}`: Values for which the histogram will be constructed.
  - `grid::SquareGrid`: A square grid.

# Returns

  - A matrix with the counts.
"""
function histogram2D(positions::Matrix{<:Number}, grid::SquareGrid)::Matrix{Int}

    # Half bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    x_borders = (grid.x_ticks[1] - h_bin_width, grid.x_ticks[end] + h_bin_width)
    y_borders = (grid.y_ticks[1] - h_bin_width, grid.y_ticks[end] + h_bin_width)

    # Allocate memory
    histogram = zeros(Int, size(grid.grid))

    for point in eachcol(positions)

        !any(isnan, point) || continue

        x = point[1]
        y = point[2]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue

        if x == x_borders[1]
            i_x = 1
        elseif x == x_borders[2]
            i_x = grid.n_bins
        else
            i_x = ceil(Int, (x - x_borders[1]) / grid.bin_width)
        end

        if y == y_borders[1]
            i_y = 1
        elseif y == y_borders[2]
            i_y = grid.n_bins
        else
            i_y = ceil(Int, (y - y_borders[1]) / grid.bin_width)
        end

        histogram[grid.n_bins - i_y + 1, i_x] += 1

    end

    return histogram

end

"""
    histogram2D(
        positions::Matrix{<:Number},
        x_edges::Vector{<:Number},
        y_edges::Vector{<:Number},
    )::Matrix{Int}

Compute a 2D histogram of `positions`.

# Arguments

  - `positions::Matrix{<:Number}`: Values for which the histogram will be constructed.
  - `x_edges::Vector{<:Number}`: A sorted list of bin edges for the x axis.
  - `y_edges::Vector{<:Number}`: A sorted list of bin edges for the y axis.

# Returns

  - A matrix with the counts.
"""
function histogram2D(
    positions::Matrix{<:Number},
    x_edges::Vector{<:Number},
    y_edges::Vector{<:Number},
)::Matrix{Int}

    issorted(x_edges) || sort!(x_edges)
    issorted(y_edges) || sort!(y_edges)

    n_x_bins = length(x_edges) - 1
    n_y_bins = length(y_edges) - 1

    x_borders = (first(x_edges), last(x_edges))
    y_borders = (first(y_edges), last(y_edges))

    # Allocate memory
    histogram = zeros(Int, (n_x_bins, n_y_bins))

    for point in eachcol(positions)

        !any(isnan, point) || continue

        x = point[1]
        y = point[2]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue

        if x == x_borders[1]
            i_x = 1
        else
            i_x = searchsortedfirst(x_edges, x) - 1
        end

        if y == y_borders[1]
            i_y = 1
        else
            i_y = searchsortedfirst(y_edges, y) - 1
        end

        histogram[n_y_bins - i_y + 1, i_x] += 1

    end

    return histogram

end

"""
    histogram3D(
        positions::Matrix{<:Number},
        values::Vector{<:Number},
        grid::CubicGrid;
        <keyword arguments>
    )::Array{<:Number,3}

Compute a 3D histogram of `values`.

# Arguments

  - `positions::Matrix{<:Number}`: Positions of the values in the grid. Each column correspond to a value and each row is a dimension. This determines to which bin each value will be added.
  - `values::Vector{<:Number}`: The values that will be added up in each square bin, according to their `positions`.
  - `grid::CubicGrid`: A cubic grid.
  - `total::Bool=true`: If the sum (default) or the mean of `values` will be computed in each bin.
  - `empty_nan::Bool=true`: If NaN will be put into empty bins, 0 is used otherwise.

# Returns

  - A 3D tensor with the histogram values.
"""
function histogram3D(
    positions::Matrix{<:Number},
    values::Vector{<:Number},
    grid::CubicGrid;
    total::Bool=true,
    empty_nan::Bool=true,
)::Array{<:Number,3}

    (
        length(values) == size(positions, 2) ||
        throw(ArgumentError("histogram3D: `values` must have as many elements as `positions` \
        has columns, but I got length(`values`) = $(length(values)) and size(`positions`, 2) = \
        $(size(positions, 2))"))
    )

    # Half bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    x_borders = (grid.x_ticks[1] - h_bin_width, grid.x_ticks[end] + h_bin_width)
    y_borders = (grid.y_ticks[1] - h_bin_width, grid.y_ticks[end] + h_bin_width)
    z_borders = (grid.z_ticks[1] - h_bin_width, grid.z_ticks[end] + h_bin_width)

    # Allocate memory
    histogram = zeros(eltype(values), size(grid.grid))
    counts = zeros(Int, size(grid.grid))

    for (i, point) in pairs(eachcol(positions))

        !isnan(values[i])  || continue
        !any(isnan, point) || continue

        x = point[1]
        y = point[2]
        z = point[3]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue
        !(z > z_borders[2] || z < z_borders[1]) || continue

        if x == x_borders[1]
            i_x = 1
        elseif x == x_borders[2]
            i_x = grid.n_bins
        else
            i_x = ceil(Int, (x - x_borders[1]) / grid.bin_width)
        end

        if y == y_borders[1]
            i_y = 1
        elseif y == y_borders[2]
            i_y = grid.n_bins
        else
            i_y = ceil(Int, (y - y_borders[1]) / grid.bin_width)
        end

        if z == z_borders[1]
            i_z = 1
        elseif z == z_borders[2]
            i_z = grid.n_bins
        else
            i_z = ceil(Int, (z - z_borders[1]) / grid.bin_width)
        end

        histogram[grid.n_bins - i_y + 1, i_x, i_z] += values[i]
        counts[grid.n_bins - i_y + 1, i_x, i_z] += 1

    end

    if empty_nan
        # Set empty bins to NaN
        nan = NaN * unit(first(values))

        for i in eachindex(histogram)
            if iszero(counts[i])
                histogram[i] = nan
            end
        end
    end

    if !total
        # Compute the mean value instead of just the sum for each bin
        for i in eachindex(histogram)
            if !iszero(counts[i])
                histogram[i] /= counts[i]
            end
        end
    end

    return histogram

end

"""
    histogram3D(positions::Matrix{<:Number}, grid::CubicGrid)::Array{Int,3}

Compute a 3D histogram of `positions`.

# Arguments

  - `positions::Matrix{<:Number}`: Values for which the histogram will be constructed.
  - `grid::CubicGrid`: A cubic grid.

# Returns

  - A 3D tensor with the counts.
"""
function histogram3D(positions::Matrix{<:Number}, grid::CubicGrid)::Array{Int,3}

    # Half bin size
    h_bin_width = grid.bin_width * 0.5

    # Compute the physical position of the grid borders
    x_borders = (grid.x_ticks[1] - h_bin_width, grid.x_ticks[end] + h_bin_width)
    y_borders = (grid.y_ticks[1] - h_bin_width, grid.y_ticks[end] + h_bin_width)
    z_borders = (grid.z_ticks[1] - h_bin_width, grid.z_ticks[end] + h_bin_width)

    # Allocate memory
    histogram = zeros(Int, size(grid.grid))

    for point in eachcol(positions)

        !any(isnan, point) || continue

        x = point[1]
        y = point[2]
        z = point[3]

        !(x > x_borders[2] || x < x_borders[1]) || continue
        !(y > y_borders[2] || y < y_borders[1]) || continue
        !(z > z_borders[2] || z < z_borders[1]) || continue

        if x == x_borders[1]
            i_x = 1
        elseif x == x_borders[2]
            i_x = grid.n_bins
        else
            i_x = ceil(Int, (x - x_borders[1]) / grid.bin_width)
        end

        if y == y_borders[1]
            i_y = 1
        elseif y == y_borders[2]
            i_y = grid.n_bins
        else
            i_y = ceil(Int, (y - y_borders[1]) / grid.bin_width)
        end

        if z == z_borders[1]
            i_z = 1
        elseif z == z_borders[2]
            i_z = grid.n_bins
        else
            i_z = ceil(Int, (z - z_borders[1]) / grid.bin_width)
        end

        histogram[grid.n_bins - i_y + 1, i_x, i_z] += 1

    end

    return histogram

end

"""
    smoothWindow(
        x_data::Vector{<:Number},
        y_data::Vector{<:Number},
        n_bins::Int;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Separate the values of `x_data` in `n_bins` bins and compute the mean value of `x_data` and `y_data` within each one.

# Arguments

  - `x_data::Vector{<:Number}`: x-axis data.
  - `y_data::Vector{<:Number}`: y-axis data.
  - `n_bins::Int`: Number of bins.
  - `scaling::Function=identity`: Scaling function for the x axis. The options are the scaling functions accepted by [Makie](https://docs.makie.org/stable/): log10, log2, log, sqrt, Makie.logit, Makie.Symlog10, Makie.pseudolog10, and identity. All the values of `x_data` must be in the domain of `scaling`.

# Returns

  - A tuple with two vectors, containing the smoothed-out x and y values.
"""
function smoothWindow(
    x_data::Vector{<:Number},
    y_data::Vector{<:Number},
    n_bins::Int;
    scaling::Function=identity,
)::NTuple{2,Vector{<:Number}}

    # Check that the input vectors have the same length
    (
        length(x_data) == length(y_data) ||
        throw(DimensionMismatch("smoothWindow: `x_data` and `y_data` must have the same length, \
        but I got length(`x_data`) = $(length(x_data)) != length(`y_data`) = $(length(y_data))"))
    )

    positions = scaling.(ustrip(x_data))
    grid = CircularGrid(maximum(positions), n_bins; shift=minimum(positions))

    smooth_x_data = histogram1D(positions, x_data, grid; total=false)
    smooth_y_data = histogram1D(positions, y_data, grid; total=false)

    # Remove empty bins
    return filter!(!isnan, smooth_x_data), filter!(!isnan, smooth_y_data)

end

"""
    cubicSplineKernel(q::Real, h::Number)::Number

2D cubic spline kernel.

# Arguments

  - `q::Real`: Relative distance to the neighbor, ``|r - r'| / h``.
  - `h::Number`: Smoothing length.

# Returns

  - The kernel function evaluated at a separation `q` * `h`, and with a smoothing length `h`.

# References

[PySPH documentation](https://pysph.readthedocs.io/en/latest/reference/kernels.html)

J. J. Monaghan (1992). *Smoothed Particle Hydrodynamics*. Annual Review of Astronomy and Astrophysics, **30**, 543-574. [doi:10.1146/annurev.aa.30.090192.002551](https://doi.org/10.1146/annurev.aa.30.090192.002551)

M.B. Liu et al. (2010). *Smoothed Particle Hydrodynamics (SPH): an Overview and Recent Developments*. Archives of Computational Methods in Engineering, **17**, 25–76. [doi:10.1007/s11831-010-9040-7](https://doi.org/10.1007/s11831-010-9040-7)
"""
function cubicSplineKernel(q::Real, h::Number)::Number

    (
        isPositive(q, h) ||
        throw(DomainError("cubicSplineKernel: `q` and `h` must be positive, \
        but I got `q` = $q and `h` = $h"))
    )

    σ3 = 10.0 / (7.0 * area(h))

    if 0 <= q <= 1
        return σ3 * (1.0 - 1.5 * q * q * (1.0 - q * 0.5))
    elseif 1 < q <= 2
        return (σ3 * 0.25) * (2.0 - q)^3.0
    else
        return zero(σ3)
    end

end

"""
    deltas(data::Vector{<:Number})::Vector{<:Number}

Compute the difference between each consecutive pair of elements in `data`.

# Arguments

  - `data::Vector{<:Number}`: Data vector. It has to have at least 2 elements.

# Returns

  - A vector with the difference between each consecutive pair of elements in `data`, the first element is 0 by convention.
"""
function deltas(data::Vector{<:Number})::Vector{<:Number}

    # Allocate memory
    Δd = similar(data)

    # Get the number of elements in `data`
    nd = length(data)

    # Check that `data` has a valid length
    (
        nd >= 2 ||
        throw(ArgumentError("deltas: `data` must have at least 2 elements, but it has $(nd)"))
    )

    # Set the first value to 0, by convention
    Δd[1] = zero(data[1])

    for i in 2:nd
        Δd[i] = data[i] - data[i - 1]
    end

    return Δd

end

"""
    reduceResolution(
        hr_matrix::Matrix{<:Number},
        factor::Int;
        <keyword arguments>
    )::Matrix{<:Number}

Reduce the number of rows and columns of `hr_matrix` by `factor`, averaging or adding up its values.

# Arguments

  - `hr_matrix::Matrix{<:Number}`: Original "high resolution" matrix. It has to be a square matrix.
  - `factor::Int`: Factor by which the number of rows and columns will be reduced. It has to divide the size of `hr_matrix` exactly.
  - `total::Bool=false`: If the sum (`total` = true) or the mean (`total` = false) of the values in each of the old pixels will be used for the new pixels.

# Returns

  - The new smaller matrix, with the average values.
"""
function reduceResolution(
    hr_matrix::Matrix{<:Number},
    factor::Int;
    total::Bool=false,
)::Matrix{<:Number}

    !isone(factor) || return hr_matrix

    r, c = size(hr_matrix)
    (
        r == c ||
        throw(ArgumentError("reduceResolution: `hr_matrix` has to be a square matrix, but it has \
        $(c) columns and $(r) rows"))
    )

    (
        factor >= 1 || throw(ArgumentError("reduceResolution: `factor` must be >= 1, \
        but I got `factor` = $(factor)"))
    )

    (
        r % factor == 0 ||
        throw(ArgumentError("reduceResolution: `factor` must divide the size of `hr_matrix` \
        exactly, but I got number of rows / `factor` = $(r / factor)"))
    )

    # Compute the size of the new matrix
    new_size = r ÷ factor

    # Allocate memory
    lr_matrix = zeros(eltype(hr_matrix), new_size, new_size)

    # Compute the number of old pixels per new pixel
    old_n_pixels = factor * factor

    for i in eachindex(lr_matrix)

        # Compute the row and column of the new matrix corresponding to index i
        row = mod1(i, new_size)
        col = ceil(Int, i / new_size)

        for j in (factor * (row - 1) + 1):(factor * row)
            for k in (factor * (col - 1) + 1):(factor * col)

                if !isnan(hr_matrix[j, k])
                    lr_matrix[i] += hr_matrix[j, k]
                end

            end
        end

    end

    return total ? lr_matrix : lr_matrix ./ old_n_pixels

end

"""
    reduceTicks(hr_ticks::Vector{<:Number}, factor::Int)::Vector{<:Number}

Reduce the number of ticks in `hr_ticks` by `factor` keeping the total length of the axis the same and assuming `hr_ticks` are regularly spaced.

# Arguments

  - `hr_ticks::Vector{<:Number}`: Original "high resolution" list of ticks.
  - `factor::Int`: Factor by which the number of ticks will be reduced. It has to divide the size of `hr_ticks` exactly.

# Returns

  - The new shorter tick list.
"""
function reduceTicks(hr_ticks::Vector{<:Number}, factor::Int)::Vector{<:Number}

    !isone(factor) || hr_ticks

    l = length(hr_ticks)
    (
        l % factor == 0 ||
        throw(ArgumentError("reduceTicks: `factor` must divide the size of `hr_ticks` \
        exactly, but I got length(`hr_ticks`) / `factor` = $(l / factor)"))
    )

    (
        factor >= 1 ||
        throw(ArgumentError("reduceTicks: `factor` must be >= 1, but I got `factor` = $(factor)"))
    )

    # Compute the size of the new vector
    new_size = l ÷ factor

    # Allocate memory
    lr_ticks = similar(hr_ticks, new_size)

    if iseven(factor)

        shift = factor ÷ 2

        for i in eachindex(lr_ticks)

            idx = (i - 1) * factor + shift

            lr_ticks[i] = (hr_ticks[idx] + hr_ticks[idx + 1]) / 2.0

        end

    else

        shift = ceil(Int, factor / 2)

        for i in eachindex(lr_ticks)

            idx = (i - 1) * factor + shift

            lr_ticks[i] = hr_ticks[idx]

        end

    end

    return lr_ticks

end

"""
    projectIntoCircularGrid(
        image::Matrix{<:Number},
        n_bins::Int;
        <keyword arguments>
    )::Vector{<:Number}

Project `image` into a circular grid, averaging the values in each concentric ring.

# Arguments

  - `image::Matrix{<:Number}`: Original matrix. It has to be a square matrix.
  - `n_bins::Int`: Number of bins for the circular grid.
  - `inscribed::Bool=true`: If the circular grid will be inscribed in `image` when doing the projection. If set to false, the matrix will be inscribed into the circular grid instead.
  - `total::Bool=false`: If the sum (`total` = true) or the mean (`total` = false) of the values in each of the old pixels that fall within each ring will be used.

# Returns

  - A vector with the averages of the values in each concentric ring.
"""
function projectIntoCircularGrid(
    image::Matrix{<:Number},
    n_bins::Int;
    inscribed::Bool=true,
    total::Bool=false,
)::Vector{<:Number}

    r, c = size(image)
    (
        r == c ||
        throw(ArgumentError("projectIntoCircularGrid: `image` has to be a square matrix, \
        but it has $(c) columns and $(r) rows"))
    )

    # Construct a square grid center in (0, 0)
    square_grid = SquareGrid(1.0, r)

    # Construct a circular grid center in (0, 0)
    circular_grid = CircularGrid(inscribed ? 0.5 : sqrt(0.5), n_bins)

    # Compute the radial distance to the origin of each pixel in the square grid
    positions = norm.(vec(square_grid.grid))

    non_nan_idxs = findall(!isnan, vec(image))

    profile = histogram1D(
        positions,
        vec(image),
        circular_grid;
        total,
        empty_nan=false,
    )

    return profile

end

"""
    density3DProjection(
        data_dict::Dict,
        grid::CubicGrid,
        quantity::Symbol,
        type::Symbol;
        <keyword arguments>
    )::Array{Float64,3}

Sample the 3D density field of a given quantity using a cubic grid

!!! note

    If the source of the field are particles, a simple 3D histogram is used. If they are Voronoi cells instead, the density of the cell that intersect each voxel is used.

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
  - `type::Symbol`: If the source of the field are `:particles` or Voronoi `:cells`.
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

  - A 3D array with the density at each point of the 3D grid.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function density3DProjection(
    data_dict::Dict,
    grid::CubicGrid,
    quantity::Symbol,
    type::Symbol;
    m_unit::Unitful.Units=u"Msun",
    l_unit::Unitful.Units=u"kpc",
    filter_function::Function=filterNothing,
)::Array{Float64,3}

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
        throw(ArgumentError("density3DProjection: I don't recognize the quantity :$(quantity)"))
    end

    # For comological simulations with comoving units, correct
    # the density so it is always in physical units
    if !PHYSICAL_UNITS && data_dict[:sim_data].cosmological
        # Correction factor for the volume
        # V [physical units] = V [comoving units] * a0^3
        physical_factor = data_dict[:snap_data].scale_factor^3
    else
        physical_factor = 1.0
    end

    # Load the cell/particle positions
    positions = filtered_dd[component]["POS "]

    # Compute the masses of the target quantity
    masses = scatterQty(filtered_dd, quantity)

    # If any of the necessary quantities are missing return an empty density field
    if any(isempty, [masses, positions])
        return fill(NaN, (grid.n_bins, grid.n_bins, grid.n_bins))
    end

    if type == :cells

        # Compute the volume of each cell
        cell_volumes = filtered_dd[component]["MASS"] ./ filtered_dd[component]["RHO "]

        # Compute the densities of the target quantity
        densities = ustrip.(m_unit * l_unit^-3, masses ./ cell_volumes)

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
        density = similar(grid.grid, Float64)

        # Compute the density in each voxel
        for i in eachindex(grid.grid)
            density[i] = densities[idxs[i]]
        end

        # Set bins with a value of 0 to NaN
        replace!(x -> iszero(x) ? NaN : x, density)

    elseif type == :particles

        # Compute the 3D histogram
        density = ustrip.(
            m_unit * l_unit^-3,
            histogram3D(positions, masses, grid; empty_nan=true) ./ grid.bin_volume,
        )

    else

        throw(ArgumentError("density3DProjection: The argument `type` must be :cells or \
        :particles, but I got :$(type)"))

    end

    if logging[]

        log_density = filter(!isnan, log10.(density))

        if isempty(log_density)

            min_max_ρ = (NaN, NaN)
            mean_ρ    = NaN
            meadian_ρ = NaN
            mode_ρ    = NaN

        else

            min_max_ρ = extrema(log_density)
            mean_ρ    = mean(log_density)
            meadian_ρ = median(log_density)
            mode_ρ    = mode(log_density)

        end

        # Print the density range
        @info(
            "\nDensity range - log₁₀(ρ [$(m_unit * l_unit^-3)]) \
            \n  Simulation: $(basename(filtered_dd[:sim_data].path)) \
            \n  Snapshot:   $(filtered_dd[:snap_data].global_index) \
            \n  Quantity:   $(quantity) \
            \n  Type:       $(type) \
            \n  Min - Max:  $(min_max_ρ) \
            \n  Mean:       $(mean_ρ) \
            \n  Median:     $(meadian_ρ) \
            \n  Mode:       $(mode_ρ)"
        )

    end

    return density

end

"""
    computeParticleProfile(
        positions::Matrix{<:Unitful.Length},
        quantity::Vector{<:Number},
        grid::CircularGrid;
        <keyword arguments>
    )::Vector{<:Number}

Compute a profile, using an 1D histogram.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `quantity::Vector{<:Number}`: The profile will be of this quantity.
  - `grid::CircularGrid`: Circular grid.
  - `norm_values::Vector{<:Number}=Number[]`: Values to normalize `quantity`.
  - `flat::Bool=true`: If the profile will be 2D, using rings, or 3D, using spherical shells.
  - `total::Bool=true`: If the sum (default) or the mean of `quantity` will be computed for each bin.
  - `cumulative::Bool=false`: If the profile will be accumulated or not.
  - `density::Bool=false`: If the profile will be of the density of `quantity`.
  - `empty_nan::Bool=true`: If empty bins will be set to NaN, 0 is used otherwise. Be carefull if `empty_nan` = true and `cumulative` = true, because every bin after the first NaN will be set to NaN.

# Returns

  - Vector with the values of the profile.
"""
function computeParticleProfile(
    positions::Matrix{<:Unitful.Length},
    quantity::Vector{<:Number},
    grid::CircularGrid;
    norm_values::Vector{<:Number}=Number[],
    flat::Bool=true,
    total::Bool=true,
    cumulative::Bool=false,
    density::Bool=false,
    empty_nan::Bool=true,
)::Vector{<:Number}

    # Return a null profile if `quantity` is empty
    if isempty(quantity)

        (
            !logging[] ||
            @warn("computeParticleProfile: The vector `quantity` is empty. The profile will be \
            filled with NaNs")
        )

        return fill(NaN, length(grid.grid))

    end

    # Compute the distances of the cells/particles to the center of the grid
    if flat
        distances = computeDistance(positions[1:2, :]; center=grid.center[1:2])
    else
        distances = computeDistance(positions; center=grid.center)
    end

    # Compute the histogram of `quantity`
    if isempty(norm_values)

        profile = histogram1D(distances, quantity, grid; total, empty_nan)

    else

        quantity_histogram = histogram1D(distances, quantity, grid; total, empty_nan)
        norm_values_histogram = histogram1D(distances, norm_values, grid; total, empty_nan=false)

        replace!(x -> iszero(x) ? oneunit(x) : x, norm_values_histogram)

        profile = quantity_histogram ./ norm_values_histogram

    end

    region = flat ? grid.bin_areas : grid.bin_volumes

    if cumulative
        return density ? cumsum(profile) ./ cumsum(region) : cumsum(profile)
    end

    return density ? profile ./ region : profile

end

"""
    computeParticleBandProfile(
        positions::Matrix{<:Unitful.Length},
        quantity::Vector{<:Number},
        grid::CircularGrid;
        <keyword arguments>
    )::NTuple{2,Vector{<:Number}}

Compute a profile of the mean and standard deviation of `quantity`, using an 1D histogram

Empty bins have NaN as mean and standard deviation.

# Arguments

  - `positions::Matrix{<:Unitful.Length}`: Positions of the cells/particles. Each column is a cell/particle and each row a dimension.
  - `quantity::Vector{<:Number}`: The profile will be of this quantity.
  - `grid::CircularGrid`: Circular grid.
  - `flat::Bool=true`: If the profile will be 2D, using rings, or 3D, using spherical shells.

# Returns

  - A tuple with two elements:

      + A vector with the mean value for each bin.
      + A vector with the standard deviation for each bin.
"""
function computeParticleBandProfile(
    positions::Matrix{<:Unitful.Length},
    quantity::Vector{<:Number},
    grid::CircularGrid;
    flat::Bool=true,
)::NTuple{3,Vector{<:Number}}

    if isempty(quantity)

        (
            !logging[] ||
            @warn("computeParticleBandProfile: The vector `quantity` is empty. The profile will be \
            filled with NaNs")
        )

        return fill(NaN, length(grid.grid))

    end

    # Compute the distances of the cells/particles to the center of the grid
    if flat
        distances = computeDistance(positions[1:2, :]; center=grid.center[1:2])
    else
        distances = computeDistance(positions; center=grid.center)
    end

    # Compute the histogram of `quantity`
    histogram = listHistogram1D(distances, quantity, grid)

    return quantile.(histogram, 0.5), quantile.(histogram, 0.25), quantile.(histogram, 0.75)

end

@doc raw"""
    energyIntegrand(a::Real, header::SnapshotHeader)::Float64

The integrand of the integral that converts the scale factor into physical time:

```math
\frac{1}{H\,\sqrt{\mathcal{E}}} \, ,
```

where

```math
\mathcal{E} = \Omega_\Lambda + (1 - \Omega_\Lambda - \Omega_m) \, a^{-2} + \Omega_m \, a^{-3} \, ,
```
```math
H = H_0 \, a \, .
```

# Arguments

  - `a::Real`: Scale factor.
  - `header::SnapshotHeader`: Header of the relevant snapshot file.

# Returns

  - The integrand evaluated at `a`, in $\mathrm{Gyr}$.
"""
function energyIntegrand(a::Real, header::SnapshotHeader)::Float64

    # Return 0 if `a` = 0, as the integrand goes to 0 in the limit a -> 0.
    !iszero(a) || return 0.0

    # Compute Ω_K (curvature)
    omega_K = 1.0 - header.omega_0 - header.omega_l

    # Compute the energy function
    E = header.omega_0 / (a * a * a) + omega_K / (a * a) + header.omega_l

    # Compute the hubble constant in Gyr^-1
    H = header.h0 * HUBBLE_CONSTANT * a

    # Return the integrand, in Gyr
    return 1.0 / (H * sqrt(E))

end

@doc raw"""
Time factor for the SF model, without the fraction factors.

τ_star(ρ_cell)    $\equiv \tau_\mathrm{star}$
τ_rec(ρ_cell)     $\equiv \tau_\mathrm{rec} \, f_i$
τ_cond(ρ_cell, Z) $\equiv \tau_\mathrm{cond} \, (1 - f_s)$

"""
τ_star(ρ_cell) = C_star / sqrt(ρ_cell)
τ_rec(ρ_cell) = C_rec / ρ_cell
τ_cond(ρ_cell, Z) = C_cond / (ρ_cell * (Z + Zeff))

"""
    flattenGrid(cubic_grid::CubicGrid)::SquareGrid

Using a `CubicGrid` construct a `SquareGrid` with the same center, number of bins, and physical side length.

# Arguments

  - `cubic_grid::CubicGrid`: Cubic grid.

# Returns

  - A square grid.
"""
function flattenGrid(cubic_grid::CubicGrid)::SquareGrid

    physical_size = cubic_grid.physical_size
    n_bins = cubic_grid.n_bins

    bin_width  = physical_size / n_bins
    shift = 0.5 * (physical_size - bin_width)

    center = [cubic_grid.x_ticks[1], cubic_grid.y_ticks[1], cubic_grid.z_ticks[1]] .+ shift

    return SquareGrid(physical_size, n_bins; center)

end

@doc raw"""
    bigiel2008(
        ΣH::Vector{<:SurfaceDensity};
        <keyword arguments>
    )::Vector{<:Number}

Kennicutt-Schmidt law for the molecular or neutral gas, taken from a set of observations of nearby galaxies.

From Bigiel et al. (2008) (Section 3.1), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{HI, H_2, gas}}{10 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```
where N is the power-law index, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $10 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `ΣH::Vector{<:SurfaceDensity}`: Values of the molecular or neutral gas surface density, with units.
  - `molecular::Bool=true`: If the x axis will be the area mass density of molecular hydrogen, or, if set to false, the area mass density of neutral hydrogen.
  - `log_output::Bool=true`: If the output will the $\log_{10}$ of the star formation area density, or the star formation area density itself (with units). If `log_output` = true, the implied unit is $\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}$

# Returns

  - The star formation area density.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function bigiel2008(
    ΣH::Vector{<:SurfaceDensity};
    molecular::Bool=true,
    log_output::Bool=true,
)::Vector{<:Number}

    log10ΣH = @. log10(uconvert(Unitful.NoUnits, ΣH / 10.0u"Msun * pc^-2"))

    if molecular
        log10Σsfr = @. A_BIGIEL2008_BF_MOLECULAR + log10ΣH * N_BIGIEL2008_BF_MOLECULAR
    else
        log10Σsfr = @. A_BIGIEL2008_NEUTRAL + log10ΣH * N_BIGIEL2008_NEUTRAL
    end

    if log_output
        return log10Σsfr
    else
        return @. exp10(log10Σsfr ) * u"Msun * yr^-1 * kpc^-2"
    end

end

@doc raw"""
    invBigiel2008(
        Σsfr ::Vector{<:MassFlowDensity};
        <keyword arguments>
    )::Vector{<:Number}

Inverse Kennicutt-Schmidt law for the molecular or neutral gas, taken from a set of observations of nearby galaxies.

From Bigiel et al. (2008) (Section 3.1, Eq. 2), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{HI, H_2, gas}}{10 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \, ,
```
where N is the power-law index, and $A = \log_{10}(a)$, where $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $10 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `Σsfr ::Vector{<:MassFlowDensity}`: Values of the star formation area density, with units.
  - `molecular::Bool=true`: If the output will be the area mass density of molecular hydrogen, or, if set to false, the area mass density of neutral hydrogen.
  - `log_output::Bool=true`: If the output will the $\log_{10}$ of the molecular or neutral gas surface density, or the molecular or neutral gas surface density itself (with units). If `log_output` = true, the implied unit is $10 \, \mathrm{M_\odot \, pc^{-2}}$

# Returns

  - The molecular or neutral gas surface density.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function invBigiel2008(
    Σsfr ::Vector{<:MassFlowDensity};
    molecular::Bool=true,
    log_output::Bool=true,
)::Vector{<:Number}

    log10Σsfr = @. log10(ustrip(u"Msun * yr^-1 * kpc^-2", Σsfr ))

    if molecular
        log10ΣH = @. (log10Σsfr - A_BIGIEL2008_BF_MOLECULAR) / N_BIGIEL2008_BF_MOLECULAR
    else
        log10ΣH = @. (log10Σsfr - A_BIGIEL2008_NEUTRAL) / N_BIGIEL2008_NEUTRAL
    end

    if log_output
        return log10ΣH
    else
        return @. exp10(log10ΣH) * 10.0u"Msun * pc^-2"
    end

end

@doc raw"""
    kennicutt1998(Σgas::Vector{<:SurfaceDensity}; <keyword arguments>)::Vector{<:Number}

Kennicutt-Schmidt law, taken from a set of observations of nearby galaxies.

From Kennicutt (1998) (Section 4, Eq. 4), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{gas}}{1 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \mathrm{M_\odot \, yr^{-1] \, kpc^{-2}} \, ,
```
where N is the power-law index and $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $1 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `Σgas::Vector{<:SurfaceDensity}`: Values of the gas mass surface density, with units.
  - `log_output::Bool=true`: If the output will the $\log_{10}$ of the star formation area density, or the star formation area density itself (with units). If `log_output` = true, the implied unit is $\mathrm{M_\odot \, yr^{-1} \, kpc^{-2}}$

# Returns

  - The star formation area density.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
function kennicutt1998(Σgas::Vector{<:SurfaceDensity}; log_output::Bool=true)::Vector{<:Number}

    log10Σgas = @. log10(ustrip(u"Msun * pc^-2", Σgas))
    log10Σsfr = @. log10(a_KS98) + log10Σgas * N_KS98

    if log_output
        return log10Σsfr
    else
        return @. exp10(log10Σsfr ) * u"Msun * yr^-1 * kpc^-2"
    end

end

@doc raw"""
    invKennicutt1998(Σsfr::Vector{<:MassFlowDensity}; <keyword arguments>)::Vector{<:Number}

Inverse Kennicutt-Schmidt law, taken from a set of observations of nearby galaxies.

From Kennicutt (1998) (Section 4, Eq. 4), we have

```math
\Sigma_\mathrm{SFR} = a \left( \frac{\Sigma_\mathrm{gas}}{1 \, \mathrm{M_\odot \, pc^{-2}}} \right)^{\!N} \mathrm{M_\odot \, yr^{-1] \, kpc^{-2}} \, ,
```
where N is the power-law index and $a$ is $\Sigma_\mathrm{SFR}$ at the fiducial gas surface density of $1 \, \mathrm{M_\odot \, pc^{-2}}$.

# Arguments

  - `Σsfr::Vector{<:MassFlowDensity}`: Values of the star formation area density, with units.
  - `log_output::Bool=true`: If the output will the $\log_{10}$ of the gas mass surface density, or the gas mass surface density itself (with units). If `log_output` = true, the implied unit is $\mathrm{M_\odot \, pc^{-2}}$

# Returns

  - The gas mass surface density.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
function invKennicutt1998(Σsfr::Vector{<:MassFlowDensity}; log_output::Bool=true)::Vector{<:Number}

    log10Σsfr = @. log10(ustrip(u"Msun * yr^-1 * kpc^-2", Σsfr))
    log10Σgas = @. (log10Σsfr - log10(a_KS98)) / N_KS98

    if log_output
        return log10Σgas
    else
        return @. exp10(log10Σgas) * u"Msun * pc^-2"
    end

end

"""
    findClosestSnapshot(simulation_path::String, times::Vector{<:Unitful.Time})::Vector{Int}

Find the global index, in the context of the simulation, of the snapshot with a physical time closest to each of the ones given in `times`.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `times::Vector{<:Unitful.Time}`: Target physical times.

# Returns

  - The indices of the snapshots with physical times closest to `times`.
"""
function findClosestSnapshot(simulation_path::String, times::Vector{<:Unitful.Time})::Vector{Int}

    # Make a simulation table
    sim_table = makeSimulationTable(simulation_path)

    # Read the physical time associated to each snapshot
    snap_times = sim_table[!, :physical_times]

    # Allocate memory
    slices = similar(times, Int)

    # Find the closest snapshot to each of the `times`
    for (i, time) in pairs(times)

        slices[i] = argmin(abs.(snap_times .- time))

    end

    return slices

end

"""
    findClosestSnapshot(simulation_path::String, time::Unitful.Time)::Int

Find the global index, in the context of the simulation, of the snapshot with a physical time closest to `time`.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `time::Unitful.Time`: Target physical time.

# Returns

  - The index of the snapshot with a physical time closest to `time`.
"""
function findClosestSnapshot(simulation_path::String, time::Unitful.Time)::Int

    return findClosestSnapshot(simulation_path, [time])[1]

end

####################################################################################################
# Unitful.jl utilities
####################################################################################################

"""
List component units of a `Unitful.FreeUnits` object.
"""
unitList(::Unitful.FreeUnits{N,D,A}) where {N,D,A} = N

"""
List component units of a `Unitful.Quantity` object.
"""
unitList(x::Unitful.Quantity) = unitList(Unitful.units(x))

####################################################################################################
# Makie.jl and plotting utilities
####################################################################################################

"""
Extract the limits of the x axis, from a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will take the limits from the current axis object.
"""
function xlimits!(axis::Makie.Axis)::NTuple{2,Float32}
    reset_limits!(axis)
    return axis.xaxis.attributes.limits[]
end
xlimits!(plot::Makie.FigureAxisPlot)::NTuple{2,Float32} = xlimits!(plot.axis)
xlimits!(fig::Makie.Figure)::NTuple{2,Float32} = xlimits!(fig.current_axis.x)

"""
Extract the limits of the y axis, from a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will take the limits from the current axis object.
"""
function ylimits!(axis::Makie.Axis)::NTuple{2,Float32}
    reset_limits!(axis)
    return axis.yaxis.attributes.limits[]
end
ylimits!(plot::Makie.FigureAxisPlot)::NTuple{2,Float32} = ylimits!(plot.axis)
ylimits!(fig::Makie.Figure)::NTuple{2,Float32} = ylimits!(fig.current_axis.x)

"""
Extract the scale function of the x axis, from a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will take the scale from the current axis object.
"""
xscale(axis::Makie.Axis)::Function = axis.xaxis.attributes.scale[]
xscale(plot::Makie.FigureAxisPlot)::Function = xscale(plot.axis)
xscale(fig::Makie.Figure)::Function = xscale(fig.current_axis.x)

"""
Extract the scale function of the y axis, from a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will take the scale from the current axis object.
"""
yscale(axis::Makie.Axis)::Function = axis.yaxis.attributes.scale[]
yscale(plot::Makie.FigureAxisPlot)::Function = yscale(plot.axis)
yscale(fig::Makie.Figure)::Function = yscale(fig.current_axis.x)

"""
Extract all the data points in a [Makie](https://docs.makie.org/stable/) plot, axis, or figure. In the case of a figure, it will only take the data from the current axis object. It only works for scatter, line and scatterline plots.
"""
function pointData(axis::Makie.Axis)::Vector{Point{2,Float32}}

    series = copy(axis.scene.plots)

    filter!(x -> isa(x, Union{Scatter,Lines,ScatterLines}), series)

    !isempty(series) || return Point{2,Float32}[]

    return collect(Iterators.flatten(serie[1][] for serie in series))

end
pointData(plot::Makie.FigureAxisPlot)::Vector{Point{2,Float32}} = pointData(plot.axis)
pointData(fig::Makie.Figure)::Vector{Point{2,Float32}} = pointData(fig.current_axis.x)

"""
    absCoor(
        plot::Union{Makie.FigureAxisPlot,Makie.Axis,Makie.Figure},
        r_x::Real,
        r_y::Real,
    )::NTuple{2,Float64}

Compute the absolute x and y coordinates of a plot, from the relative ones.

# Arguments

  - `plot::Union{Makie.FigureAxisPlot,Makie.Axis,Makie.Figure}`: Plot, axis, or figure for which the absolute coordinates will be calculated. In the case of a figure, it will use the limits from the current axis object.
  - `r_x::Real`: Relative x coordinate.
  - `r_y::Real`: Relative y coordinate.

# Returns

  - A tuple with the absolute coordinates, (x, y).

# Examples

```julia-repl
julia> absCoor(lines(rand(100)), 0.5, 0.5)
(50.50000071525574, 0.48792968317866325)
```
"""
function absCoor(
    plot::Union{Makie.FigureAxisPlot,Makie.Axis,Makie.Figure},
    r_x::Real,
    r_y::Real,
)::NTuple{2,Float64}

    # Get the scaling functions
    x_scale = xscale(plot)
    y_scale = yscale(plot)

    # Get the limits of the axes
    x_limits = x_scale.(xlimits!(plot))
    y_limits = y_scale.(ylimits!(plot))

    # Compute the absolute coordinates
    a_x = Makie.inverse_transform(x_scale).(x_limits[1] + abs(r_x) * (x_limits[2] - x_limits[1]))
    a_y = Makie.inverse_transform(y_scale).(y_limits[1] + abs(r_y) * (y_limits[2] - y_limits[1]))

    return sign(r_x) * a_x, sign(r_y) * a_y

end

"""
    cleanPlot!(figure::Makie.Figure)::Nothing

Delete all the legends of a figure and empty all its axes.

# Arguments

  - `figure::Makie.Figure`: Figure to be cleaned.
"""
function cleanPlot!(figure::Makie.Figure)::Nothing

    # Compute the number of elements (axes and legends) in the figure
    n_elements = length(figure.content)
    i = 1

    for _ in 1:n_elements
        i += cleanPlot!(figure.content[i])
    end

    return nothing

end

"""
    cleanPlot!(ax::Makie.Axis)::Bool

Empty an axis.

# Arguments

  - `ax::Makie.Axis`: Axis to be emptied.

# Returns

  - Flag to indicate that an axis has been emptied.
"""
function cleanPlot!(ax::Makie.Axis)::Bool

    empty!(ax)

    return true

end

"""
    cleanPlot!(legend::Union{Makie.Legend,Makie.Colorbar})::Bool

Delete a legend or colorbar.

# Arguments

  - `legend::Union{Makie.Legend,Makie.Colorbar}`: Legend or colorbar to be deleted.

# Returns

  - Flag to indicate that a legend or colorbar has been deleted.
"""
function cleanPlot!(legend::Union{Makie.Legend,Makie.Colorbar})::Bool

    delete!(legend)

    return false

end

"""
Default function to end `cleanPlot!` recursion if an unknown type is encountered.
"""
cleanPlot!(default) = error("cleanPlot!: I cannot clean elements of type $(typeof(default))")

"""
    getUnitLabel(factor::Int, unit::Unitful.Units; <keyword arguments>)::AbstractString

Construct the unit part of an axis label.

# Arguments

  - `factor::Int`: Exponential factor to scale down the units. If different from 0, a term of the form 10^`factor` will be added to the label.
  - `unit::Unitful.Units`: Unit of the axis.
  - `latex::Bool=true`: If the output will be a `LaTeXString`, or a plain `String`.

# Returns

  - The `LaTeXString` or `String`: "10^`factor` `unit`". The `factor` term only appears if `factor` != 0, the unit term only appears if `unit` != `Unitful.NoUnits`.
"""
function getUnitLabel(factor::Int, unit::Unitful.Units; latex::Bool=true)::AbstractString

    if latex

        # Replace special characters for their LaTeX counterpart
        str_unit = replace(
            string(unit),
            "M⊙" => raw"M_{\odot}",
            r"\^(-?[0-9]+)" => s"^{\1}",
            " " => raw"\,",
        )

        if iszero(factor)
            if isempty(str_unit)
                out_str = ""
            else
                out_str = L"\mathrm{%$str_unit}"
            end
        else
            if isempty(str_unit)
                out_str = L"10^{%$factor}"
            else
                out_str = L"10^{%$factor} \, \mathrm{%$str_unit}"
            end
        end

    else

        str_unit = string(unit)

        if iszero(factor)
            if isempty(str_unit)
                out_str = ""
            else
                out_str = str_unit
            end
        else
            if isempty(str_unit)
                out_str = "10^$factor"
            else
                out_str = "10^$factor $str_unit"
            end
        end

    end

    return out_str

end

"""
    getLabel(
        label::AbstractString,
        factor::Int,
        unit::Unitful.Units;
        <keyword arguments>
    )::AbstractString

Construct an axis label.

# Arguments

  - `label::AbstractString`: Variable name.
  - `factor::Int`: Exponential factor to scale down the units. If different from 0, a term of the form 10^`factor` will be added to the label.
  - `unit::Unitful.Units`: Unit of the axis.
  - `latex::Bool=true`: If the output will be a `LaTeXString`, or a plain `String`.

# Returns

  - The `LaTeXString` or `String`: "`label` [10^`factor` `unit`]". If `label` is "", an empty string is returned. The `factor` term only appears if `factor` != 0, the unit term only appears if `unit` != `Unitful.NoUnits`, and the brackets only appears if there are a factor and/or a unit term.
"""
function getLabel(
    label::AbstractString,
    factor::Int,
    unit::Unitful.Units;
    latex::Bool=true,
)::AbstractString

    !isempty(label) || return ""

    unit_label = getUnitLabel(factor, unit; latex)

    if isempty(unit_label)
        return latex ? L"%$label" : label
    end

    return latex ? L"%$label [%$unit_label]" : "$label [$unit_label]"

end

"""
    barPlotLabelFormater(x::Number)::LaTeXString

Format a number to be a barplot label.

For values between 0 and 0.01 the label will be "< 0.01", otherwise it will be the value itself with 2 digits.

# Arguments

  - `x::Number`: Value to be formated.

# Returns

  - The bar label.
"""
function barPlotLabelFormater(x::Number)::LaTeXString

    if 0 < x < 0.01
        return L"< \, 0.01"
    end

    return latexstring(round(x; digits=2))

end

"""
    barPlotLabelFormater(x::LaTeXString)::LaTeXString

Format a number to be a barplot label.

Method for compatibility with the barplot! function of [Makie](https://docs.makie.org/stable/).

# Arguments

  - `x::Number`: Value to be formated.

# Returns

  - The bar label.
"""
barPlotLabelFormater(x::LaTeXString)::LaTeXString = x

"""
    formatError(q_mean::Number, q_error::Number)::NTuple{2,<:Number}

Nicely format a magnitude with uncertainty.

It follows the traditional rules for error presentation: the error has only one significant digit, unless such digit is a one, in which case two significant digits are used. The mean will have as many digits as to match the last significant position of the error. An error equal to 0 will leave the mean unchanged.

# Arguments

  - `q_mean::Number`: Mean value.
  - `q_error::Number`: Error value. It must be positive.

# Returns

  - A tuple with the formatted mean and error values.

# Examples

```julia-repl
julia> formatError(69.42069, 0.038796)
(69.42, 0.04)

julia> formatError(69.42069, 0.018796)
(69.421, 0.019)

julia> formatError(15.42, 0.00004)
(15.42, 4.0e-5)

julia> formatError(69.42069, 0.0)
(69.42069, 0.0)

julia> formatError(69.42069, 93.4)
(70.0, 90.0)

julia> formatError(69.42069, 123.4)
(70.0, 120.0)

julia> formatError(15.42069, 16.4)
(15.0, 16.0)
```
"""
function formatError(q_mean::Number, q_error::Number)::NTuple{2,<:Number}

    # Positive error check
    (
        q_error >= zero(q_error) ||
        throw(ArgumentError("formatError: `q_error` must be positive, but I got \
        `q_error` = $(q_error)"))
    )

    # Use the values without units
    mean = ustrip(q_mean)
    error = ustrip(q_error)

    if iszero(error)

        round_mean = mean
        round_error = error

    else

        sigdigit_pos = abs(log10(error))

        if error < 1.0
            first_digit = trunc(error * 10.0^(floor(sigdigit_pos) + 1.0))
            extra = first_digit == 1.0 ? 1 : 0
            digits = ceil(Int, sigdigit_pos) + extra
            round_mean = round(mean; digits)
            round_error = round(error, sigdigits=1 + extra)
        else
            first_digit = trunc(error * 10.0^(-floor(sigdigit_pos)))
            extra = first_digit == 1.0 ? 2 : 1
            sigdigits = ceil(Int, log10(abs(mean))) - ceil(Int, sigdigit_pos) + extra
            round_mean = round(mean; sigdigits)
            round_error = round(error, sigdigits=extra)
        end

    end

    return round_mean * unit(q_mean), round_error * unit(q_error)

end

"""
    plotParams(quantity::Symbol)::PlotParams

Select the plotting parameters for a given `quantity`.

# Arguments

  - `quantity::Symbol`: The options are:

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
      + `:generic_mass`                -> Parameters for plots with several diferent masses.
      + `:stellar_number`              -> Number of stellar particles.
      + `:gas_number`                  -> Number of gas cells.
      + `:dm_number`                   -> Number of dark matter particles.
      + `:bh_number`                   -> Number of black hole particles.
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
      + `:generic_fraction`            -> Parameters for plots with several diferent fraction.
      + `:gas_mass_density`            -> Gas mass density.
      + `:hydrogen_mass_density`       -> Hydrogen mass density.
      + `:gas_number_density`          -> Gas number density.
      + `:molecular_number_density`    -> Molecular hydrogen number density.
      + `:br_molecular_number_density` -> Molecular hydrogen number density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_number_density`       -> Atomic hydrogen number density.
      + `:ionized_number_density`      -> Ionized hydrogen number density.
      + `:neutral_number_density`      -> Neutral hydrogen number density.
      + `:stellar_area_density`        -> Stellar area mass density, for a radius of `DISK_R`.
      + `:gas_area_density`            -> Gas mass surface density, for a radius of `DISK_R`.
      + `:molecular_area_density`      -> Molecular mass surface density, for a radius of `DISK_R`.
      + `:br_molecular_area_density`   -> Molecular mass surface density, for a radius of `DISK_R`, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_area_density`         -> Atomic hydrogen area mass density, for a radius of `DISK_R`.
      + `:ionized_area_density`        -> Ionized hydrogen area mass density, for a radius of `DISK_R`.
      + `:neutral_area_density`        -> Neutral mass surface density, for a radius of `DISK_R`.
      + `:sfr_area_density`            -> Star formation rate area density, for the last `AGE_RESOLUTION` and a radius of `DISK_R`.
      + `:generic_area_density`        -> Parameters for plots with several diferent area densities.
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
      + `:mass_accretion`              -> Gas accretion rate. Positive values mean gas infall into the virial radius ``R_{200}``, and negative values mean outflow.
      + `:stellar_specific_am`         -> Norm of the stellar specific angular momentum.
      + `:gas_specific_am`             -> Norm of the gas specific angular momentum.
      + `:dm_specific_am`              -> Norm of the dark matter specific angular momentum.
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
      + `:scale_factor`                -> Scale factor.
      + `:redshift`                    -> Redshift.
      + `:physical_time`               -> Physical time since the Big Bang.
      + `:lookback_time`               -> Physical time left to reach the last snapshot.
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

  - A [`PlotParams`](@ref) object, with entries:

      + `request::Dict{Symbol,Vector{String}}` -> Data request for [`readSnapshot`](@ref).
      + `var_name::AbstractString`             -> Name of the quantity for the plot axis.
      + `exp_factor::Int`                      -> Numerical exponent to scale down the axis.
      + `unit::Unitful.Units`                  -> Target unit for the axis.
      + `axis_label::AbstractString`           -> Label for the axis.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)
"""
function plotParams(quantity::Symbol)::PlotParams

    if quantity == :stellar_mass

        plot_params = PlotParams(;
            request    = Dict(:stars => ["MASS", "POS "]),
            var_name   = L"M_\star",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :gas_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name   = L"M_\mathrm{gas}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :hydrogen_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name   = L"M_\mathrm{H}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :dm_mass

        plot_params = PlotParams(;
            request    = Dict(:halo => ["MASS", "POS "]),
            var_name   = L"M_\mathrm{DM}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :bh_mass

        plot_params = PlotParams(;
            request  = Dict(:black_hole => ["MASS", "POS "]),
            var_name = L"M_\mathrm{BH}",
            unit     = u"Msun",
        )

    elseif quantity == :molecular_mass

        plot_params = PlotParams(;
            request    = Dict(
                :gas => [
                    "MASS",
                    "POS ",
                    "FRAC",
                    "RHO ",
                    "TAUS",
                    "PRES",
                    "NH  ",
                    "NHP ",
                ],
            ),
            var_name   = L"M_\mathrm{H_2}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :br_molecular_mass

        plot_params = PlotParams(;
            request    = Dict(
                :gas => ["MASS", "PRES", "RHO ", "NH  ", "NHP "],
            ),
            var_name   = L"M_\mathrm{H_2^{BR}}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :atomic_mass

        plot_params = PlotParams(;
            request    = Dict(
                :gas => [
                    "MASS",
                    "POS ",
                    "FRAC",
                    "NH  ",
                    "NHP ",
                    "RHO ",
                    "TAUS",
                    "PRES",
                ],
            ),
            var_name   = L"M_\mathrm{HI}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :ionized_mass

        plot_params = PlotParams(;
            request    = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "TAUS"],
            ),
            var_name   = L"M_\mathrm{HII}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :neutral_mass

        plot_params = PlotParams(;
            request    = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "TAUS"],
            ),
            var_name   = L"M_\mathrm{HI + H_2}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :stellar_gas_mass

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "FRAC", "RHO "]),
            var_name   = L"M_\star^\mathrm{gas}",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :generic_mass

        plot_params = PlotParams(;
            request    = Dict(
                :stars   => ["MASS", "POS "],
                :gas     => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO ", "TAUS"],
                :dm_mass => ["MASS", "POS "],
                :bh_mass => ["MASS", "POS "],
            ),
            var_name   = L"M",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :stellar_number

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS "]),
            var_name = L"\mathrm{Number \,\, of \,\, stellar \,\, particles}",
        )

    elseif quantity == :gas_number

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS "]),
            var_name = L"\mathrm{Number \,\, of \,\, gas \,\, cells}",
        )

    elseif quantity == :dm_number

        plot_params = PlotParams(;
            request  = Dict(:halo => ["MASS", "POS "]),
            var_name = L"\mathrm{Number \,\, of \,\, DM \,\, particles}",
        )

    elseif quantity == :bh_number

        plot_params = PlotParams(;
            request  = Dict(:black_hole => ["MASS", "POS "]),
            var_name = L"\mathrm{Number \,\, of \,\, BH \,\, particles}",
        )

    elseif quantity == :molecular_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["RHO ", "MASS", "POS ", "FRAC", "TAUS", "PRES", "NH  ", "NHP "],
            ),
            var_name = L"f_\mathrm{H_2}",
        )

    elseif quantity == :br_molecular_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "PRES", "NH  ", "NHP "],
            ),
            var_name = L"f_\mathrm{H_2}^\mathrm{BR}",
        )

    elseif quantity == :atomic_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "TAUS", "PRES"],
            ),
            var_name = L"f_\mathrm{HI}",
        )

    elseif quantity == :ionized_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "TAUS"],
            ),
            var_name = L"f_\mathrm{HII}",
        )

    elseif quantity == :neutral_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "TAUS"],
            ),
            var_name = L"f_\mathrm{H_I + H_2}",
        )

    elseif quantity == :molecular_neutral_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO ", "TAUS"],
            ),
            var_name = L"f_\mathrm{H_2} \, / f_\mathrm{n}",
        )

    elseif quantity == :ionized_neutral_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO ", "TAUS"],
            ),
            var_name = L"f_\mathrm{HII} \, / f_\mathrm{n}",
        )

    elseif quantity == :stellar_gas_fraction

        plot_params = PlotParams(;
            request  = Dict(:gas => ["RHO ", "FRAC"]),
            var_name = L"f_\star",
        )

    elseif quantity == :mol_eq_quotient

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["ETAD", "FRAC", "RHO ", "PARZ"],
            ),
            var_name = L"\log_{10} \, \mathrm{LS^{H_2} / RS^{H_2}}",
        )

    elseif quantity == :ion_eq_quotient

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["ETAI", "PARR", "FRAC", "RHO "],
            ),
            var_name = L"\log_{10} \, \mathrm{LS^{HII} / RS^{HII}}",
        )

    elseif quantity == :generic_fraction

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["RHO ", "MASS", "POS ", "FRAC", "NH  ", "NHP ", "PRES", "RHO ", "TAUS"],
            ),
            var_name = L"f",
        )

    elseif quantity == :gas_mass_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name = L"\rho_\mathrm{gas}",
            unit     = u"Msun*kpc^-3",
        )

    elseif quantity == :hydrogen_mass_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name = L"\rho_\mathrm{H}",
            unit     = u"Msun*kpc^-3",
        )

    elseif quantity == :gas_number_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name = L"n_\mathrm{gas}",
            unit     = u"cm^-3",
        )

    elseif quantity == :molecular_number_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => [
                    "MASS",
                    "POS ",
                    "FRAC",
                    "RHO ",
                    "TAUS",
                    "PRES",
                    "NH  ",
                    "NHP ",
                ],
            ),
            var_name = L"n_\mathrm{H_2}",
            unit     = u"cm^-3",
        )

    elseif quantity == :br_molecular_number_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "PRES", "RHO ", "NH  ", "NHP "],
            ),
            var_name = L"n_\mathrm{H_2}^{BR}",
            unit     = u"cm^-3",
        )

    elseif quantity == :atomic_number_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => [
                    "MASS",
                    "POS ",
                    "FRAC",
                    "NH  ",
                    "NHP ",
                    "RHO ",
                    "TAUS",
                    "PRES",
                ],
            ),
            var_name = L"n_\mathrm{HI}",
            unit     = u"cm^-3",
        )

    elseif quantity == :ionized_number_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "TAUS"],
            ),
            var_name = L"n_\mathrm{HII}",
            unit     = u"cm^-3",
        )

    elseif quantity == :neutral_number_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "TAUS"],
            ),
            var_name = L"n_\mathrm{HI + H_2}",
            unit     = u"cm^-3",
        )

    elseif quantity == :stellar_area_density

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS "]),
            var_name = L"\Sigma_\star",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :gas_area_density

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "RHO "]),
            var_name = L"\Sigma_\mathrm{gas}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :molecular_area_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => [
                    "MASS",
                    "POS ",
                    "FRAC",
                    "RHO ",
                    "TAUS",
                    "PRES",
                    "NH  ",
                    "NHP ",
                ],
            ),
            var_name = L"\Sigma_\mathrm{H_2}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :br_molecular_area_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "PRES", "RHO ", "NH  ", "NHP "],
            ),
            var_name = L"\Sigma_\mathrm{H_2}^\mathrm{BR}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :atomic_area_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => [
                    "MASS",
                    "POS ",
                    "FRAC",
                    "NH  ",
                    "NHP ",
                    "RHO ",
                    "TAUS",
                    "PRES",
                ],
            ),
            var_name = L"\Sigma_\mathrm{HI}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :ionized_area_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "TAUS"],
            ),
            var_name = L"\Sigma_\mathrm{HII}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :neutral_area_density

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "TAUS"],
            ),
            var_name = L"\Sigma_\mathrm{HI + H_2}",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :sfr_area_density

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "GAGE"]),
            var_name = L"\Sigma_\mathrm{SFR}",
            unit     = u"Msun*yr^-1*kpc^-2",
        )

    elseif quantity == :generic_area_density

        plot_params = PlotParams(;
            request  = Dict(
                :stars => ["MASS", "POS ", "GAGE"],
                :gas => [
                    "MASS", "POS ", "FRAC", "NH  ", "NHP ", "RHO ", "PRES", "TAUS",
                ],
            ),
            var_name = L"\Sigma",
            unit     = u"Msun*pc^-2",
        )

    elseif quantity == :gas_td

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "SFR "]),
            var_name = L"t_d^\mathrm{gas}",
            unit     = u"Gyr",
        )

    elseif quantity == :molecular_td

        plot_params = PlotParams(;
            request  = Dict(
                :gas => [
                    "MASS",
                    "POS ",
                    "FRAC",
                    "RHO ",
                    "TAUS",
                    "PRES",
                    "NH  ",
                    "NHP ",
                    "SFR ",
                ],
            ),
            var_name = L"t_d^\mathrm{H_2}",
            unit     = u"Gyr",
        )

    elseif quantity == :br_molecular_td

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "PRES", "RHO ", "NH  ", "NHP ", "SFR "],
            ),
            var_name = L"t_d^\mathrm{H_2^\mathrm{BR}}",
            unit     = u"Gyr",
        )

    elseif quantity == :atomic_td

        plot_params = PlotParams(;
            request  = Dict(
                :gas => [
                    "MASS",
                    "POS ",
                    "FRAC",
                    "NH  ",
                    "NHP ",
                    "RHO ",
                    "TAUS",
                    "PRES",
                    "SFR ",
                ],
            ),
            var_name = L"t_d^\mathrm{HI}",
            unit     = u"Gyr",
        )

    elseif quantity == :ionized_td

        plot_params = PlotParams(;
            request  = Dict(
                :gas => [
                    "MASS",
                    "POS ",
                    "FRAC",
                    "NH  ",
                    "NHP ",
                    "RHO ",
                    "TAUS",
                    "SFR ",
                ],
            ),
            var_name = L"t_d^\mathrm{HII}",
            unit     = u"Gyr",
        )

    elseif quantity == :neutral_td

        plot_params = PlotParams(;
            request  = Dict(
                :gas => [
                    "MASS",
                    "POS ",
                    "FRAC",
                    "NH  ",
                    "NHP ",
                    "RHO ",
                    "TAUS",
                    "SFR ",
                ],
            ),
            var_name = L"t_d^\mathrm{HI + H_2}",
            unit     = u"Gyr",
        )

    elseif quantity == :gas_metallicity

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "GMET", "GZ  "]),
            var_name = L"Z_\mathrm{gas} \, [\mathrm{Z_\odot}]",
        )

    elseif quantity == :stellar_metallicity

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "GME2", "GZ2 "]),
            var_name = L"Z_\star \, [\mathrm{Z_\odot}]",
        )

    elseif quantity ∈ GAS_ABUNDANCE

        element_string = first(split(string(quantity), "_"))

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "GMET"]),
            axis_label = L"12 + \log_{10}(\mathrm{%$element_string} \, / \, \mathrm{H})",
        )

    elseif quantity ∈ STELLAR_ABUNDANCE

        element_string = first(split(string(quantity), "_"))

        plot_params = PlotParams(;
            request    = Dict(:stars => ["MASS", "POS ", "GME2"]),
            axis_label = L"12 + \log_{10}(\mathrm{%$element_string} \, / \, \mathrm{H})",
        )

    elseif quantity == :stellar_radial_distance

        plot_params = PlotParams(;
            request  = Dict(:stars => ["POS "]),
            var_name = L"r",
            unit     = u"kpc",
        )

    elseif quantity == :gas_radial_distance

        plot_params = PlotParams(;
            request  = Dict(:gas => ["POS "]),
            var_name = L"r",
            unit     = u"kpc",
        )

    elseif quantity == :dm_radial_distance

        plot_params = PlotParams(;
            request  = Dict(:halo => ["POS "]),
            var_name = L"r",
            unit     = u"kpc",
        )

    elseif quantity == :stellar_xy_distance

        plot_params = PlotParams(;
            request  = Dict(:stars => ["POS "]),
            var_name = L"r_{xy}",
            unit     = u"kpc",
        )

    elseif quantity == :gas_xy_distance

        plot_params = PlotParams(;
            request  = Dict(:gas => ["POS "]),
            var_name = L"r_{xy}",
            unit     = u"kpc",
        )

    elseif quantity == :dm_xy_distance

        plot_params = PlotParams(;
            request  = Dict(:halo => ["POS "]),
            var_name = L"r_{xy}",
            unit     = u"kpc",
        )

    elseif quantity == :gas_sfr

        plot_params = PlotParams(;
            request  = Dict(:gas => ["SFR "]),
            var_name = L"\mathrm{SFR_{gas}}",
            unit     = u"Msun*yr^-1",
        )

    elseif quantity == :mass_accretion

        plot_params = PlotParams(;
            request  = Dict(
                :gas         => ["ID  ", "MASS"],
                :stars       => ["ID  ", "MASS"],
                :black_hole  => ["ID  ", "MASS"],
                :group       => ["G_R_Crit200", "G_M_Crit200"],
                :tracer      => ["PAID", "TRID"],
            ),
            var_name = "Mass accretion",
            unit     = u"Msun*yr^-1",
        )

    elseif quantity == :stellar_specific_am

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "VEL"]),
            var_name = L"j_\star",
            unit     = u"kpc^2*s^-1",
        )

    elseif quantity == :gas_specific_am

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "POS ", "VEL"]),
            var_name = L"j_\mathrm{gas}",
            unit     = u"kpc^2*s^-1",
        )

    elseif quantity == :dm_specific_am

        plot_params = PlotParams(;
            request  = Dict(:halo => ["MASS", "POS ", "VEL"]),
            var_name = L"j_\mathrm{DM}",
            unit     = u"kpc^2*s^-1",
        )

    elseif quantity == :stellar_circularity

        # `daBandProfile` expects that the first element in the request dictionary is for the stars
        plot_params = PlotParams(;
            request  = Dict(
                :stars      => ["MASS", "POS ", "VEL "],
                :gas        => ["MASS", "POS "],
                :halo       => ["MASS", "POS "],
                :black_hole => ["MASS", "POS "],
            ),
            var_name = L"\epsilon",
        )

    elseif quantity == :stellar_vcirc

        # `daBandProfile` expects that the first element in the request dictionary is for the stars
        plot_params = PlotParams(;
            request  = Dict(
                :stars      => ["MASS", "POS "],
                :gas        => ["MASS", "POS "],
                :halo       => ["MASS", "POS "],
                :black_hole => ["MASS", "POS "],
            ),
            var_name = L"v_\mathrm{circ}",
            unit     = u"km*s^-1",
        )

    elseif quantity == :stellar_vradial

        plot_params = PlotParams(;
            request  = Dict(:stars => ["POS ", "VEL "]),
            var_name = L"v_r",
            unit     = u"km*s^-1",
        )

    elseif quantity == :stellar_vtangential

        plot_params = PlotParams(;
            request  = Dict(:stars => ["POS ", "VEL "]),
            var_name = L"v_\theta",
            unit     = u"km*s^-1",
        )

    elseif quantity == :stellar_vzstar

        plot_params = PlotParams(;
            request  = Dict(:stars => ["POS ", "VEL "]),
            var_name = L"v_z \,\, \mathrm{sign}(z)",
            unit     = u"km*s^-1",
        )

    elseif quantity == :stellar_age

        plot_params = PlotParams(;
            request  = Dict(:stars => ["GAGE"]),
            var_name = L"\mathrm{Stellar \,\, age}",
            unit     = u"Gyr",
        )

    elseif quantity == :sfr

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "GAGE"]),
            var_name = L"SFR",
            unit     = u"Msun*yr^-1",
        )

    elseif quantity == :ssfr

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "GAGE"]),
            var_name = L"sSFR",
            unit     = u"yr^-1",
        )

    elseif quantity == :observational_sfr

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "GAGE"]),
            var_name = L"SFR",
            unit     = u"Msun*yr^-1",
        )

    elseif quantity == :observational_ssfr

        plot_params = PlotParams(;
            request  = Dict(:stars => ["MASS", "POS ", "GAGE"]),
            var_name = L"sSFR",
            unit     = u"yr^-1",
        )

    elseif quantity == :stellar_eff

        plot_params = PlotParams(;
            request  = Dict(:stars => ["GMAS", "GSFR", "RHOC"]),
            var_name = L"\epsilon_\mathrm{ff}^\star",
        )

    elseif quantity == :gas_eff

        plot_params = PlotParams(;
            request  = Dict(:gas => ["MASS", "SFR ", "RHO "]),
            var_name = L"\epsilon_\mathrm{ff}^\mathrm{gas}",
        )

    elseif quantity == :molecular_eff

        plot_params = PlotParams(;
            request  = Dict(
                :gas => [
                    "MASS",
                    "POS ",
                    "FRAC",
                    "RHO ",
                    "TAUS",
                    "PRES",
                    "NH  ",
                    "NHP ",
                    "SFR ",
                ],
            ),
            var_name = L"\epsilon_\mathrm{ff}^\mathrm{H_2}",
        )

    elseif quantity == :br_molecular_eff

        plot_params = PlotParams(;
            request  = Dict(
                :gas => ["MASS", "PRES", "RHO ", "NH  ", "NHP ", "SFR "],
            ),
            var_name = L"\epsilon_\mathrm{ff}^\mathrm{H_2^\mathrm{BR}}",
        )

    elseif quantity == :atomic_eff

        plot_params = PlotParams(;
            request  = Dict(
                :gas => [
                    "MASS",
                    "POS ",
                    "FRAC",
                    "NH  ",
                    "NHP ",
                    "RHO ",
                    "TAUS",
                    "PRES",
                    "SFR ",
                ],
            ),
            var_name = L"\epsilon_\mathrm{ff}^\mathrm{HI}",
        )

    elseif quantity == :ionized_eff

        plot_params = PlotParams(;
            request  = Dict(
                :gas => [
                    "MASS",
                    "POS ",
                    "FRAC",
                    "NH  ",
                    "NHP ",
                    "RHO ",
                    "TAUS",
                    "SFR ",
                ],
            ),
            var_name = L"\epsilon_\mathrm{ff}^\mathrm{HII}",
        )

    elseif quantity == :neutral_eff

        plot_params = PlotParams(;
            request  = Dict(
                :gas => [
                    "MASS",
                    "POS ",
                    "FRAC",
                    "NH  ",
                    "NHP ",
                    "RHO ",
                    "TAUS",
                    "SFR ",
                ],
            ),
            var_name = L"\epsilon_\mathrm{ff}^\mathrm{HI + H_2}",
        )

    elseif quantity == :temperature

        plot_params = PlotParams(;
            request    = Dict(:gas => ["MASS", "POS ", "TEMP"]),
            axis_label = L"\log_{10} \, T \, [\mathrm{K}]",
        )

    elseif quantity == :pressure

        plot_params = PlotParams(;
            request    = Dict(:gas => ["PRES"]),
            var_name   = L"P",
            unit       = u"Pa",
        )

    elseif quantity == :scale_factor

        plot_params = PlotParams(;
            var_name = L"a",
        )

    elseif quantity == :redshift

        plot_params = PlotParams(;
            var_name = L"z",
        )

    elseif quantity == :physical_time

        plot_params = PlotParams(;
            var_name = L"t",
            unit = u"Gyr",
        )

    elseif quantity == :lookback_time

        plot_params = PlotParams(;
            var_name = L"\mathrm{Lookback \,\, time}",
            unit     = u"Gyr",
        )

    elseif quantity == :time_step

        plot_params = PlotParams(;
            var_name = L"\mathrm{Number \,\, of \,\, time \,\, steps}",
        )

    elseif quantity == :clock_time_s

        plot_params = PlotParams(;
            var_name = L"\mathrm{Wallclock \,\, time}",
            unit     = u"s",
        )

    elseif quantity == :clock_time_percent

        plot_params = PlotParams(;
            axis_label = L"\mathrm{Wallclock \,\, time \, [\%]}",
        )

    elseif quantity == :tot_clock_time_s

        plot_params = PlotParams(;
            var_name = L"\mathrm{Cumulative \,\, wallclock \,\, time}",
            unit     = u"s",
        )

    elseif quantity == :tot_clock_time_percent

        plot_params = PlotParams(;
            axis_label = L"\mathrm{Cumulative \,\, wallclock \,\, time \, [\%]}",
        )

    elseif quantity == :ode_gas_it

        plot_params = PlotParams(;
            request  = Dict(:gas => ["ODIT"]),
            var_name = L"\mathrm{Integration\,\, time}",
            unit     = u"Myr",
        )

    elseif quantity == :ode_gas_accu_it

        plot_params = PlotParams(;
            request  = Dict(:gas => ["ACIT"]),
            var_name = L"\mathrm{Accumulated \,\, integration \,\, time}",
            unit     = u"Gyr",
        )

    elseif quantity == :ode_gas_tau_s

        plot_params = PlotParams(;
            request  = Dict(:gas => ["TAUS"]),
            var_name = L"\tau_\mathrm{S}",
            unit     = u"Myr",
        )

    elseif quantity == :ode_gas_eta_d

        plot_params = PlotParams(;
            request  = Dict(:gas => ["ETAD"]),
            var_name = L"\eta_\mathrm{diss}",
        )

    elseif quantity == :ode_gas_eta_i

        plot_params = PlotParams(;
            request  = Dict(:gas => ["ETAI"]),
            var_name = L"\eta_\mathrm{ion}",
        )

    elseif quantity == :ode_gas_r

        plot_params = PlotParams(;
            request  = Dict(:gas => ["PARR"]),
            var_name = L"R",
        )

    elseif quantity == :ode_gas_cold_mf

        plot_params = PlotParams(;
            request  = Dict(:gas => ["COLF"]),
            var_name = L"c_f",
        )

    elseif quantity == :ode_stellar_it

        plot_params = PlotParams(;
            request  = Dict(:stars => ["ODIT"]),
            var_name = L"\mathrm{it}^\star",
            unit     = u"Myr",
        )

    elseif quantity == :ode_stellar_accu_it

        plot_params = PlotParams(;
            request  = Dict(:stars => ["ACIT"]),
            var_name = L"\mathrm{ait}^\star",
            unit     = u"Gyr",
        )

    elseif quantity == :ode_stellar_tau_s

        plot_params = PlotParams(;
            request  = Dict(:stars => ["TAUS"]),
            var_name = L"\tau_\mathrm{S}^\star",
            unit     = u"Myr",
        )

    elseif quantity == :ode_stellar_eta_d

        plot_params = PlotParams(;
            request  = Dict(:stars => ["ETAD"]),
            var_name = L"\eta_\mathrm{diss}^\star",
        )

    elseif quantity == :ode_stellar_eta_i

        plot_params = PlotParams(;
            request  = Dict(:stars => ["ETAI"]),
            var_name = L"\eta_\mathrm{ion}^\star",
        )

    elseif quantity == :ode_stellar_r

        plot_params = PlotParams(;
            request  = Dict(:stars => ["PARR"]),
            var_name = L"R^\star",
        )

    elseif quantity == :ode_stellar_cold_mf

        plot_params = PlotParams(;
            request  = Dict(:stars => ["COLF"]),
            var_name = L"c_f^\star",
        )

    elseif quantity == :ode_stellar_gas_rho

        plot_params = PlotParams(;
            request  = Dict(:stars => ["RHOC"]),
            var_name = L"\rho_\mathrm{gas}^\star",
            unit     = u"Msun*kpc^-3",
        )

    elseif quantity == :ode_stellar_gas_Z

        plot_params = PlotParams(;
            request  = Dict(:stars => ["PARZ"]),
            var_name = L"Z_\mathrm{gas}^\star \, [\mathrm{Z_\odot}]",
        )

    elseif quantity == :ode_stellar_gas_mass

        plot_params = PlotParams(;
            request    = Dict(:stars => ["GMAS"]),
            var_name   = L"M_\mathrm{gas}^\star",
            exp_factor = 10,
            unit       = u"Msun",
        )

    elseif quantity == :ode_stellar_gas_sfr

        plot_params = PlotParams(;
            request  = Dict(:stars => ["GSFR"]),
            var_name = L"\mathrm{SFR}_\mathrm{gas}^\star",
            unit     = u"Msun*yr^-1",
        )

    elseif quantity == :ode_stellar_gas_P

        plot_params = PlotParams(;
            request    = Dict(:stars => ["GPRE"]),
            var_name   = L"P^\star",
            unit       = u"Pa",
        )

    else

        throw(ArgumentError("plotParams: I don't recognize the quantity :$(quantity)"))

    end

    return plot_params

end
