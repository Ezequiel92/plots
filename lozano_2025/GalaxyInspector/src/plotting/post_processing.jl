####################################################################################################
# Post-processing functions
####################################################################################################
#
# A post-processing function must take a Makie figure, add something to it, and return how to label
# the additions (or `nothing` when no new labels should be drawn).
#
# Expected signature:
#
#   post_processing(figure, args...; kwargs...) -> ([marker, ...], [label, ...]) or `nothing`
#
# where:
#
#   - figure::Makie.Figure
#   - marker::LegendElement
#   - label::AbstractString
#
####################################################################################################

"""
    ppVerticalFlags!(
        figure::Makie.Figure,
        positions::Vector{<:Real};
        <keyword arguments>
    )::Nothing

Draw vertical lines.

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `positions::Vector{<:Real}`: The x coordinates of the lines.
  - `colors::Vector{<:ColorType}=[:red]`: Colors of the lines.
  - `line_styles::Vector{<:LineStyleType}=[nothing]`: Styles of the lines. `nothing` will produce a solid line.
"""
function ppVerticalFlags!(
    figure::Makie.Figure,
    positions::Vector{<:Real};
    colors::Vector{<:ColorType}=[:red],
    line_styles::Vector{<:LineStyleType}=[nothing],
)::Nothing

    # Filter out values outside the range of the original plot
    rangeCut!(positions, xlimits!(figure); keep_edges=false)

    if isempty(positions)

        # Don't draw anything if all the flags are outside the original plot window
        !logging[] || @warn("ppVerticalFlags!: All vertical lines lie outside the plot range")

    else

        # Draw the vertical lines
        for (i, position) in pairs(positions)
            color = ring(colors, i)
            linestyle = ring(line_styles, i)
            vlines!(figure.current_axis.x, position; color, linestyle)
        end

    end

    return nothing

end

"""
    ppFillBelowLine!(
        figure::Makie.Figure;
        <keyword arguments>
    )::Nothing

Fill the space below a line plot, up to a `lower_limit`, with a solid `color`.

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `lower_limit::Number=0.0`: Lower bound.
  - `color::ColorType=:red`: Color.
  - `alpha::Float64=0.3`: Level of transparency. 1.0 is completely opaque and 0.0 is completely transparent.
"""
function ppFillBelowLine!(
    figure::Makie.Figure;
    lower_limit::Number=0.0,
    color::ColorType=:red,
    alpha::Float64=0.3,
)::Nothing

    # Read the data points in the plot
    points = pointData(figure)

    if isempty(points)
        !logging[] || @warn("ppFillBelowLine!: There are no points in the figure")
        return nothing
    end

    xs = Vector{Float64}(undef, length(points))
    ys = Vector{Float64}(undef, length(points))

    for (i, point) in pairs(points)
        xs[i] = point[1]
        ys[i] = point[2]
    end

    ys_low = fill(lower_limit, length(xs))

    band!(figure.current_axis.x, xs, ys_low, ys; color, alpha)

    return nothing

end

"""
    ppHorizontalFlags!(
        figure::Makie.Figure,
        positions::Vector{<:Real};
        <keyword arguments>
    )::Nothing

Draw horizontal lines.

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `positions::Vector{<:Real}`: The y coordinates of the lines.
  - `colors::Vector{<:ColorType}=[:red]`: Colors of the lines.
  - `line_styles::Vector{<:LineStyleType}=[nothing]`: Styles of the lines. `nothing` will produce a solid line.
"""
function ppHorizontalFlags!(
    figure::Makie.Figure,
    positions::Vector{<:Real};
    colors::Vector{<:ColorType}=[:red],
    line_styles::Vector{<:LineStyleType}=[nothing],
)::Nothing

    # Filter out values outside the range of the original plot
    rangeCut!(positions, ylimits!(figure); keep_edges=false)

    if isempty(positions)

        # Don't draw anything if all the flags are outside the original plot window
        !logging[] || @warn("ppHorizontalFlags!: All horizontal lines lie outside the plot range")

    else

        # Draw the horizontal lines
        for (i, position) in pairs(positions)
            color = ring(colors, i)
            linestyle = ring(line_styles, i)
            hlines!(figure.current_axis.x, position; color, linestyle)
        end

    end

    return nothing

end

"""
    ppCross!(
        figure::Makie.Figure,
        cross_point::Tuple{<:Real,<:Real};
        <keyword arguments>
    )::Nothing

Draw two lines, one horizontal and one vertical.

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `cross_point::Tuple{<:Real,<:Real}`: Crossing point of the lines.
  - `color::ColorType=Makie.wong_colors()[6]`: Color of the lines.
  - `linestyle::LineStyleType=nothing`: Style of the lines. `nothing` will produce a solid line.
"""
function ppCross!(
    figure::Makie.Figure,
    cross_point::Tuple{<:Real,<:Real};
    color::ColorType=Makie.wong_colors()[6],
    linestyle::LineStyleType=nothing,
)::Nothing

    x_limits = xlimits!(figure)
    y_limits = ylimits!(figure)

    if x_limits[1] < cross_point[1] < x_limits[2]
        # Draw the vertical line
        vlines!(figure.current_axis.x, cross_point[1]; color, linestyle)
    else
        !logging[] || @warn("ppCross!: The vertical line lies outside the plot range")
    end

    if y_limits[1] < cross_point[2] < y_limits[2]
        # Draw the horizontal line
        hlines!(figure.current_axis.x, cross_point[2]; color, linestyle)
    else
        !logging[] || @warn("ppCross!: The horizontal line lies outside the plot range")
    end

    return nothing

end

"""
    ppArrows!(
        figure::Makie.Figure,
        positions::Vector{NTuple{4,Float64}};
        <keyword arguments>
    )::Nothing

Draw arrows.

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `positions::Vector{NTuple{4,Float64}}`: The x, y, u, and v coordinates of the arrows. x and y indicate the position of the tails, and u and v the position of the point of the arrows.
  - `colors::Vector{<:ColorType}=[:red]`: Colors of the arrows.
"""
function ppArrows!(
    figure::Makie.Figure,
    positions::Vector{NTuple{4,Float64}};
    colors::Vector{<:ColorType}=[:red],
)::Nothing

    for (i, position) in pairs(positions)

        x, y, u, v = position
        color = ring(colors, i)

        # Draw the arrows
        arrows!(figure.current_axis.x, [x], [y], [u], [v]; color)

    end

    return nothing

end

"""
    ppAnnotation!(figure::Makie.Figure, text::AbstractString; <keyword arguments>)::Nothing

Add an annotation to the plot.

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `text::AbstractString`: Text to be written.
  - `position::Tuple{<:Real,<:Real}=(0.04, 0.98)`: Relative position of the top left corner of the text box.
  - `color=:black`: Text color.
  - `fontsize::Int=35`: Font size.
"""
function ppAnnotation!(
    figure::Makie.Figure,
    text::AbstractString;
    position::Tuple{<:Real,<:Real}=(0.04, 0.98),
    color=:black,
    fontsize::Int=35,
)::Nothing

    text!(
        figure.current_axis.x,
        position[1],
        position[2];
        text,
        align=(:left, :top),
        color,
        space=:relative,
        fontsize,
    )

    return nothing

end

@doc raw"""
    ppFitLine!(
        figure::Makie.Figure;
        <keyword arguments>
    )::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

Draw a linear fit for the data in `figure`.

An annotation with the equation $y = a \, x + b$, and the fitted values for $a$ and $b$, will be positioned in the upper right corner of the plot.

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `top_position::Tuple{<:Real,<:Real}=(0.04, 0.98)`: Relative position of the top left corner of the legend box.
  - `wts::Union{Vector{Float64},Nothing}=nothing`: Weights for the fits. Set to `nothing` for a non-weighted fit.
  - `error_formating::Symbol=:std_error`: Error format for the annotation. The options are:

      + `:std_error`     -> mean ± standard_error.
      + `:conf_interval` -> mean ± max(upper$_{95\%}$ - mean, mean - lower$_{95\%}$).
  - `color::ColorType=Makie.wong_colors()[6],`: Color of the line.
  - `linestyle::LineStyleType=nothing`: Style of the line. `nothing` will produce a solid line.
  - `linewidth::Int=3`: Line width.

# Returns

  - A tuple with the elements for the legend:

      + A `LineElement` to be used as the marker.
      + The label string.
"""
function ppFitLine!(
    figure::Makie.Figure;
    top_position::Tuple{<:Real,<:Real}=(0.04, 0.98),
    wts::Union{Vector{Float64},Nothing}=nothing,
    error_formating::Symbol=:std_error,
    color::ColorType=Makie.wong_colors()[6],
    linestyle::LineStyleType=nothing,
    linewidth::Int=3,
)::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

    # Read the data points in the plot
    points = pointData(figure)

    if isempty(points)
        !logging[] || @warn("ppFitLine!: There are no points in the figure")
        return nothing
    end

    sort!(points; by=x -> x[1])

    # Get the scaling of each axis
    x_scaling = xscale(figure)
    y_scaling = yscale(figure)

    ################################################################################################
    # Linear fit
    ################################################################################################

    # Get the x coordinates of the points
    x_points = x_scaling.(Float64[point[1] for point in points])

    # Get the y coordinates of the points
    y_points = y_scaling.(Float64[point[2] for point in points])

    # Find NaN points
    xnan_idxs = map(isnan, x_points)
    ynan_idxs = map(isnan, y_points)

    # Delete NaN points
    deleteat!(x_points, xnan_idxs ∪ ynan_idxs)
    deleteat!(y_points, xnan_idxs ∪ ynan_idxs)

    # Compute the linear fit in scaled (e.g. log10) space
    X = [ones(length(x_points)) x_points]
    if isnothing(wts)
        linear_model = lm(X, y_points)
    else
        linear_model = lm(X, y_points; wts)
    end

    # Read the fitted coeficients
    coeff = coef(linear_model)
    intercept_mean = coeff[1]
    slope_mean = coeff[2]

    y_limits = [extrema(y_points)...]
    x_limits = @. (y_limits - intercept_mean) / slope_mean

    # Plot the linear fit
    lp = lines!(
        figure.current_axis.x,
        Makie.inverse_transform(x_scaling).(x_limits),
        Makie.inverse_transform(y_scaling).(y_limits);
        color,
        linestyle,
        linewidth,
    )

    ################################################################################################
    # Annotation
    ################################################################################################

    # Compute the fitting errors
    if error_formating == :conf_interval

        conf_int = confint(linear_model)
        intercept_error = max(conf_int[1, 2] - intercept_mean, intercept_mean - conf_int[1, 1])
        slope_error = max(conf_int[2, 2] - slope_mean, slope_mean - conf_int[2, 1])

    elseif error_formating == :std_error

        intercept_error, slope_error = stderror(linear_model)

    else

        throw(ArgumentError("ppFitLine!: `error_formating` can only be :conf_interval or \
        :std_error, but I got :$(error_formating)"))

    end

    # Format the mean and error values
    intercept, δintercept = formatError(intercept_mean, intercept_error)
    slope, δslope = formatError(slope_mean, slope_error)

    # Set the position of the top left corner of the legend box
    x_tp = top_position[1]
    y_tp = top_position[2]

    # Draw the annotation
    text!(
        figure.current_axis.x,
        [(x_tp, y_tp), (x_tp, y_tp - 0.05), (x_tp, y_tp - 0.1)];
        text=[L"y = a \, x + b", L"a = %$slope \pm %$δslope", L"b = %$intercept \pm %$δintercept"],
        align=(:left, :top),
        color,
        space=:relative,
    )

    ################################################################################################
    # Band
    ################################################################################################

    # Compute the intercept and slope as numbers with uncertainties
    band_intercept = intercept_mean ± intercept_error
    band_slope     = slope_mean ± slope_error

    # Compute the values of the y axis as numbers with uncertainties
    x_band = collect(range(x_limits[1], x_limits[2], 100))
    y_band = @. Makie.inverse_transform(y_scaling)(x_band * band_slope + band_intercept)

    values        = Measurements.value.(y_band)
    uncertainties = Measurements.uncertainty.(y_band)

    # Draw the uncertainty band
    bp = band!(
        figure.current_axis.x,
        x_band,
        values .- uncertainties,
        values .+ uncertainties;
        color=(color, 0.3),
    )

    # Put the post processing elements at the back of the plot
    translate!(Accum, bp, 0, 0, -10)
    translate!(Accum, lp, 0, 0, -9)

    return nothing

end

"""
    ppKennicutt1998!(
        figure::Makie.Figure;
        <keyword arguments>
    )::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

Draw a line plot with the fit for the KS relation in Kennicutt (1998).

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `x_unit::Unitful.Units=u"Msun * pc^-2"`: Unit for the area density of gas used in `figure`.
  - `y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2"`: Unit for the area density of star formation rate used in `figure`.
  - `x_log::Bool=true`: If the x axis is ``\\log_{10}(\\Sigma_\\mathrm{HI + H_2})`` (`x_log` = true) or just ``\\Sigma_\\mathrm{HI + H_2}`` (`x_log` = false).
  - `y_log::Bool=true`: If the y axis is ``\\log_{10}(\\Sigma_\\mathrm{SFR})`` (`y_log` = true) or just ``\\Sigma_\\mathrm{SFR}`` (`y_log` = false).
  - `extend::Float64=0.0`: By default the y axis limits of the line will be the vertical range of point in the plot. This can be extended by the fraction `extend` of the vertical range.
  - `colors::Vector{<:ColorType}=[Makie.wong_colors()[6], Makie.wong_colors()[7]]`: Colors for the line. The first color will indicate the range for which there are experimental data, and the second color will be for the extrapolation.
  - `linestyle::LineStyleType=nothing`: Style of the line. `nothing` will produce a solid line.
  - `linewidth::Int=3`: Line width.

# Returns

  - A tuple with the elements for the legend:

      + A `LineElement` to be used as the marker.
      + The label string.

# References

R. C. Kennicutt (1998). *The Global Schmidt Law in Star-forming Galaxies*. The Astrophysical Journal, **498(2)**, 541-552. [doi:10.1086/305588](https://doi.org/10.1086/305588)
"""
function ppKennicutt1998!(
    figure::Makie.Figure;
    x_unit::Unitful.Units=u"Msun * pc^-2",
    y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2",
    x_log::Bool=true,
    y_log::Bool=true,
    extend::Float64=0.0,
    colors::Vector{<:ColorType}=[Makie.wong_colors()[6], Makie.wong_colors()[7]],
    linestyle::LineStyleType=nothing,
    linewidth::Int=3,
)::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

    # Read the data points in the plot
    points = pointData(figure)

    if isempty(points)
        !logging[] || @warn("ppKennicutt1998!: There are no points in the figure")
        return nothing
    end

    # Get the extrema of the y coordinates
    y_limits = [extrema(Float64[point[2] for point in points])...]

    # Compute the vertical range, and if it is too small extend it
    y_range = y_limits[2] - y_limits[1]

    if iszero(extend) && (abs(y_range / y_limits[1]) < 0.1)
        extend = 0.4 / abs(y_range / y_limits[1])
    end

    # Extend the y limits the required amount
    if extend > 0.0

        extension = y_range * extend

        y_limits[1] -= extension
        y_limits[2] += extension

    end

    # Set the correct unit and scale for the range of experimental values
    kennicutt_range = ustrip.(y_unit, KS98_SFR_RANGE)

    if y_log
        kennicutt_range = log10.(kennicutt_range)
    end

    # Set the y ranges to be plotted
    y_ranges = [kennicutt_range]

    # If there is extrapolation add new ranges
    if y_limits[2] > kennicutt_range[2]
        push!(y_ranges, [kennicutt_range[2], y_limits[2]])
    end
    if y_limits[1] < kennicutt_range[1]
        push!(y_ranges, [y_limits[1], kennicutt_range[1]])
    end

    for (y_zone, color) in zip(y_ranges, [colors..., colors[2]])

        # Compute the extrema of the star formation area density
        if y_log
            Σsfr = @. exp10(y_zone) * y_unit
        else
            Σsfr = @. y_zone * y_unit
        end

        # Compute the extrema of the gas mass area density
        Σgas     = invKennicutt1998(Σsfr; log_output=false)
        x_limits = Measurements.value.(Σgas)

        # Compute the values for the x axis
        x_points = collect(range(x_limits[1], x_limits[2], 100))

        if x_log
            x_axis = @. log10(ustrip(x_unit, x_points))
        else
            x_axis = @. ustrip(x_unit, x_points)
        end

        # Compute the values for the y axis
        y_points = kennicutt1998(x_points; log_output=false)

        if y_log
            y_axis = @. log10(ustrip(y_unit, y_points))
        else
            y_axis = @. ustrip(y_unit, y_points)
        end

        values        = Measurements.value.(y_axis)
        uncertainties = Measurements.uncertainty.(y_axis)

        bp = band!(
            figure.current_axis.x,
            x_axis,
            values .- uncertainties,
            values .+ uncertainties;
            color=(color, 0.3),
        )

        lp = lines!(figure.current_axis.x, x_axis, values; color, linestyle, linewidth)

        # Put the post processing elements at the back of the plot
        translate!(Accum, bp, 0, 0, -10)
        translate!(Accum, lp, 0, 0, -9)

    end

    return [LineElement(; color=colors[1], linestyle, linewidth)], ["Kennicutt (1998)"]

end

"""
    ppBigiel2008!(
        figure::Makie.Figure,
        molecular::Bool;
        <keyword arguments>
    )::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

Draw a line plot with the fit for the KS law, taken from Bigiel et al. (2008).

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `molecular::Bool`: If the x axis will be the area mass density of molecular hydrogen, or, if set to `false`, the area mass density of neutral hydrogen.
  - `x_unit::Unitful.Units=u"Msun * pc^-2"`: Unit for the area density of gas used in `figure`.
  - `y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2"`: Unit for the area density of star formation rate used in `figure`.
  - `x_log::Bool=true`: If the x axis is ``\\log_{10}(\\Sigma_\\mathrm{H})`` (`x_log` = true) or just ``\\Sigma_\\mathrm{H}`` (`x_log` = false).
  - `y_log::Bool=true`: If the y axis is ``\\log_{10}(\\Sigma_\\mathrm{SFR})`` (`y_log` = true) or just ``\\Sigma_\\mathrm{SFR}`` (`y_log` = false).
  - `extend::Float64=0.0`: By default the y axis limits will be the vertical range of the point in the plot. This can be extended by the multiplicative factor `extend` of the vertical range.
  - `colors::Vector{<:ColorType}=[Makie.wong_colors()[6], Makie.wong_colors()[7]]`: Colors for the line. The first color will indicate the range for which there are experimental data, and the second color will be for the extrapolation.
  - `linestyle::LineStyleType=nothing`: Style of the line. `nothing` will produce a solid line.
  - `linewidth::Int=3`: Line width.

# Returns

  - A tuple with the elements for the legend:

      + A `LineElement` to be used as the marker.
      + The label string.

# References

F. Bigiel et al. (2008). *THE STAR FORMATION LAW IN NEARBY GALAXIES ON SUB-KPC SCALES*. The Astrophysical Journal, **136(6)**, 2846. [doi:10.1088/0004-6256/136/6/2846](https://doi.org/10.1088/0004-6256/136/6/2846)
"""
function ppBigiel2008!(
    figure::Makie.Figure,
    molecular::Bool;
    x_unit::Unitful.Units=u"Msun * pc^-2",
    y_unit::Unitful.Units=u"Msun * yr^-1 * kpc^-2",
    x_log::Bool=true,
    y_log::Bool=true,
    extend::Float64=0.0,
    colors::Vector{<:ColorType}=[Makie.wong_colors()[6], Makie.wong_colors()[7]],
    linestyle::LineStyleType=nothing,
    linewidth::Int=3,
)::Union{Tuple{Vector{<:LegendElement},Vector{AbstractString}},Nothing}

    # Read the data points in the plot
    points = pointData(figure)

    if isempty(points)
        !logging[] || @warn("ppBigiel2008!: There are no points in the figure")
        return nothing
    end

    # Get the extrema of the y coordinates
    y_limits = [extrema(Float64[point[2] for point in points])...]

    # Compute the vertical range, and if it is too small extend it
    y_range = y_limits[2] - y_limits[1]

    if iszero(extend) && (abs(y_range / y_limits[1]) < 0.1)
        extend = 0.4 / abs(y_range / y_limits[1])
    end

    # Extend the y limits the required amount
    if extend > 0.0

        extension = y_range * extend

        y_limits[1] -= extension
        y_limits[2] += extension

    end

    # Set the correct unit and scale for the range of experimental values
    bigiel_range = ustrip.(y_unit, BIGIEL2008_SFR_RANGE)

    if y_log
        bigiel_range = log10.(bigiel_range)
    end

    # Set the y ranges to be plotted
    y_ranges = [bigiel_range]

    # If there is extrapolation add new ranges
    if y_limits[2] > bigiel_range[2]
        push!(y_ranges, [bigiel_range[2], y_limits[2]])
    end
    if y_limits[1] < bigiel_range[1]
        push!(y_ranges, [y_limits[1], bigiel_range[1]])
    end

    for (y_zone, color) in zip(y_ranges, [colors..., colors[2]])

        # Compute the extrema of the star formation area density
        if y_log
            Σsfr = @. exp10(y_zone) * y_unit
        else
            Σsfr = @. y_zone * y_unit
        end

        # Compute the extrema of the gas mass area density
        ΣH       = invBigiel2008(Σsfr; molecular, log_output=false)
        x_limits = Measurements.value.(ΣH)

        # Compute the values for the x axis
        x_points = collect(range(x_limits[1], x_limits[2], 100))

        if x_log
            x_axis = @. log10(ustrip(x_unit, x_points))
        else
            x_axis = @. ustrip(x_unit, x_points)
        end

        # Compute the values for the y axis
        y_points = bigiel2008(x_points; molecular, log_output=false)

        if y_log
            y_axis = @. log10(ustrip(y_unit, y_points))
        else
            y_axis = @. ustrip(y_unit, y_points)
        end

        values        = Measurements.value.(y_axis)
        uncertainties = Measurements.uncertainty.(y_axis)

        bp = band!(
            figure.current_axis.x,
            x_axis,
            values .- uncertainties,
            values .+ uncertainties;
            color=(color, 0.25),
        )

        lp = lines!(figure.current_axis.x, x_axis, values; color, linestyle, linewidth)

        # Put the post processing elements at the back of the plot
        translate!(Accum, bp, 0, 0, -10)
        translate!(Accum, lp, 0, 0, -9)

    end

    return [LineElement(; color=colors[1], linestyle, linewidth)], ["Bigiel et al. (2008)"]

end

"""
    ppBigiel2010!(
	    figure::Makie.Figure;
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Draw a scatter plot of the SFR surface density vs gas surface density (Kennicutt-Schmidt law) for a given galaxy, using the data of Bigiel et al. (2010).

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `galaxy::Union{String,Symbol}="NGC 628"`: Target galaxy. The options are:

      + With molecular and and atomic data (Table 2): "NGC 628", "NGC 3184", "NGC 3521", "NGC 4736", "NGC 5055", "NGC 5194", "NGC 6946".
      + With only atomic data (Table 3): "NGC 925", "NGC 2403", "NGC 2841", "NGC 2903", "NGC 3198", "NGC 3351", "NGC 3621", "NGC 3627", "NGC 5236", "NGC 5457", "NGC 7331", "NGC 7793".
      + :all: Every galaxy that is available for the given `quantity`.
    For more information on each galaxy see Bigiel et al. (2010).
  - `quantity::Symbol=:molecular`: Gas quantity for the x axis. The options are:

      + `:molecular` -> Surface density of molecular gas.
      + `:neutral`   -> Surface density of neutral gas.
      + `:atomic`    -> Surface density of atomic gas.
  - `x_log::Bool=true`: If the x axis will be plotted as the log10 of the gas surface density.
  - `y_log::Bool=true`: If the y axis will be plotted as the log10 of the SFR surface density.
  - `x_unit::Unitful.Units=u"Msun * kpc^-2"`: Unit for the x axis.
  - `y_unit::Unitful.Units=u"Msun * yr^-1 *  kpc^-2"`: Unit for the y axis.
  - `color::ColorType=Makie.wong_colors()[2]`: Color of the markers.

# Returns

  - A tuple with the elements for the legend:

      + A `MarkerElement` to be used as the marker.
      + The label string.

# References

F. Bigiel et al. (2010). *EXTREMELY INEFFICIENT STAR FORMATION IN THE OUTER DISKS OF NEARBY GALAXIES*. The Astrophysical Journal, **140(5)**, 1194. [doi:10.1088/0004-6256/140/5/1194](https://doi.org/10.1088/0004-6256/140/5/1194)
"""
function ppBigiel2010!(
	figure::Makie.Figure;
	galaxy::Union{String,Symbol}="NGC 628",
	quantity::Symbol=:molecular,
	x_log::Bool=true,
	y_log::Bool=true,
	x_unit::Unitful.Units=u"Msun * kpc^-2",
	y_unit::Unitful.Units=u"Msun * yr^-1 *  kpc^-2",
    color::ColorType=Makie.wong_colors()[2],
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

	################################################################################################
	# Load table 2 from Bigiel et al. 2010
	################################################################################################

	raw_data_2 = readdlm(BIGIEL2010_TABLE_2, '\t', skipstart=48, header=true)[1]

	table_2 = DataFrame(
		gtype=String[],
		name=String[],
		logHI=Union{Float64,Missing}[],
		logH2=Union{Float64,Missing}[],
		logSFR=Union{Float64,Missing}[],
	)

	for row in eachrow(raw_data_2)

		data = row[1]

		gtype  = strip(data[1:13])
		name   = strip(data[14:24])
		logHI  = parserWS(data[26:29])
		logH2  = parserWS(data[36:39])
		logSFR = parserWS(data[46:50])

		push!(table_2, [gtype name logHI logH2 logSFR])

	end

	spirals_2 = filter(:gtype => isequal("Spirals"), table_2)

	################################################################################################
	# Load table 3 from Bigiel et al. 2010
	################################################################################################

	raw_data_3 = readdlm(BIGIEL2010_TABLE_3, '\t', skipstart=46, header=true)[1]

	table_3 = DataFrame(
		gtype=String[],
		name=String[],
		logHI=Union{Float64,Missing}[],
		SFR=Union{Float64,Missing}[],
	)

	for row in eachrow(raw_data_3)

		data = row[1]

		gtype = strip(data[1:8])
		name  = strip(data[9:19])
		logHI = parserWS(data[21:25])
		SFR   = parserWS(data[32:37]) .* exp10(-5.0)

		push!(table_3, [gtype name logHI SFR])

	end

	spirals_3 = filter(:gtype => isequal("Spirals"), table_3)

	################################################################################################
	# Find the target galaxy and read its values
	################################################################################################

	if quantity == :molecular

        if galaxy == :all
            target_galaxy = spirals_2
        else
		    target_galaxy = filter(:name => isequal(galaxy), spirals_2)
        end

		if isempty(target_galaxy)

			throw(ArgumentError("ppBigiel2010!: `galaxy` = $(galaxy) is not a spiral galaxy in \
            Table 2 of Bigiel et al. 2010 (notice that you have `quantity` = $(quantity))."))

		end

		Σg   = exp10.(target_galaxy[!, :logH2]) .* u"Msun * pc^-2"
		Σsfr = exp10.(target_galaxy[!, :logSFR]) .* u"Msun * yr^-1 * kpc^-2"

    elseif quantity == :neutral

        if galaxy == :all
            target_galaxy = spirals_2
        else
		    target_galaxy = filter(:name => isequal(galaxy), spirals_2)
        end

		if isempty(target_galaxy)

			throw(ArgumentError("ppBigiel2010!: `galaxy` = $(galaxy) is not a spiral galaxy in \
            Table 2 of Bigiel et al. 2010 (notice that you have `quantity` = $(quantity))."))

		end

        ΣH2   = exp10.(target_galaxy[!, :logH2]) .* u"Msun * pc^-2"
        ΣHI   = exp10.(target_galaxy[!, :logHI]) .* u"Msun * pc^-2"

		Σg   = ΣH2 .+ ΣHI
		Σsfr = exp10.(target_galaxy[!, :logSFR]) .* u"Msun * yr^-1 * kpc^-2"

	elseif quantity == :atomic

		if galaxy == :all

            # Find the galaxies that are exclusive to table 3
            galaxies_t3 = setdiff(unique(spirals_3[!, :name]), unique(spirals_2[!, :name]))

            # Join the data in table 2 and table 3, ignoring repeated galaxies
            target_galaxy = vcat(
                spirals_2[!, [:logHI, :SFR]],
                filter(:name => in(galaxies_t3), spirals_3)[!, [:logHI, :SFR]],
            )

        else

		    target_galaxy = filter(:name => isequal(galaxy), spirals_2)

        end

		if isempty(target_galaxy)

			target_galaxy = filter(:name => isequal(galaxy), spirals_3)

			if isempty(target_galaxy)

				throw(ArgumentError("ppBigiel2010!: `galaxy` $(galaxy) is not a spiral galaxy in \
                Table 2 or 3 of Bigiel et al. 2010."))

			end

		end

		Σg   = exp10.(target_galaxy[!, :logHI]) .* u"Msun * pc^-2"
		Σsfr = target_galaxy[!, :SFR] .* u"Msun * yr^-1 * kpc^-2"

	else

		throw(ArgumentError("ppBigiel2010!: `quantity` can only be :molecular, :neutral or \
        :atomic, but I got :$(quantity)."))

	end

	# Delete missing data
	idxs = map(ismissing, Σg) ∪ map(ismissing, Σsfr)

    deleteat!(Σg, idxs)
	deleteat!(Σsfr, idxs)

	# Set the correct scale and units
	if x_log
		Σg = log10.(ustrip.(x_unit, Σg))
	else
		Σg = ustrip.(x_unit, Σg)
	end

	if y_log
		Σsfr = log10.(ustrip.(y_unit, Σsfr))
	else
		Σsfr = ustrip.(y_unit, Σsfr)
	end

	################################################################################################
	# Plot the galactic data
	################################################################################################

	sp = scatter!(
        figure.current_axis.x,
        Σg,
        Σsfr;
        color=(color, 0.5),
        marker=:star4,
        markersize=10,
    )

    # Put the post processing elements at the back of the plot
    translate!(Accum, sp, 0, 0, -10)

    if galaxy == :all
        label = "Bigiel et al. 2010"
    else
        label = "$(galaxy) - Bigiel et al. 2010"
    end

    return ([MarkerElement(; color, marker=:star4)], [label])

end

@doc raw"""
    ppSun2023!(
	    figure::Makie.Figure;
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Draw a scatter plot of the SFR surface density vs molecular surface density (molecular Kennicutt-Schmidt law) for a given galaxy, using the data of Sun et al. (2023).

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `galaxy::Union{String,Symbol}=:main`: Target galaxy. The options are:

      + One of the 80 galaxies in the dataset, e.g. "ESO097-013", "IC1954", "IC5273", "NGC1546", "NGC1559", "NGC1566", etc. For a full list see the reference below.
      + :main: Every galaxy with $-2.0 < \log_{10}(t_\mathrm{dep} \, / \, \mathrm{Gyr}) < 2.0$, where $t_\mathrm{dep} = \Sigma_\mathrm{H_2} / \Sigma_\mathrm{SFR}$ is the depletion time.
      + :all: All 80 galaxies in the dataset.
    For more information on each galaxy see Sun et al. (2023).
  - `sfr_calibration::Symbol=:Halpha`: SFR calibration for combining UV, optical, and/or IR data. The options are: `:Halpha`, `:FUV`, and `:AV_corrected_Halpha`. For an explanation of each one see section 2 of Sun et al. (2023).
  - `h2_prescription::Symbol=:S20`: Prescription for the CO-to-H₂ conversion factor. The options are: `:S20`, `:Mw`, `:B13`, and `:G20`. For an explanation of each one see section 2 of Sun et al. (2023).
  - `x_log::Bool=true`: If the x axis will be plotted as the log10 of the gas surface density.
  - `y_log::Bool=true`: If the y axis will be plotted as the log10 of the SFR surface density.
  - `x_unit::Unitful.Units=u"Msun * pc^-2"`: Unit for the x axis.
  - `y_unit::Unitful.Units=u"Msun * yr^-1 *  kpc^-2"`: Unit for the y axis.
  - `color::ColorType=Makie.wong_colors()[2]`: Color of the markers.

# Returns

  - A tuple with the elements for the legend:

      + A `MarkerElement` to be used as the marker.
      + The label string.

# References

J. Sun et al. (2023). *Star Formation Laws and Efficiencies across 80 Nearby Galaxies*. The Astrophysical Journal Letters, **945(2)**, L19. [doi:10.3847/2041-8213/acbd9c](https://doi.org/10.3847/2041-8213/acbd9c)
"""
function ppSun2023!(
	figure::Makie.Figure;
	galaxy::Union{String,Symbol}=:main,
    sfr_calibration::Symbol=:Halpha,
	h2_prescription::Symbol=:S20,
	x_log::Bool=true,
	y_log::Bool=true,
	x_unit::Unitful.Units=u"Msun * pc^-2",
	y_unit::Unitful.Units=u"Msun * yr^-1 *  kpc^-2",
    color::ColorType=Makie.wong_colors()[2],
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

	################################################################################################
	# Load table A1 from Sun et al. 2023
	################################################################################################

	raw_data = readdlm(SUN2023_TABLE, skipstart=57, header=false)

    clean_data = DataFrame(replace(raw_data, "" => Inf), :auto)

    # Shift in the column indices correspondig to each SFR calibration
    sfr_calibrations = Dict(:Halpha => 0, :FUV => 2, :AV_corrected_Halpha => 4)

    # Shift in the column indices correspondig to each CO-to-H₂ prescription
	h2_prescriptions = Dict(:S20 => 0, :Mw => 2, :B13 => 4, :G20 => 6)

    # List of available galaxies
    galaxies = unique(raw_data[:, 1])

	################################################################################################
	# Find the target galaxy and read its values
	################################################################################################

    if isa(galaxy, String)

        (
            galaxy ∈ galaxies ||
            throw(ArgumentError("ppSun2023!: `galaxy` = $(galaxy) is not a valid galaxy"))
        )

        data = filter(:x1 => isequal(galaxy), clean_data)

    else

        (
            galaxy ∈ [:main, :all] ||
            throw(
                ArgumentError("ppSun2023!: `galaxy` can noly be :main or :all but I got :$(galaxy)")
            )
        )

        data = clean_data

    end

    Σsfr = data[:, 4 + sfr_calibrations[sfr_calibration]] .* u"Msun * yr^-1 * kpc^-2"
    Σh2  = data[:, 10 + h2_prescriptions[h2_prescription]] .* u"Msun * pc^-2"

    # Set the correct scale and units
	if x_log
        x_data = @. log10(ustrip(x_unit, Σh2))
    else
        x_data = @. ustrip(x_unit, Σh2)
    end

    if y_log
        y_data = @. log10(ustrip(y_unit, Σsfr))
    else
        y_data = @. ustrip(y_unit, Σsfr)
    end

    # Delete missing data
    filter = x -> isnan(x) || isinf(x)
    idxs   = map(filter, x_data) ∪ map(filter, y_data)

    if galaxy == :main
        # Compute the depletion time
        tdep = @. log10(ustrip(u"Gyr", Σh2 ./ Σsfr))

        # Filter galaxies with tdep outside the range [-2.0, 2.0]
        tdep_filter = x -> isnan(x) || isinf(x) || x < -2.0 || x > 2.0

        idxs = idxs ∪ map(tdep_filter, tdep)
    end

    deleteat!(x_data, idxs)
    deleteat!(y_data, idxs)

	################################################################################################
	# Plot the galactic data
	################################################################################################

	sp = scatter!(
        figure.current_axis.x,
        x_data,
        y_data;
        color=(color, 0.5),
        marker=:star4,
        markersize=10,
    )

    # Put the post processing elements at the back of the plot
    translate!(Accum, sp, 0, 0, -10)

    if isa(galaxy, String)
        label = "$(galaxy) - Sun et al. 2023"
    elseif galaxy == :all
        label = "Sun et al. 2023 (all galaxies)"
    else
        label = "Sun et al. 2023"
    end

    return ([MarkerElement(; color, marker=:star4)], [label])

end

"""
    ppMolla2015!(
        figure::Makie.Figure,
        quantity::Symbol,
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Draw a profile for the Milky Way using the data compiled by Mollá et al. (2015).

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_area_density`      -> Stellar area mass density.
      + `:molecular_area_density`    -> Molecular mass surface density.
      + `:br_molecular_area_density` -> Molecular mass surface density, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_area_density`       -> Atomic hydrogen area mass density.
      + `:sfr_area_density`          -> Star formation rate area density.
      + `:X_stellar_abundance`       -> Stellar abundance of element ``\\mathrm{X}``, as ``12 + \\log_{10}(\\mathrm{X \\, / \\, H})``. ``\\mathrm{X}`` can be O (oxygen), N (nitrogen), or C (carbon).
  - `color::ColorType=Makie.wong_colors()[6]`: Color of the line.
  - `linestyle::LineStyleType=nothing`: Style of the line. `nothing` will produce a solid line.
  - `error_bars::Bool=true`: If the error bars will be plotted.

# Returns

  - A tuple with the elements for the legend:

      + A `MarkerElement` to be used as the marker.
      + The label string.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)

M. Mollá et al. (2015). *Galactic chemical evolution: stellar yields and the initial mass function*. Monthly Notices of the Royal Astronomical Society **451(4)**, 3693–3708. [doi:10.1093/mnras/stv1102](https://doi.org/10.1093/mnras/stv1102)
"""
function ppMolla2015!(
    figure::Makie.Figure,
    quantity::Symbol;
    color::ColorType=Makie.wong_colors()[6],
    linestyle::LineStyleType=nothing,
    error_bars::Bool=true,
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

    # Read the file with the compiled data
    raw = CSV.read(
        MOLLA2015_DATA_PATH,
        DataFrame;
        header=[
            "R",
            "ΣHI",
            "ΣHI error",
            "ΣH2",
            "ΣH2 error",
            "logΣ*",
            "logΣ* error",
            "logΣsfr",
            "logΣsfr error",
            "C/H",
            "ΔC/H",
            "N/H",
            "ΔN/H",
            "O/H",
            "ΔO/H",
        ],
    )

    x_values = raw[!, "R"]

    # Select the quantity for the y axis
    if quantity == :stellar_area_density
        y_data = 10 .^ (raw[!, "logΣ*"] .± raw[!, "logΣ* error"])
    elseif quantity ∈ [:molecular_area_density, :br_molecular_area_density]
        y_data = raw[!, "ΣH2"] .± raw[!, "ΣH2 error"]
    elseif quantity == :atomic_area_density
        y_data = raw[!, "ΣHI"] .± raw[!, "ΣHI error"]
    elseif quantity == :sfr_area_density
        y_data = 10 .^ (raw[!, "logΣsfr"] .± raw[!, "logΣsfr error"])
    elseif quantity == :O_stellar_abundance
        y_data = raw[!, "O/H"] .± raw[!, "ΔO/H"]
    elseif quantity == :N_stellar_abundance
        y_data = raw[!, "N/H"] .± raw[!, "ΔN/H"]
    elseif quantity == :C_stellar_abundance
        y_data = raw[!, "C/H"] .± raw[!, "ΔC/H"]
    else
        throw(ArgumentError("ppMolla2015: `x_quantity` can only be  :stellar_area_density, \
        :molecular_area_density, :br_molecular_area_density, :atomic_area_density, \
        :sfr_area_density, :O_stellar_abundance, :N_stellar_abundance, \
        or :C_stellar_abundance, but I got :$(quantity)"))
    end

    y_values = Measurements.value.(y_data)
    y_uncertainties = Measurements.uncertainty.(y_data)

    # Plot the mean values
    scatterlines!(figure.current_axis.x, x_values, y_values; color, linestyle, marker=:utriangle)

    if error_bars
        # Plot the error bars
        errorbars!(figure.current_axis.x, x_values, y_values, y_uncertainties; color)
    end

    return [MarkerElement(; color, marker=:utriangle)], ["Mollá et al. (2015)"]

end

"""
    ppFeldmann2020!(
        figure::Makie.Figure,
        x_quantity::Symbol,
        y_quantity::Symbol;
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Draw a line, or scatter, plot using the experimental data from the xGASS and xCOLD GASS collaborations, processed by Feldmann (2020).

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `x_quantity::Symbol`: Quantity for the x axis. The options are:

      + `:stellar_mass`      -> Stellar mass.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:observational_sfr` -> The star formation rate of the last `AGE_RESOLUTION`.
  - `y_quantity::Symbol`: Quantity for the y axis. The options are:

      + `:stellar_mass`      -> Stellar mass.
      + `:molecular_mass`    -> Molecular hydrogen (``\\mathrm{H_2}``) mass.
      + `:br_molecular_mass` -> Molecular hydrogen (``\\mathrm{H_2}``) mass, computed using the pressure relation in Blitz et al. (2006).
      + `:atomic_mass`       -> Atomic hydrogen (``\\mathrm{HI}``) mass.
      + `:observational_sfr` -> The star formation rate of the last `AGE_RESOLUTION`.
  - `scatter::Bool=false`: If the data will be presented as a line plot with error bands (default), or alternatively, a scatter plot.

# Returns

  - A tuple with the elements for the legend:

      + A `MarkerElement` to be used as the marker.
      + The label string.

# References

L. Blitz et al. (2006). *The Role of Pressure in GMC Formation II: The H2-Pressure Relation*. The Astrophysical Journal, **650(2)**, 933. [doi:10.1086/505417](https://doi.org/10.1086/505417)

R. Feldmann (2020). *The link between star formation and gas in nearby galaxies*. Communications Physics **3(226)**. [doi:10.1038/s42005-020-00493-0](https://doi.org/10.1038/s42005-020-00493-0)
"""
function ppFeldmann2020!(
    figure::Makie.Figure,
    x_quantity::Symbol,
    y_quantity::Symbol;
    scatter::Bool=false,
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

    # Read the CSV file with the raw data
    raw = CSV.read(FELDMANN2020_DATA_PATH, DataFrame, comment="#")

    # Select the quantity for the x axis
    if x_quantity == :stellar_mass

        x_data = exp10.(raw[!, "lgMstar"])

    elseif x_quantity ∈ [:molecular_mass, :br_molecular_mass]

        x_data = raw[!, "MH2"]

    elseif x_quantity == :atomic_mass

        x_data = raw[!, "MHI"]

    elseif x_quantity == :observational_sfr

        x_data = raw[!, "SFR"]

    else

        throw(ArgumentError("ppFeldmann2020!: `x_quantity` can only be :stellar_mass, \
        :molecular_mass, :br_molecular_mass, :atomic_mass, or :observational_sfr, \
        but I got :$(x_quantity)"))

    end

    # Select the quantity for the y axix; with its uncertainty
    if y_quantity == :stellar_mass

        # Compute the mean "error" for the stellar mass
        err_low = raw[!, "lgMstar"] .- raw[!, "lgMstar_p16"]
        err_high = raw[!, "lgMstar_p84"] .- raw[!, "lgMstar"]
        err_mean = @. (err_low + err_high) / 2.0

        y_data = exp10.(raw[!, "lgMstar"] .+ err_mean)

    elseif y_quantity ∈ [:molecular_mass, :br_molecular_mass]

        y_data = raw[!, "MH2"] .± raw[!, "e_MH2"]

    elseif y_quantity == :atomic_mass

        y_data = raw[!, "MHI"] .± raw[!, "e_MHI"]

    elseif y_quantity == :observational_sfr

        y_data = raw[!, "SFR"] .± raw[!, "e_SFR"]

    else

        throw(ArgumentError("ppFeldmann2020!: `y_quantity` can only be :stellar_mass, \
        :molecular_mass, :br_molecular_mass, :atomic_mass, or :observational_sfr, \
        but I got :$(y_quantity)"))

    end

    if scatter
        y_mean_values = Measurements.value.(y_data)

        # Plot the selected values as a scatter plot
        scatter!(
            figure.current_axis.x,
            x_data,
            y_mean_values;
            color=Makie.wong_colors()[6],
            marker=:star4,
            markersize=12,
        )

        return [MarkerElement(; color=Makie.wong_colors()[6], marker=:star4)], ["Feldmann (2020)"]
    end

    # Get the scaling of each axis
    x_scaling = xscale(figure)
    y_scaling = yscale(figure)

    # Scale the raw data
    scaled_xs = x_scaling.(x_data)
    scaled_ys = y_scaling.(y_data)

    # Set up the grid
    x_min, x_max = extrema(skipmissing(scaled_xs))
    grid         = CircularGrid(x_max, 30; shift=x_min)
    n_bins       = length(grid.grid)
    bin_width    = (x_max - x_min) / n_bins

    # Allocate memory for the histogram of indices
    histogram = Vector{Int}[Vector{Int}[] for _ in 1:n_bins]

    # Compute the histogram; ignoring missings and values outside the grid range
    for (x_idx, x_value) in pairs(scaled_xs)

        if ismissing(x_value)
            continue
        elseif x_value < x_min || x_max < x_value
            continue
        elseif x_value == x_min
            hist_idx = 1
        elseif x_value == x_max
            hist_idx = n_bins
        else
            hist_idx = ceil(Int, (x_value - x_min) / bin_width)
        end

        push!(histogram[hist_idx], x_idx)

    end

    ################################################################################################
    # Values for the y axis
    ################################################################################################

    # Allocate memory
    y_axis = Vector{Union{Measurement{Float64},Missing}}(undef, n_bins)

    for (i, hist_idxs) in pairs(histogram)

        # For each bin select the corresponding y measurements
        binned_y_measurements = skipmissing(scaled_ys[hist_idxs])

        if isempty(binned_y_measurements)

            y_axis[i] = missing

        else

            # Compute the weighted mean of all the y measurements in this bin
            weighted_mean = weightedmean(binned_y_measurements)

            # Compute the y mean
            y_mean = Measurements.value(weighted_mean)

            # Compute y uncertainty
            spread = std(Measurements.value.(binned_y_measurements); corrected=false)
            stderr = Measurements.uncertainty(weighted_mean)
            y_err  = sqrt(spread^2 + stderr^2)

            y_axis[i] = y_mean ± y_err

        end

    end

    # Find the empty bins
    missing_idxs = map(ismissing, y_axis)

    # Delete the empty bins
    deleteat!(y_axis, missing_idxs)

    ################################################################################################
    # Values for the x axis
    ################################################################################################

    x_axis = Makie.inverse_transform(x_scaling).(grid.grid)

    # Delete the empty bins
    deleteat!(x_axis, missing_idxs)

    ################################################################################################
    # Plots
    ################################################################################################

    y_axis_value       = Measurements.value.(y_axis)
    y_axis_uncertainty = Measurements.uncertainty.(y_axis)

    # Plot the mean line
    lp = lines!(
        figure.current_axis.x,
        x_axis,
        Makie.inverse_transform(y_scaling).(y_axis_value);
        color=(:red, 0.4),
    )

    # Construct the 1σ band
    yupper_1σ = Makie.inverse_transform(y_scaling).(y_axis_value .+ y_axis_uncertainty)
    ylower_1σ = Makie.inverse_transform(y_scaling).(y_axis_value .- y_axis_uncertainty)
    # Plot the 1σ band
    band!(figure.current_axis.x, x_axis, ylower_1σ, yupper_1σ; color=(:red, 0.15))

    # Construct the 2σ upper band
    high_band_2σ = Makie.inverse_transform(y_scaling).(y_axis_value .+ 2 * y_axis_uncertainty)
    # Plot the 2σ upper band
    band!(figure.current_axis.x, x_axis, yupper_1σ, high_band_2σ; color=(:orange, 0.15))

    # Construct the 2σ lower band
    low_band_2σ = Makie.inverse_transform(y_scaling).(y_axis_value .- 2 * y_axis_uncertainty)
    # Plot the 2σ lower band
    bp = band!(figure.current_axis.x, x_axis, low_band_2σ, ylower_1σ; color=(:orange, 0.15))

    # Put the post processing elements at the back of the plot
    translate!(Accum, bp, 0, 0, -10)
    translate!(Accum, lp, 0, 0, -9)

    return (
        [PolyElement(; color=(:red, 0.3)), PolyElement(; color=(:orange, 0.3))],
        [
            L"\mathrm{Feldmann \,\, (2020)} \,\, 1 \, \sigma",
            L"\mathrm{Feldmann \,\, (2020)} \,\, 2 \, \sigma",
        ],
    )

end

"""
    ppBarPlotLabels(
        ::Makie.Figure,
        include_stars::Bool;
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Return the legend elements for the plot made by [`gasBarPlot`](@ref).

# Arguments

  - `::Makie.Figure`: Makie figure to be drawn over.
  - `include_stars::Bool=false`: If the stars will be included as one of the gas phases.
  - `colors=Makie.wong_colors()`: Colors for the bars.

# Returns

  - A tuple with the elements for the legend:

      + A `PolyElement` to be used in the legend.
      + The label strings.
"""
function ppBarPlotLabels(
    ::Makie.Figure,
    include_stars::Bool;
    colors=Makie.wong_colors(),
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

    if include_stars
        labels = ["Molecular gas", "Atomic gas", "Ionized gas", "Stars"]
    else
        labels = ["Molecular gas", "Atomic gas", "Ionized gas"]
    end

    return [PolyElement(polycolor=colors[i]) for i in 1:length(labels)], labels

end

"""
    ppAgertz2021!(
        figure::Makie.Figure;
        <keyword arguments>
    )::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

Draw a line plot using the experimental data and fits from McMillan (2011) and Leroy et al. (2008), show in Agertz et al. (2021). for the stellar density profile.

# Arguments

  - `figure::Makie.Figure`: Makie figure to be drawn over.
  - `galaxies::Vector{String}=["MW"]`: Target galaxies. The options are:

      + One of the 23 galaxies in Leroy et al. (2008), e.g. "DDO154", "HOI", "HOII", "IC2574", "NGC0628", "NGC0925", etc. For a full list see the reference below.
      + "MW": The Milky Way fits from McMillan (2011).
  - `x_unit::Unitful.Units=u"kpc"`: Unit for the x axis.
  - `y_unit::Unitful.Units=u"Msun * kpc^-2"`: Unit for the y axis.
  - `x_log::Bool=false`: If the x axis will be plotted as the log10 of the galactocentric radius.
  - `y_log::Bool=true`: If the y axis will be plotted as the log10 of the stellar surface density.
  - `error_band::Bool=true`: If the error band will be plotted.
  - `colors::Vector{<:ColorType}=[Makie.wong_colors()[2]]`: Colors for the lines.
  - `linestyle::LineStyleType=nothing`: Style of the lines. `nothing` will produce a solid line.
  - `linewidth::Int=3`: Width of the lines.

# Returns

  - A tuple with the elements for the legend:

      + A vector of `LineElement`.
      + A vector of labels.

# References

O. Agertz et al. (2021). *VINTERGATAN – I. The origins of chemically, kinematically, and structurally distinct discs in a simulated Milky Way-mass galaxy*. Monthly Notices of the Royal Astronomical Society. **503(4), 5826–5845. [doi:10.1093/mnras/stab322](https://doi.org/10.1093/mnras/stab322)

A. K. Leroy et al. (2008). *THE STAR FORMATION EFFICIENCY IN NEARBY GALAXIES: MEASURING WHERE GAS FORMS STARS EFFECTIVELY*. The Astronomical Journal **136(6)**, 2782–2845. [doi:10.1088/0004-6256/136/6/2782](https://doi.org/10.1088/0004-6256/136/6/2782)

P. J. McMillan (2011). *Mass models of the Milky Way*. Monthly Notices of the Royal Astronomical Society **414(3)**, 2446–2457. [doi:10.1111/j.1365-2966.2011.18564.x](https://doi.org/10.1111/j.1365-2966.2011.18564.x)
"""
function ppAgertz2021!(
    figure::Makie.Figure;
    galaxies::Vector{String}=["MW"],
    x_unit::Unitful.Units=u"kpc",
    y_unit::Unitful.Units=u"Msun * kpc^-2",
    x_log::Bool=false,
    y_log::Bool=true,
    error_band::Bool=true,
    colors::Vector{<:ColorType}=[Makie.wong_colors()[2]],
    linestyle::LineStyleType=nothing,
    linewidth::Int=3,
)::Tuple{Vector{<:LegendElement},Vector{AbstractString}}

    ################################################################################################
	# Load the data from McMillan (2011) and Leroy et al. (2008)
	################################################################################################

	mcmillan2011 = load(MCMILLAN2011_DATA_PATH)["dataframe"]
    leroy2008 = load(LEROY2008_DATA_PATH)["dataframe"]

    # List of available galaxies in Leroy et al. (2008)
    leroy_galaxies = unique(leroy2008[!, "Name"])

	################################################################################################
	# Find the target galaxies and read its values
	################################################################################################

    line_elements = Vector{LineElement}(undef, length(galaxies))
    labels = Vector{String}(undef, length(galaxies))

    for (i, galaxy) in pairs(galaxies)

        if galaxy ∈ leroy_galaxies

            leroy_data = filter(:Name => isequal(galaxy), leroy2008)
            r  = leroy_data[!, "Rad"]
            Σs = leroy_data[!, "Sigma*"] ± leroy_data[!, "e_Sigma*"]

        elseif galaxy == "MW"

            r  = mcmillan2011[!, :r]
            Σs = mcmillan2011[!, :Sigma]

        else

            throw(ArgumentError("ppAgertz2021!: $(galaxy) is not a valid galaxy"))

        end

        # Set the correct scale and units for the x axis
        if x_log
            x_data = @. log10(ustrip(x_unit, r))
        else
            x_data = @. ustrip(x_unit, r)
        end

        # Set the correct scale and units for the x axis
        if y_log
            y_data = @. Measurements.value(log10(ustrip(y_unit, Σs)))
            y_error = @. Measurements.uncertainty(log10(ustrip(y_unit, Σs)))
        else
            y_data = @. Measurements.value(ustrip(y_unit, Σs))
            y_error = @. Measurements.uncertainty(ustrip(y_unit, Σs))
        end

        color = ring(colors, i)

        if error_band
            bp = band!(
                figure.current_axis.x,
                x_data,
                y_data .- y_error,
                y_data .+ y_error;
                color=(color, 0.3),
            )
        end

        lp = lines!(figure.current_axis.x, x_data, y_data; color, linestyle, linewidth)

        # Put the post processing elements at the back of the plot
        if error_band
            translate!(Accum, bp, 0, 0, -10)
        end
        translate!(Accum, lp, 0, 0, -9)

        line_elements[i] = LineElement(; color, linestyle, linewidth)
        labels[i] = "$(galaxy)"

    end

    return line_elements, labels

end
