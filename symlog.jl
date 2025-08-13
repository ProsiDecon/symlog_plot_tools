
log_signed(x) = sign(x) * log10(abs(x) + eps())
function symlog(x; linthresh=1e-7, linscale=1.0)
    return sign.(x) .* (log10.(abs.(x) ./ linthresh .+ 1) .* linscale)
end


function round_to_nice(x; base=10, direction=:nearest)
    """
        round_to_nice(x; base=10, direction=:nearest)

        Rounds `x` to the nearest (or lower/upper) power of `base`.
        `direction` can be `:nearest`, `:up`, or `:down`.
    """
    
    if x == 0
        return 0.0
    end
    exponent = log(base, abs(x))
    rounded_exp = direction == :nearest ? round(exponent) :
                  direction == :down    ? floor(exponent) :
                  direction == :up      ? ceil(exponent)  :
                  error("Invalid direction")
    return sign(x) * base^rounded_exp
end


function symlog_ticks(bounds::Tuple{Real, Real};
                             threshold=1e-7,
                             min_exp_step=1)
    """
    symmetric_log_ticks(bounds::Tuple{Real, Real};
                        threshold=1e-2,
                        min_exp_step=1)

    Returns log-scale-like ticks between bounds, symmetric around zero,
    but excludes values smaller than `threshold` in magnitude.
    Includes zero only if bounds contain it.
    """
    a, b = bounds
    ticks = Float64[]

    if a > b
        a, b = b, a
    end

    log_thresh = log10(threshold)

    # Negative side
    if a < -threshold
        neg_min_exp = floor(log10(abs(a)))
        neg_max_exp = ceil(log_thresh)
        neg_ticks = -1.0 .* (10.0 .^ (neg_min_exp:min_exp_step:neg_max_exp))
        append!(ticks, neg_ticks)
    end

    # Zero (only if within bounds)
    if a < 0 && b > 0
        push!(ticks, 0.0)
    end

    # Positive side
    if b > threshold
        pos_min_exp = ceil(log_thresh)
        pos_max_exp = floor(log10(b))
        pos_ticks = 10.0 .^ (pos_min_exp:min_exp_step:pos_max_exp)
        append!(ticks, pos_ticks)
    end

    return sort(ticks)
end

function symlog_ticks(x::Vector;
                             threshold=1e-7,
                             min_exp_step=1)
    """
    symmetric_log_ticks(bounds::Tuple{Real, Real};
                        threshold=1e-2,
                        min_exp_step=1)

    Returns log-scale-like ticks between bounds, symmetric around zero,
    but excludes values smaller than `threshold` in magnitude.
    Includes zero only if bounds contain it.
    """
    a = round_to_nice(minimum(x); direction = :down)
    b = round_to_nice(maximum(x); direction = :up)
    ticks = Float64[]

    if a > b
        a, b = b, a
    end

    log_thresh = log10(threshold)

    # Negative side
    if a < -threshold
        neg_min_exp = floor(log10(abs(a)))
        neg_max_exp = ceil(log_thresh)
        neg_ticks = -1.0 .* (10.0 .^ (neg_max_exp:min_exp_step:neg_min_exp))
        append!(ticks, sort(neg_ticks))
    end

    # Zero (only if within bounds)
    if a < 0 && b > 0
        push!(ticks, 0.0)
    end

    # Positive side
    if b > threshold
        pos_min_exp = ceil(log_thresh)
        pos_max_exp = floor(log10(b))
        pos_ticks = 10.0 .^ (pos_min_exp:min_exp_step:pos_max_exp)
        append!(ticks, pos_ticks)
    end

    return sort(ticks)
end

"""
    This function plots a scatter-plot with symmetric log axes.
"""
function symlog_plot(x::Vector,
                        y::Vector;
                        group::Union{Nothing, Vector} = nothing,
                        groupcolors::Union{Nothing, Vector} = nothing,
                        symlog_axes::Symbol = :both,
                        threshold::Float64 = 1e-7,
                        size_tuple::Tuple{Int64,Int64} =  (800, 500),
                        fortyfive_line::Bool = true)
    
    @assert length(x) == length(y) "Input vectors x and y not of same length"
    if !isnothing(group)
        @assert length(group) == length(x) "Vector of grouping variable not of same length as input vectors"
    end

    if symlog_axes ∈ [:both, :x]
        # define the symlog ticks for, labels and coordinate positions of observations for x-axis
        x_ticks = round.(symlog_ticks(x; threshold = threshold); digits = convert(Int64, abs(log10(threshold))) +1)
        x_pos = symlog.(x_ticks; linthresh = threshold)
        x_ticks = string.(x_ticks)
        x_tick_input = (x_pos, x_ticks)

        x_plot = symlog.(x; linthresh = threshold)
    else # make linear ticks, labels and plots
        x_tick_input = :automatic   # for linear axes, use automatic ticking

        x_plot = x
    end

    if symlog_axes ∈ [:both, :y] 
        y_ticks = round.(symlog_ticks(y; threshold = threshold); digits = convert(Int64, abs(log10(threshold))) +1)
        y_pos = symlog.(y_ticks; linthresh = threshold)
        y_ticks = string.(y_ticks)
        y_tick_input = (y_pos, y_ticks)

        y = symlog.(y; linthresh = threshold)
    else
        y_tick_input = :automatic       # for linear axes, use automatic ticking

        y_plot = y
    end

    if isnothing(group) || length(unique(group)) == 1
        p = Plots.scatter(x_plot, y_plot, xticks=x_tick_input, yticks=y_tick_input, msize = .5, markershape = :+, markercolor = :blue, alpha = .3; size = size_tuple)
    else
        allgroups = unique(group)
        p = Plots.scatter(x_plot[group .== allgroups[1]], y_plot[group .== allgroups[1]], xticks=x_tick_input, yticks=y_tick_input, msize = .5, markershape = :+, alpha = .3, label = string(allgroups[1]); size = size_tuple)
        for nowgroup ∈ allgroups[2:end]
            p = Plots.scatter!(x_plot[group .== nowgroup], y_plot[group .== nowgroup],  xticks=x_tick_input, yticks=y_tick_input, msize = .5, markershape = :+, alpha = .3, label = string(nowgroup); size = size_tuple)
        end
    end
    
    if fortyfive_line
        p = Plots.plot!(x_plot, x_plot, alpha = .8, linecolor = :black, linestyle = :solid, linewidth = 1., label = "45-degree line")
    end

    return p
end



###
# a specific utility for a personal application
function symlog_plot(d_out::DataFrame,
                        x::Symbol,
                        y::Symbol;
                        variety::Symbol = :both,        ### XXX tbd replace with grouping variable
                        threshold::Float64 = 1e-7,
                        size_tuple::Tuple{Int64,Int64} =  (800, 500))
        
    x_ticks = round.(symlog_ticks(d_out[d_out.variety .== string(variety), x]; threshold = threshold); digits = convert(Int64, abs(log10(threshold)))+1)
    x_pos = symlog.(x_ticks; linthresh = threshold)
    x_ticks = string.(x_ticks)

    y_ticks = round.(symlog_ticks(d_out[d_out.variety .== string(variety), y]; threshold = threshold); digits = convert(Int64, abs(log10(threshold)))+1)
    y_pos = symlog.(y_ticks; linthresh = threshold)
    y_ticks = string.(y_ticks)

    p = Plots.scatter(symlog.(d_out[d_out.for_any .== false .&& d_out.variety .== string(variety), x]; linthresh = threshold), symlog.(d_out[d_out.for_any .== false .&& d_out.variety .== string(variety), y]; linthresh = threshold), xticks=(x_pos, x_ticks), yticks=(y_pos, y_ticks), msize = .5, markershape = :+, markercolor = :blue, alpha = .3, label = "domestic"; size = size_tuple)
    p = Plots.scatter!(symlog.(d_out[d_out.for_any .== true .&& d_out.variety .== string(variety), x]; linthresh = threshold), symlog.(d_out[d_out.for_any .== true .&& d_out.variety .== string(variety), y]; linthresh = threshold), msize = .5, markershape = :x, markercolor = :red, alpha = .3, label = "foreign")
    p = Plots.plot!(symlog.(d_out[d_out.for_any .== false .&& d_out.variety .== string(variety), x]; linthresh = threshold), symlog.(d_out[d_out.for_any .== false .&& d_out.variety .== string(variety), x]; linthresh = threshold), alpha = .8, linecolor = :black, linestyle = :solid, linewidth = 1., label = "45-degree line")
    
    return p
end
