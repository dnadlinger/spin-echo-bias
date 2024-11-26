"""
Generates lists of tick locations and tick labels for an axis labelled in (fractional)
multiples of π, e.g. ([-π / 2, 0.0, π / 2], ["-π / 2", "0", "π / 2"])
"""
function pi_axis_labels(min, max, step = 1//2)
    first, last = round.([min, max] / (pi * step))
    locs = Float64[]
    labels = String[]
    for k in first:last
        factor = Integer(k) * step
        push!(locs, factor * π)

        sign = factor < 0 ? "-" : ""

        abs_numerator = abs(numerator(factor))
        numerator_string =
            if abs_numerator == 0
                "0"
            elseif abs_numerator == 1
                "π"
            else
                "$abs_numerator π"
            end

        if denominator(factor) == 1
            push!(labels, "$sign $numerator_string")
        else
            push!(labels, "$sign $numerator_string / $(denominator(factor))")
        end
    end

    locs, labels
end

"""
Places ticks and labels on both axes in the given range, spaced by multiples of π.
"""
function set_pi_axis_labels!(axis, min, max, step = 1//2)
    locs, labels = pi_axis_labels(min, max, step)
    axis.xticks = (locs, labels)
    axis.yticks = (locs, labels)
end
