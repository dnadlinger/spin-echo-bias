using Revise

includet("pulses.jl")

"""
A sequence of pulses to spend time accumulating a perturbation, along with a potential
transformation of the measurement phases in a "wrapping" Ramsey experiment to achieve
the same measurement outcome as in a straightforward Ramsey experiment.
"""
struct WaitSequence
    pulses::Vector{Pulse}
    final_phase_transform::Function
end

"""
Constructs a trivial wait sequence consisting of a single delay of the given time.
"""
function trivial_wait_sequence(perturbation_duration::Float64)::WaitSequence
    WaitSequence([Delay(perturbation_duration, true)], identity)
end

"""
Constructs a spin echo wait sequence with the given phase for the spin-echo pulse.
`on_in_first` determines whether the perturbation is on during the first or second
half.
"""
function spin_echo_wait_sequence(
    perturbation_duration::Float64,
    echo_phase::Float64;
    on_in_first::Bool=true
)::WaitSequence
    WaitSequence(
        [
            Delay(perturbation_duration, on_in_first),
            XYRotation(π, echo_phase),
            Delay(perturbation_duration, !on_in_first)
        ],
        ϕ -> 2echo_phase + π + (-1)^on_in_first * ϕ
    )
end

"""
Constructs a "multi-spin echo" wait sequence consisting of multiple π pulses with the
given phases, where the total time spent with the perturbation on is given by
`total_perturbation_duration`.

The perturbation is always on during alternating intervals; `on_in_first` determines
in which half it is.

For an even number of π pulses, the number of intervals between the bracketing π / 2
pulses is odd. As such, the first/last delays are made half the duration to keep the
overall time spent "in each parity" constant.
"""
function multi_spin_echo_wait_sequence(
    total_perturbation_duration::Float64,
    echo_phases::AbstractVector{Float64};
    on_in_first::Bool=true
)::WaitSequence
    perturbation_duration = total_perturbation_duration / ((length(echo_phases) + 1) ÷ 2)
    function echo_pulses(on_in_next)
        pulses = Pulse[XYRotation(π, echo_phases[1])]
        for ϕ in echo_phases[2:end]
            push!(pulses, Delay(perturbation_duration, on_in_next))
            on_in_next = !on_in_next
            push!(pulses, XYRotation(π, ϕ))
        end
        pulses
    end
    pulses = if iseven(length(echo_phases))
        [
            Delay(perturbation_duration / 2, on_in_first);
            echo_pulses(!on_in_first);
            Delay(perturbation_duration / 2, on_in_first)
        ]
    else
        [
            Delay(perturbation_duration, on_in_first);
            echo_pulses(!on_in_first);
            Delay(perturbation_duration, !on_in_first)
        ]
    end

    # Track what state is expected at the end of the sequence if there is no
    # perturbation.
    null_phase = foldl(
        (prev, pulse) -> mod(2pulse - prev, 2π),
        echo_phases;
        init=-π/2
    )
    # TODO: Sign of even-length sequences correct?
    WaitSequence(pulses, ϕ -> null_phase + π / 2 +
        (-1)^(on_in_first ⊻ iseven(length(echo_phases))) * ϕ
    )
end


"""
Constructs a Ramsey sequence with the given wait sequence and final phase.

The final phase is the phase of the second π/2 pulse, which is the phase
of the qubit state at the end of the sequence.
"""
function ramsey_sequence(
    wait_sequence::AbstractVector{Pulse}, final_phase::Float64
)::Vector{Pulse}
    [
        XYRotation(π / 2, 0.0);
        wait_sequence;
        XYRotation(π / 2, final_phase)
    ]
end
