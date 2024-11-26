using Revise

includet("sequences.jl")
includet("systems.jl")


"""
Measures the phase accumulation over the given `WaitSequence` using a simple Ramsey
experiment by first doing an X π/2 pulse (which makes the state |-Y⟩), and then
measuring in the X basis (via a ±π / 2 pulse) to determine the amount the Bloch vector
has rotated by (if any).

This works only in the range [-π / 2, π / 2], as the expectation value is the sin() of
the accumulated angle, so ambiguous beyond that.

If `flip_analysis` is true, the measurement is performed in the -X direction instead,
and the result is negated (asin(⟨X⟩) = -asin(⟨-X⟩)).
"""
function measure_simple_ramsey(system::System, wait_sequence::WaitSequence; flip_analysis=false)
    final_phase = wait_sequence.final_phase_transform(flip_analysis ? π / 2 : -π / 2)
    pulses = Pulse[
        XYRotation(π / 2, 0);
        wait_sequence.pulses;
        XYRotation(π / 2, final_phase)
    ]
    
    x_expected = 1 - 2 * measure(system, pulses)
    asin(x_expected) * (-1)^flip_analysis
end

"""
Measures the phase accumulation over the given `WaitSequence` using a "four-point"
Ramsey experiment, where first an X π/2 pulse prepares the |-Y⟩ state, and after the
wait sequence, the expectation value is measured along all four ±X, ±Y directions.

This has several effects: First, it allows us to reconstruct the entire xy part of the
Bloch vector including signs, and as works across such the [-π, π] range. It is also
robust against decoherence (if the length of the Bloch vector reduces, this would cause
the simple Ramsey method to under-report the accumulated angle).

Furthermore, using both the positive and negative signs helps reject bias in the
readout: For instance, if |0⟩ always correctly results in the 0 outcome, but |1⟩ only
95% of the time results in the 1 outcome and sometimes as 0 instead, this would bias a
simple Ramsey measurement, but that effect is averaged out in this measurement (such
effects are very common in physical systems, but not currently modelled in BasicSystem;
we should add them).

Lastly, using four different bases will certainly change how the measurement responds to
coherent pulse angle errors, but it doesn't seem too obvious how (which is part of what
we want to figure out here).
"""
function measure_four_point_ramsey(system::System, wait_sequence::WaitSequence)
    function meas(final_phase)
        pulses = Pulse[
            XYRotation(π / 2, 0);
            wait_sequence.pulses;
            XYRotation(π / 2, wait_sequence.final_phase_transform(final_phase))
        ]
        measure(system, pulses)
    end
    x_expected = meas(π / 2) - meas(-π / 2)
    y_expected = meas(π) - meas(0)
    atan(x_expected, -y_expected)
end
