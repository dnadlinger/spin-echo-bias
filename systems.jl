using Revise
using StaticArrays

includet("pulses.jl")

"""
A description of a physical system.

Should implement `measure()` to give the measured spin flip probability for a given
ideal pulse sequence under the conditions described, i.e if |0⟩ is the initial state,
the probability of measuring |1⟩ (which is |⟨1|ψ⟩|² for the final state |ψ⟩).
"""
abstract type System end


# Some helpers for single-qubit Pauli matrices/unitaries.

const id = @SMatrix [
    1 0;
    0 1
]
const x = @SMatrix [
    0 1;
    1 0
]
const y = @SMatrix [
    0 -1im;
    1im 0
]

"""
The single-qubit unitary for a rotation by amount θ around the axis in the xy plane of
the Bloch sphere with polar angle ϕ.
"""
xy_rot(θ, ϕ) = cos(θ / 2) * id - 1im * sin(θ / 2) * (cos(ϕ) * x + sin(ϕ) * y)


# BasicSystem with straightforward unitaries.

"""
A simple system, where the description of all operations is unitary, and the gate pulses
are (apart from a configurable area error) ideal.
"""
@kwdef struct BasicSystem <: System
    """
    The strength of the perturbation to measure, in rad/unit time; active during some of
    the wait operations.
    """
    perturbation_strength::Float64

    """
    The strength of a constant (unexpected) detuning of the system, in rad/unit time,
    which leads to phase accumulation during all wait operations (irrespective of
    whether the perturbation is supposed to be on or not).
    """
    detuning_strength::Float64

    """
    A rabi frequency scale factor that leads to coherent area errors on all the pulses.
    """
    rabi_frequency_scale::Float64 = 1.0
end

function unitary(system::BasicSystem, delay::Delay)
    freq = system.detuning_strength
    if delay.perturbation_on
        freq += system.perturbation_strength
    end
    @SMatrix [
        1 0;
        0 exp(1im * freq * delay.duration)
    ]
end

function unitary(system::BasicSystem, rot::XYRotation)
    xy_rot(
        rot.θ * system.rabi_frequency_scale,
        rot.ϕ
    )
end

function measure(system::BasicSystem, pulses::Vector{Pulse})
    ψ = @SVector ComplexF64[1, 0]
    for pulse in pulses
        ψ = unitary(system, pulse) * ψ
    end
    # Normalise just for numerical stability to avoid values a few epsilon above 1.
    # (Ought to clip to 0:1 instead?)
    abs2(ψ[2]) / sum(abs2, ψ)
end

