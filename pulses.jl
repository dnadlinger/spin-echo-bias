"""
One of a series of oferations on a two-level system.

"Pulse" is used in a very loose sense here, as this includes delays, etc.
"""
abstract type Pulse end

struct Delay <: Pulse
    """The duration of the delay.
    
    The amount is only meaningful in conjunction with some externa information giving
    the strength of a perturbation/static qubit detuning, etc., as there is ideally no
    time evolution during a wait time.
    """
    duration::Float64

    "Whether to switch on the perturbation to measure during this delay."
    perturbation_on::Bool
end

struct XYRotation <: Pulse
    "The (target) rotation angle (gate area), in radians."
    θ::Float64

    "The phase of the gate (angle of rotation axis in the xy plane), in radians."
    ϕ::Float64
end

