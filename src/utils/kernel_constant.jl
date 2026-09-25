"""
    KernelConstant(value)

Explicitly carry a read-only scalar or isbits parameter bundle through a
KernelAbstractions launch. The Enzyme extension marks this wrapper inactive so
CUDA reverse mode does not reinterpret fixed floating-point kernel arguments
as `Active` values.
"""
struct KernelConstant{T}
    value::T
end

@inline kernel_constant(value) = KernelConstant(value)
@inline kernel_constant(reference, value) =
    _kernel_constant(KernelAbstractions.get_backend(reference), value)
@inline _kernel_constant(::KernelAbstractions.GPU, value) = KernelConstant(value)
@inline _kernel_constant(::KernelAbstractions.CPU, value) = value

@inline kernel_value(value) = value
@inline kernel_value(parameter::KernelConstant) = parameter.value

Adapt.@adapt_structure KernelConstant
