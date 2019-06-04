module StencilLab

export # kernels
    ∇²_kernel_naive!,
    diffusion_kernel_naive!,

    # physics convenience
    no_flux!

using
    GPUifyLoops,
    OffsetArrays

# Import CUDA utilities if cuda is detected.
const HAVE_CUDA = try
    using CUDAdrv, CUDAnative, CuArrays
    true
catch
    false
end

macro hascuda(ex)
    return HAVE_CUDA ? :($(esc(ex))) : :(nothing)
end

@hascuda CuArrays, CUDAnative

include("operators.jl")
include("kernels.jl")
include("physics.jl")

end # module
