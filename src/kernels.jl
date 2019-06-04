const c = 1
const f = -1

function ∇²_kernel_naive!(q, ψ, L)

    cuindex() = CartesianIndex(
        (blockIdx().x - 1) * blockDim().x + threadIdx().x,
        (blockIdx().y - 1) * blockDim().y + threadIdx().y,
        (blockIdx().z - 1) * blockDim().z + threadIdx().z
    )

    @loop for I in (CartesianIndices(ψ); cuindex())
        @inbounds q[I] = ∇²(I, L, ψ)
    end
    return nothing
end

∇²_kernel_naive_ccc!(q, ψ) = ∇²_kernel_naive!(q, ψ, (c, c, c))
∇²_kernel_naive_fcc!(q, ψ) = ∇²_kernel_naive!(q, ψ, (f, c, c))
∇²_kernel_naive_cfc!(q, ψ) = ∇²_kernel_naive!(q, ψ, (c, f, c))
∇²_kernel_naive_ccf!(q, ψ) = ∇²_kernel_naive!(q, ψ, (c, c, f))
