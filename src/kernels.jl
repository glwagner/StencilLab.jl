@inline OneToMinusTwoPlus(n) = Base.OneTo(-2+n)

function interioraxes(a)
    Base.@_inline_meta
    return map(OneToMinusTwoPlus, size(a))
end

InteriorCartesianIndices(a) = CartesianIndices(interioraxes(a))

const c = 1
const f = -1

# Some 3D niceties
const ccc = (c, c, c)
const fcc = (f, c, c)
const cfc = (c, f, c)
const ccf = (c, c, f)

function ∇²_kernel_naive!(q, ψ, L)

    cuindex() = CartesianIndex(
        (blockIdx().x - 1) * blockDim().x + threadIdx().x,
        (blockIdx().y - 1) * blockDim().y + threadIdx().y,
        (blockIdx().z - 1) * blockDim().z + threadIdx().z
    )

    @loop for I in (InteriorCartesianIndices(ψ); cuindex())
        @inbounds q[I] = ∇²(I, L, ψ)
    end

    return nothing
end

∇²_kernel_naive_ccc!(q, ψ) = ∇²_kernel_naive!(q, ψ, ccc)
∇²_kernel_naive_fcc!(q, ψ) = ∇²_kernel_naive!(q, ψ, fcc)
∇²_kernel_naive_cfc!(q, ψ) = ∇²_kernel_naive!(q, ψ, cfc)
∇²_kernel_naive_ccf!(q, ψ) = ∇²_kernel_naive!(q, ψ, ccf)

function diffusion_kernel_naive!(ϕ, Δt)

    cuindex() = CartesianIndex(
        (blockIdx().x - 1) * blockDim().x + threadIdx().x,
        (blockIdx().y - 1) * blockDim().y + threadIdx().y,
        (blockIdx().z - 1) * blockDim().z + threadIdx().z
    )

    @loop for I in (InteriorCartesianIndices(ϕ); cuindex())
        @inbounds ϕ[I] += Δt * ∇²(I, ccc, ϕ)
    end

    return nothing
end
