function ∇²_kernel_naive!(q, ψ, N)
    @loop for k in (1:N.z; (blockIdx().z - 1) * blockDim().z + threadIdx().z)
        @loop for j in (1:N.y; (blockIdx().y - 1) * blockDim().y + threadIdx().y)
            @loop for i in (1:N.x; (blockIdx().x - 1) * blockDim().x + threadIdx().x)

                @inbounds q[i, j, k] = ∇²_ccc(i, j, k, ψ)

            end
        end
    end
    return nothing
end
