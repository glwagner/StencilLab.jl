using StencilLab, GPUifyLoops, OffsetArrays, PyPlot, Printf

Δt = 0.01
 n = 32
dev = CPU()

ϕ = OffsetArray(rand(n+2, n+2, n+2), 0:n+1, 0:n+1, 0:n+1)
no_flux!(ϕ)

threads = (4, 4, 4)
blocks = map(x->Int(n/x), threads)
@show threads blocks

ϕ₀ = deepcopy(ϕ)
t = 0.0

for i = 1:200
    global t
    @launch dev threads=threads blocks=blocks diffusion_kernel_naive!(ϕ, Δt)
    no_flux!(ϕ)
    t += Δt
end

fig, axs = subplots(ncols=2)

sca(axs[1])
imshow(view(ϕ₀, :, :, 1))
axs[1].axis("off")
title(L"t = 0")

sca(axs[2])
imshow(view(ϕ, :, :, 1))
axs[2].axis("off")
title(@sprintf("\$ t = %.1f \$", t))
gcf()

savefig("diffusion.png", dpi=480)
