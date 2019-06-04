using StencilLab, GPUifyLoops, OffsetArrays, PyPlot

n = 32
ϕ = OffsetArray(rand(n+2, n+2, n+2), 0:n+1, 0:n+1, 0:n+1)
no_flux!(ϕ)

threads = (4, 4, 4)
blocks = map(x->Int(n/x), threads)

@show threads blocks

dev = CPU()

Δt = 0.01

fig, axs = subplots()
imshow(view(ϕ, :, :, 1))
gcf()

for i = 1:10
    for i = 1:100
        @launch dev threads=threads blocks=blocks diffusion_kernel_naive!(ϕ, Δt)
        no_flux!(ϕ)
    end

    imshow(view(ϕ, :, :, 1))
    gcf()
end
