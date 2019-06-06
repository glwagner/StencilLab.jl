@inline shift(I::CartesianIndex{N}, dim, direction) where N =
    I + CartesianIndex(ntuple(i->ifelse(i==dim, direction, 0), Val(N)))

δ(I, dim, dir, u::T) where T<:Number = -zero(T)
δ(I, dim, dir, u::AbstractArray) = dir * u[shift(I, dim, dir)] - dir * u[I]
δ(I, dim, dir, u::Function, args...) = dir * u(shift(I, dim, dir), args...) - dir * u(I, args...)

▶(I, dim, dir, u::Number) = u
▶(I, dim, dir, u::AbstractArray) = (u[I] + u[shift(I, dim, dir)]) / 2
▶(I, dim, dir, u::Function, args...) = (u(I, args...) + u(shift(I, dim, dir), args...)) / 2

δ²(I, dim, dir, u, args...) = δ(I, dim, -dir, δ, dim, dir, u, args...)

∇²(I, L, u) = sum(δ²(I, dim, loc, u) for (dim, loc) in enumerate(L))

"""
    ▶(I, to, from, u)

Interpolate `u` from `from` to `to` at `I`.
"""
▶(I, to, from, u) = u[I] + (to - from) * δ(I, to, u) / 4


"""
    a_b(I, Lab, a, La, b, Lb)

Interpolate `a` and `b` from `La` and `Lb` to `Lab`
and multiply them there.
"""
a_b(I, Lab, La, a, Lb, b) = ▶(I, Lab, La, a) * ▶(I, Lab, Lb, b)
