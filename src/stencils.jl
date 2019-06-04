# 
# The simplest stencils are formed by single operators.
# Even more interesting stencils arise in operator composition.
# 

@inline δx²_caa(i, j, k, u, args...) = δx_caa(i, j, k, δx_faa, u, args...)
@inline δy²_aca(i, j, k, u, args...) = δy_aca(i, j, k, δy_afa, u, args...)
@inline δz²_aac(i, j, k, u, args...) = δz_aac(i, j, k, δz_aaf, u, args...)

@inline δx²_faa(i, j, k, u, args...) = δx_faa(i, j, k, δx_caa, u, args...)
@inline δy²_afa(i, j, k, u, args...) = δy_aca(i, j, k, δy_aca, u, args...)
@inline δz²_aaf(i, j, k, u, args...) = δz_aaf(i, j, k, δz_aac, u, args...)

@inline function ∇²_ccc(i, j, k, u, args...)
    return (  ∂x²_caa(i, j, k, u, args...)
            + ∂y²_aca(i, j, k, u, args...)
            + ∂z²_aac(i, j, k, u, args...)
    )
end

@inline function ∇²_fcc(i, j, k, u, args...)
    return (  ∂x²_faa(i, j, k, u, args...)
            + ∂y²_aca(i, j, k, u, args...)
            + ∂z²_aac(i, j, k, u, args...)
    )
end

@inline function ∇²_cfc(i, j, k, u, args...)
    return (  ∂x²_caa(i, j, k, u, args...)
            + ∂y²_afa(i, j, k, u, args...)
            + ∂z²_aac(i, j, k, u, args...)
    )
end

@inline function ∇²_ccf(i, j, k, u, args...)
    return (  ∂x²_caa(i, j, k, u, args...)
            + ∂y²_aca(i, j, k, u, args...)
            + ∂z²_aaf(i, j, k, u, args...)
    )
end

#
# Double interpolation
#

@inline ▶xy_fca(i, j, k, u, args...) = ▶x_faa(i, j, k, ▶y_aca, u, args...)
@inline ▶xy_cfa(i, j, k, u, args...) = ▶x_caa(i, j, k, ▶y_afa, u, args...)
@inline ▶xy_cca(i, j, k, u, args...) = ▶x_caa(i, j, k, ▶y_aca, u, args...)
@inline ▶xy_ffa(i, j, k, u, args...) = ▶x_faa(i, j, k, ▶y_afa, u, args...)

@inline ▶xz_fac(i, j, k, u, args...) = ▶x_faa(i, j, k, ▶z_aac, u, args...)
@inline ▶xz_caf(i, j, k, u, args...) = ▶x_caa(i, j, k, ▶z_aaf, u, args...)
@inline ▶xz_cac(i, j, k, u, args...) = ▶x_caa(i, j, k, ▶z_aac, u, args...)
@inline ▶xz_faf(i, j, k, u, args...) = ▶x_faa(i, j, k, ▶z_aaf, u, args...)

@inline ▶yz_afc(i, j, k, u, args...) = ▶y_faa(i, j, k, ▶z_aac, u, args...)
@inline ▶yz_acf(i, j, k, u, args...) = ▶y_caa(i, j, k, ▶z_aaf, u, args...)
@inline ▶yz_acc(i, j, k, u, args...) = ▶y_caa(i, j, k, ▶z_aac, u, args...)
@inline ▶yz_aff(i, j, k, u, args...) = ▶y_faa(i, j, k, ▶z_aaf, u, args...)

# squaring operator
@inline ²(i, j, k, u) = @inbounds u[i, j, k]^2

