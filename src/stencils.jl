# 
# Composed operators
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
