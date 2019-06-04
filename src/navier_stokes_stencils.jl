# 
# Stencils associated with advection operators.
#

# Velocity * tracer
@inline u_ϕ_fcc(i, j, k, u, ϕ) = @inbounds u[i, j, k] * ▶x_faa(i, j, k, ϕ)
@inline v_ϕ_cfc(i, j, k, u, ϕ) = @inbounds v[i, j, k] * ▶y_afa(i, j, k, ϕ)
@inline w_ϕ_ccf(i, j, k, u, ϕ) = @inbounds w[i, j, k] * ▶z_aaf(i, j, k, ϕ)

# Velocity * velocity: u*v at ccc and ffc
@inline u_v_ccc(i, j, k, u, v) = ▶x_caa(i, j, k, u) * ▶y_aca(i, j, k, v)
@inline u_v_ffc(i, j, k, u, v) = ▶x_faa(i, j, k, v) * ▶y_afa(i, j, k, u)

# Velocity * velocity: u*w at ccc and fcf
@inline u_w_ccc(i, j, k, u, w) = ▶x_caa(i, j, k, u) * ▶z_aac(i, j, k, w)
@inline u_w_fcf(i, j, k, u, w) = ▶x_faa(i, j, k, w) * ▶z_aaf(i, j, k, u)

# Velocity * velocity: v*w at ccc and cff
@inline v_w_ccc(i, j, k, v, w) = ▶y_aca(i, j, k, v) * ▶z_aac(i, j, k, w)
@inline v_w_cff(i, j, k, v, w) = ▶y_afa(i, j, k, w) * ▶z_aaf(i, j, k, v)

# Advection operators
@inline function ∇_U_u(i, j, k, u, v, w)
    return (
              δx_faa(i, j, k, ▶x_caa,  ², u)
            + δy_aca(i, j, k, u_v_ffc, u, v)
            + δz_aac(i, j, k, u_w_fcf, u, w)
           )
end

@inline function ∇_U_v(i, j, k, u, v, w)
    return (
              δx_caa(i, j, k, u_v_ffc, u, v)
            + δy_afa(i, j, k, ▶y_aca,  ², u)
            + δz_aac(i, j, k, v_w_fcf, v, w)
           )
end

@inline function ∇_U_w(i, j, k, u, v, w)
    return (
              δx_caa(i, j, k, u_w_fcf, u, w)
            + δy_aca(i, j, k, v_w_cff, v, w)
            + δz_aaf(i, j, k, ▶z_aaf,  ², w)
           )
end

@inline function ∇_U_ϕ(i, j, k, u, v, w)
    return (
            δx_caa(i, j, k, u_ϕ_fcc, u, ϕ)
            δy_aca(i, j, k, v_ϕ_cfc, v, ϕ)
            δz_aac(i, j, k, w_ϕ_ccf, w, ϕ)
           )
end
