#
# Differential operators for regular grids
#

@inline δx_caa(i, j, k, u::Number) = u
@inline δy_aca(i, j, k, u::Number) = u
@inline δz_aac(i, j, k, u::Number) = u

@inline δx_faa(i, j, k, u::Number) = u
@inline δy_afa(i, j, k, u::Number) = u
@inline δz_aaf(i, j, k, u::Number) = u

@inline δx_caa(i, j, k, u::AbstractArray) = @inbounds u[i+1, j, k] - u[i, j, k]
@inline δy_aca(i, j, k, u::AbstractArray) = @inbounds u[i, j+1, k] - u[i, j, k]
@inline δz_aac(i, j, k, u::AbstractArray) = @inbounds u[i, j, k+1] - u[i, j, k]

@inline δx_faa(i, j, k, u::AbstractArray) = @inbounds u[i, j, k] - u[i-1, j, k]
@inline δy_afa(i, j, k, u::AbstractArray) = @inbounds u[i, j, k] - u[i, j-1, k]
@inline δz_aaf(i, j, k, u::AbstractArray) = @inbounds u[i, j, k] - u[i, j, k-1]

@inline δx_caa(i, j, k, F::Function, args...) = F(i+1, j, k, args...) - F(i, j, k, args...)
@inline δy_aca(i, j, k, F::Function, args...) = F(i, j+1, k, args...) - F(i, j, k, args...)
@inline δz_aac(i, j, k, F::Function, args...) = F(i, j, k+1, args...) - F(i, j, k, args...)

@inline δx_faa(i, j, k, F::Function, args...) = F(i, j, k, args...) - F(i-1, j, k, args...)
@inline δy_afa(i, j, k, F::Function, args...) = F(i, j, k, args...) - F(i, j-1, k, args...)
@inline δz_aaf(i, j, k, F::Function, args...) = F(i, j, k, args...) - F(i, j, k-1, args...)

@inline ▶x_faa(i, j, k, u::Number, args...) = u 
@inline ▶y_afa(i, j, k, u::Number, args...) = u 
@inline ▶z_aaf(i, j, k, u::Number, args...) = u 

@inline ▶x_caa(i, j, k, u::Number, args...) = u 
@inline ▶y_aca(i, j, k, u::Number, args...) = u 
@inline ▶z_aac(i, j, k, u::Number, args...) = u 

@inline ▶x_faa(i, j, k, u::AbstractArray, args...) = @inbounds (u[i, j, k] + u[i-1, j, k]) / 2
@inline ▶y_afa(i, j, k, u::AbstractArray, args...) = @inbounds (u[i, j, k] + u[i, j-1, k]) / 2
@inline ▶z_aaf(i, j, k, u::AbstractArray, args...) = @inbounds (u[i, j, k] + u[i, j, k-1]) / 2

@inline ▶x_caa(i, j, k, u::AbstractArray, args...) = @inbounds (u[i+1, j, k] + u[i, j, k]) / 2
@inline ▶y_aca(i, j, k, u::AbstractArray, args...) = @inbounds (u[i, j+1, k] + u[i, j, k]) / 2
@inline ▶z_aac(i, j, k, u::AbstractArray, args...) = @inbounds (u[i, j, k+1] + u[i, j, k]) / 2

@inline ▶x_faa(i, j, k, F::Function, args...) = (F(i, j, k, args...) + F(i-1, j, k, args...)) / 2
@inline ▶y_afa(i, j, k, F::Function, args...) = (F(i, j, k, args...) + F(i, j-1, k, args...)) / 2
@inline ▶z_aaf(i, j, k, F::Function, args...) = (F(i, j, k, args...) + F(i, j, k-1, args...)) / 2

@inline ▶x_caa(i, j, k, F::Function, args...) = (F(i+1, j, k, args...) + F(i, j, k, args...)) / 2
@inline ▶y_aca(i, j, k, F::Function, args...) = (F(i, j+1, k, args...) + F(i, j, k, args...)) / 2
@inline ▶z_aac(i, j, k, F::Function, args...) = (F(i, j, k+1, args...) + F(i, j, k, args...)) / 2
