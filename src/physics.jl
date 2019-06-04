function no_flux!(u)

    @views @. u.parent[1, :, :]   = u.parent[2, :, :]
    @views @. u.parent[:, 1, :]   = u.parent[:, 2, :]
    @views @. u.parent[:, :, 1]   = u.parent[:, :, 2]

    @views @. u.parent[end, :, :] = u.parent[end-1, :, :]
    @views @. u.parent[:, end, :] = u.parent[:, end-1, :]
    @views @. u.parent[:, :, end] = u.parent[:, :, end-1]

    return nothing
end
