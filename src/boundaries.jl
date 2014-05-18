#Different boundary conditions

#Standard continuous boundaries (allows outflow)
function apply_bcs(hyd::data2d)

    #arrays starting from zero
    #       |g                  |n-g #
    #[0 1 2 x x x  .....  x x x 7 8 9]

    #arrays starting from 1
    #     |g                  |n-g    #
    #[1 2 3 x x x  .....  x x x 8 9 10]

    #setup x boundaries
    for j = 1:hyd.g
        hyd.rho[:, j] = hyd.rho[:, hyd.g+1]
        hyd.velx[:, j] = hyd.velx[:, hyd.g+1]
        hyd.vely[:, j] = hyd.vely[:, hyd.g+1]
        hyd.eps[:, j] = hyd.eps[:, hyd.g+1]
        hyd.press[:, j] = hyd.press[:, hyd.g+1]
    end

    for j = (hyd.nx-hyd.g+1) : hyd.nx
        hyd.rho[:, j] = hyd.rho[:, hyd.nx-hyd.g]
        hyd.velx[:, j] = hyd.velx[:, hyd.nx-hyd.g]
        hyd.vely[:, j] = hyd.vely[:, hyd.nx-hyd.g]
        hyd.eps[:, j] = hyd.eps[:, hyd.nx-hyd.g]
        hyd.press[:, j] = hyd.press[:, hyd.nx-hyd.g]
    end

    #y boundaries
    for j = 1:hyd.g
        hyd.rho[j, :] = hyd.rho[hyd.g+1, :]
        hyd.velx[j, :] = hyd.velx[hyd.g+1, :]
        hyd.vely[j, :] = hyd.vely[hyd.g+1, :]
        hyd.eps[j, :] = hyd.eps[hyd.g+1, :]
        hyd.press[j, :] = hyd.press[hyd.g+1, :]
    end

    for j = (hyd.ny-hyd.g+1) : hyd.ny
        hyd.rho[j, :] = hyd.rho[hyd.ny-hyd.g, :]
        hyd.velx[j, :] = hyd.velx[hyd.ny-hyd.g, :]
        hyd.vely[j, :] = hyd.vely[hyd.ny-hyd.g, :]
        hyd.eps[j, :] = hyd.eps[hyd.ny-hyd.g, :]
        hyd.press[j, :] = hyd.press[hyd.ny-hyd.g, :]
    end

    #TODO: check corner boundaries

    return hyd
end

#reflective boundaries
function apply_refl_bcs(hyd::data2d)

    #arrays starting from zero
    #       |g                  |n-g #
    #[0 1 2 x x x  .....  x x x 7 8 9]

    #arrays starting from 1
    #     |g                  |n-g    #
    #[1 2 3 x x x  .....  x x x 8 9 10]

    #setup x boundaries
    for j = 1:hyd.g
        hyd.rho[:, j] = hyd.rho[:, hyd.g+1]
        hyd.velx[:, j] = -hyd.velx[:, hyd.g+1]
        hyd.vely[:, j] = -hyd.vely[:, hyd.g+1]
        hyd.eps[:, j] = hyd.eps[:, hyd.g+1]
        hyd.press[:, j] = hyd.press[:, hyd.g+1]
    end

    for j = (hyd.nx-hyd.g+1) : hyd.nx
        hyd.rho[:, j] = hyd.rho[:, hyd.nx-hyd.g]
        hyd.velx[:, j] = -hyd.velx[:, hyd.nx-hyd.g]
        hyd.vely[:, j] = -hyd.vely[:, hyd.nx-hyd.g]
        hyd.eps[:, j] = hyd.eps[:, hyd.nx-hyd.g]
        hyd.press[:, j] = hyd.press[:, hyd.nx-hyd.g]
    end

    #y boundaries
    for j = 1:hyd.g
        hyd.rho[j, :] = hyd.rho[hyd.g+1, :]
        hyd.velx[j, :] = -hyd.velx[hyd.g+1, :]
        hyd.vely[j, :] = -hyd.vely[hyd.g+1, :]
        hyd.eps[j, :] = hyd.eps[hyd.g+1, :]
        hyd.press[j, :] = hyd.press[hyd.g+1, :]
    end

    for j = (hyd.ny-hyd.g+1) : hyd.ny
        hyd.rho[j, :] = hyd.rho[hyd.ny-hyd.g, :]
        hyd.velx[j, :] = -hyd.velx[hyd.ny-hyd.g, :]
        hyd.vely[j, :] = -hyd.vely[hyd.ny-hyd.g, :]
        hyd.eps[j, :] = hyd.eps[hyd.ny-hyd.g, :]
        hyd.press[j, :] = hyd.press[hyd.ny-hyd.g, :]
    end

    #TODO: check corner boundaries

    return hyd
end

#periodic boundaries (toroidal)
function apply_per_bcs(hyd::data2d)

    #arrays starting from zero
    #       |g                  |n-g #
    #[0 1 2 x x x  .....  x x x 7 8 9]

    #arrays starting from 1
    #     |g                  |n-g    #
    #[1 2 3 x x x  .....  x x x 8 9 10]

    #setup x boundaries
    for j = 1:hyd.g
        hyd.rho[:, j] = hyd.rho[:, hyd.nx-2hyd.g+j-1]
        hyd.velx[:, j] = hyd.velx[:, hyd.nx-2hyd.g+j-1]
        hyd.vely[:, j] = hyd.vely[:, hyd.nx-2hyd.g+j-1]
        hyd.eps[:, j] = hyd.eps[:, hyd.nx-2hyd.g+j-1]
        hyd.press[:, j] = hyd.press[:, hyd.nx-2hyd.g+j-1]
    end

    i=1
    for j = (hyd.nx-hyd.g+1) : hyd.nx
        hyd.rho[:, j] = hyd.rho[:, hyd.g+i]
        hyd.velx[:, j] = hyd.velx[:, hyd.g+i]
        hyd.vely[:, j] = hyd.vely[:, hyd.g+i]
        hyd.eps[:, j] = hyd.eps[:, hyd.g+i]
        hyd.press[:, j] = hyd.press[:, hyd.g+i]
        i+=1
    end

    #y boundaries
    for j = 1:hyd.g
        hyd.rho[j, :] = hyd.rho[hyd.ny-2hyd.g+j-1, :]
        hyd.velx[j, :] = hyd.velx[hyd.ny-2hyd.g+j-1, :]
        hyd.vely[j, :] = hyd.vely[hyd.ny-2hyd.g+j-1, :]
        hyd.eps[j, :] = hyd.eps[hyd.ny-2hyd.g+j-1, :]
        hyd.press[j, :] = hyd.press[hyd.ny-2hyd.g+j-1, :]
    end

    i=1
    for j = (hyd.ny-hyd.g+1) : hyd.ny
        hyd.rho[j, :] = hyd.rho[hyd.g+i, :]
        hyd.velx[j, :] = hyd.velx[hyd.g+i, :]
        hyd.vely[j, :] = hyd.vely[hyd.g+i, :]
        hyd.eps[j, :] = hyd.eps[hyd.g+i, :]
        hyd.press[j, :] = hyd.press[hyd.g+i, :]
        i += 1
    end

    #TODO: check corner boundaries

    return hyd
end
