#Grid setup

type data1d
    x::Array{Float64,1} #cell centers
    xi::Array{Float64,1} #cell LEFT interfaces

    rho::Array{Float64,1}
    rhop::Array{Float64,1}
    rhom::Array{Float64,1}

    vel::Array{Float64,1}
    velp::Array{Float64,1}
    velm::Array{Float64,1}

    eps::Array{Float64,1}
    epsp::Array{Float64,1}
    epsm::Array{Float64,1}

    press::Array{Float64,1}
    pressp::Array{Float64,1}
    pressm::Array{Float64,1}

    q::Array{Float64,2} #conserved quantities
    qp::Array{Float64,2}
    qm::Array{Float64,2}

    n::Int64
    g::Int64 #ghost cells

    function data1d(nzones::Int64)

        x = zeros(nzones)
        xi = zeros(nzones)

        rho = zeros(nzones)
        rhop = zeros(nzones)
        rhom = zeros(nzones)

        vel = zeros(nzones)
        velp = zeros(nzones)
        velm = zeros(nzones)

        eps = zeros(nzones)
        epsp = zeros(nzones)
        epsm = zeros(nzones)

        press = zeros(nzones)
        pressp = zeros(nzones)
        pressm = zeros(nzones)

        q = zeros(nzones,3)
        qp = zeros(nzones,3)
        qm = zeros(nzones,3)
        n = nzones
        g = 3

        new(x,
            xi,
            rho,
            rhop,
            rhom,
            vel,
            velp,
            velm,
            eps,
            epsp,
            epsm,
            press,
            pressp,
            pressm,
            q,
            qp,
            qm,
            n,
            g)
    end
end


type data2d
    x::Array{Float64, 1} #cell centers for X
    y::Array{Float64, 1} #cell centers for Y
    xi::Array{Float64, 1} #cell LEFT interfaces
    yi::Array{Float64, 1} #cell LOWER interfaces

    rho::Array{Float64, 2}
    rhop::Array{Float64, 3}
    rhom::Array{Float64, 3}

    velx::Array{Float64, 2}
    velxp::Array{Float64, 3}
    velxm::Array{Float64, 3}
    vely::Array{Float64, 2}
    velyp::Array{Float64, 3}
    velym::Array{Float64, 3}

    eps::Array{Float64, 2}
    epsp::Array{Float64, 3}
    epsm::Array{Float64, 3}

    press::Array{Float64, 2}
    pressp::Array{Float64, 3}
    pressm::Array{Float64, 3}

    q::Array{Float64, 3} #conserved quantities
    qp::Array{Float64, 4}
    qm::Array{Float64, 4}

    nx::Int64
    ny::Int64
    g::Int64 #ghost cells

    function data2d(nx::Int64, ny::Int64)

        x = zeros(nx)
        y = zeros(ny)
        xi = zeros(nx)
        yi = zeros(ny)

        rho = zeros(ny, nx)
        rhop = zeros(ny, nx, 2)
        rhom = zeros(ny, nx, 2)

        velx = zeros(ny, nx)
        velxp = zeros(ny, nx, 2)
        velxm = zeros(ny, nx, 2)
        vely = zeros(ny, nx)
        velyp = zeros(ny, nx, 2)
        velym = zeros(ny, nx, 2)

        eps = zeros(ny, nx)
        epsp = zeros(ny, nx, 2)
        epsm = zeros(ny, nx, 2)

        press = zeros(ny, nx)
        pressp = zeros(ny, nx, 2)
        pressm = zeros(ny, nx, 2)

        q = zeros(ny, nx, 4)
        qp = zeros(ny, nx, 4, 2)
        qm = zeros(ny, nx, 4, 2)
        g = 3

        new(x,
            y,
            xi,
            yi,
            rho,
            rhop,
            rhom,
            velx,
            velxp,
            velxm,
            vely,
            velyp,
            velym,
            eps,
            epsp,
            epsm,
            press,
            pressp,
            pressm,
            q,
            qp,
            qm,
            nx,
            ny,
            g)
    end
end



#2-dim grid
function grid_setup(self, xmin, xmax, ymin, ymax)
    dx = (xmax - xmin) / (self.nx -self.g*2 - 1)
    dy = (ymax - ymin) / (self.ny -self.g*2 - 1)

    xmin = xmin - self.g*dx
    xmax = xmax + self.g*dx
    ymin = ymin - self.g*dy
    ymax = ymax + self.g*dy

    #x setup
    for i = 1:self.nx
        self.x[i] = xmin + i*dx
    end

    for i = 1:self.nx
        self.xi[i] = self.x[i] - 0.5dx
    end

    #y setup
    for i = 1:self.ny
        self.y[i] = ymin + i*dy
    end

    for i = 1:self.ny
        self.yi[i] = self.y[i] - 0.5dy
    end


    return self
end

#Shoctube initial data
#function setup_tube(self)
#    rchange = 0.5(self.x[self.n - self.g] - self.x[self.g + 1])

#    rho1 = 1.0
#    rho2 = 0.125
#    press1 = 1.0
#    press2 = 0.1

#    for i = 1:self.n
#        if self.x[i] < rchange
#            self.rho[i] = rho1
#            self.press[i] = press1
#            self.eps[i] = press1/rho1/(gamma - 1.0)
#            self.vel[i] = 0.0
#        else
#            self.rho[i] = rho2
#            self.press[i] = press2
#            self.eps[i] = press2/rho2/(gamma - 1.0)
#            self.vel[i] = 0.0
#        end
#    end

#    return self
#end


#Shoctube initial data
function setup_blast(self)
    #rchange = 0.5(self.x[self.n - self.g] - self.x[self.g + 1])

    rho1 = 0.1
    press1 = 0.1
    press2 = 1.0

    self.rho[:,:] = rho1*ones(self.ny, self.nx)
    self.press[:,:] = press1*ones(self.ny, self.nx)
    self.eps[:,:] = press1./rho1./(gamma - 1.0)
    self.velx[:,:] = zeros(self.ny, self.nx)
    self.vely[:,:] = zeros(self.ny, self.nx)

    midx = int(self.nx/2)
    midy = int(self.ny/2)
    self.press[midy-2:midy+2, midx-2:midx+2] = press2
    self.eps[midy-2:midy+2, midx-2:midx+2] = press2./rho1./(gamma - 1.0)

    return self
end

#Shoctube initial data
function setup_taylor(self)
    rchange = int(0.8self.ny)

    rho1 = 0.8
    rho2 = 0.1
    press1 = 0.1

    #println("$(self.y[self.ny - self.g])")
    #println("$(self.y[self.g + 1])")

    #println("rchange= ", rchange)
    #println("$((rchange+1):self.ny)")

    self.rho[1:rchange, :] = rho2*ones(rchange, self.nx)
    self.rho[(rchange+1):(self.ny), :] = rho1*ones((self.ny-rchange), self.nx)

    self.press[:,:] = press1*ones(self.ny, self.nx)
    self.eps[:,:] = press1./rho1./(gamma - 1.0)

    self.velx[:,:] = zeros(self.ny, self.nx)
    self.vely[:,:] = zeros(self.ny, self.nx)

    dx = self.x[self.nx] - self.x[1]
    for j = rchange:self.ny
        self.vely[j,10:(self.nx-10-1)] = Float64[-0.5sin(pi*self.x[n]/dx) for n = 10:(self.nx-10-1)]
    end

    return self
end




function apply_bcs(hyd)

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

    for j = (hyd.nx-hyd.g+1) : hyd.nx
        hyd.rho[j, :] = hyd.rho[hyd.nx-hyd.g, :]
        hyd.velx[j, :] = hyd.velx[hyd.nx-hyd.g, :]
        hyd.vely[j, :] = hyd.vely[hyd.nx-hyd.g, :]
        hyd.eps[j, :] = hyd.eps[hyd.nx-hyd.g, :]
        hyd.press[j, :] = hyd.press[hyd.nx-hyd.g, :]
    end

    #TODO: check corner boundaries

    return hyd
end

#splice 2dim data into x-axis pencils
#to save time we copy only the edge values
function splice(hyd::data2d, j::Int64, dim::Int64)

    if dim == 1 #x-dim
        #println("x j=$j $(hyd.ny)")
        @assert j < hyd.ny
        nn = hyd.nx
        xs = 1:nn
        ys = j
        hyd1 = data1d(nn)
        hyd1.x = hyd.x
        hyd1.xi = hyd.xi

        hyd1.velp[:] = vec(hyd.velxp[ys, xs, dim])
        hyd1.velm[:] = vec(hyd.velxm[ys, xs, dim])

        hyd1.qm[:, 2] = vec(hyd.qm[ys, xs, 2, dim])
        hyd1.qp[:, 2] = vec(hyd.qp[ys, xs, 2, dim])

    elseif dim == 2 #y-dim
        #println("y j=$j $(hyd.nx)")
        @assert j < hyd.nx
        nn = hyd.ny
        xs = j
        ys = 1:nn
        hyd1 = data1d(nn)
        hyd1.x = hyd.y
        hyd1.xi = hyd.yi

        hyd1.velp[:] = vec(hyd.velyp[ys, xs, dim])
        hyd1.velm[:] = vec(hyd.velym[ys, xs, dim])

        hyd1.qm[:, 3] = vec(hyd.qm[ys, xs, 3, dim])
        hyd1.qp[:, 3] = vec(hyd.qp[ys, xs, 3, dim])
    else
        error("dim = $dim is not valid")
    end


    #hyd1.g = hyd.g
    #hyd1.rho[:] = vec(hyd.rho[i, :])
    hyd1.rhom[:] = vec(hyd.rhom[ys, xs, dim])
    hyd1.rhop[:] = vec(hyd.rhop[ys, xs, dim])

    hyd1.epsm[:] = vec(hyd.epsm[ys, xs, dim])
    hyd1.epsp[:] = vec(hyd.epsp[ys, xs, dim])

    hyd1.pressm[:] = vec(hyd.pressm[ys, xs, dim])
    hyd1.pressp[:] = vec(hyd.pressp[ys, xs, dim])

    hyd1.qm[:, 1] = vec(hyd.qm[ys, xs, 1, dim])
    hyd1.qp[:, 1] = vec(hyd.qp[ys, xs, 1, dim])
    hyd1.qm[:, 3] = vec(hyd.qm[ys, xs, 4, dim])
    hyd1.qp[:, 3] = vec(hyd.qp[ys, xs, 4, dim])

    return hyd1
end

