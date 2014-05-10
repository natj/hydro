#Grid setup

type data1d

    #positions
    x::Array{Float64,1} #cell centers
    xi::Array{Float64,1} #cell LEFT interfaces

    #density
    rho::Array{Float64,1}
    rhop::Array{Float64,1}
    rhom::Array{Float64,1}

    #normal velocity (e.g. x-dir)
    vel::Array{Float64,1}
    velp::Array{Float64,1}
    velm::Array{Float64,1}

    #transverse velocity (e.g. y-dir)
    #NOTE: in reality we have also second transverse velocity in z-dir
    #velt::Array{Float64,1}
    veltp::Array{Float64,1}
    veltm::Array{Float64,1}
    
    #internal energy
    eps::Array{Float64,1}
    epsp::Array{Float64,1}
    epsm::Array{Float64,1}

    #pressure
    press::Array{Float64,1}
    pressp::Array{Float64,1}
    pressm::Array{Float64,1}

    #conserved quantities
    q::Array{Float64,2}
    qp::Array{Float64,2}
    qm::Array{Float64,2}

    #other
    n::Int64 #grid points
    g::Int64 #ghost cells

    function data1d(nzones::Int64)

        #initialize all arrays

        x = zeros(nzones)
        xi = zeros(nzones)

        rho = zeros(nzones)
        rhop = zeros(nzones)
        rhom = zeros(nzones)

        vel = zeros(nzones)
        velp = zeros(nzones)
        velm = zeros(nzones)

        velt = zeros(nzones)
        veltp = zeros(nzones)
        veltm = zeros(nzones)

        eps = zeros(nzones)
        epsp = zeros(nzones)
        epsm = zeros(nzones)

        press = zeros(nzones)
        pressp = zeros(nzones)
        pressm = zeros(nzones)

        q = zeros(nzones, 4)
        qp = zeros(nzones, 4)
        qm = zeros(nzones, 4)
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
            veltp,
            veltm,
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

    #positions
    x::Array{Float64, 1} #cell centers for X
    y::Array{Float64, 1} #cell centers for Y
    xi::Array{Float64, 1} #cell LEFT interfaces
    yi::Array{Float64, 1} #cell LOWER interfaces

    #density
    rho::Array{Float64, 2}
    rhop::Array{Float64, 3}
    rhom::Array{Float64, 3}

    #velocities
    velx::Array{Float64, 2}
    velxp::Array{Float64, 3}
    velxm::Array{Float64, 3}
    vely::Array{Float64, 2}
    velyp::Array{Float64, 3}
    velym::Array{Float64, 3}

    #internal energy
    eps::Array{Float64, 2}
    epsp::Array{Float64, 3}
    epsm::Array{Float64, 3}

    #pressure
    press::Array{Float64, 2}
    pressp::Array{Float64, 3}
    pressm::Array{Float64, 3}

    #conserved quantities
    q::Array{Float64, 3}
    qp::Array{Float64, 4}
    qm::Array{Float64, 4}

    #other
    nx::Int64 #x-grid points
    ny::Int64 #y-grid points
    g::Int64 #ghost cells

    function data2d(nx::Int64, ny::Int64)

        #initialize all arrays
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
function grid_setup(self::data2d, xmin, xmax, ymin, ymax)
    dx = (xmax - xmin) / (self.nx -self.g*2 - 1)
    dy = (ymax - ymin) / (self.ny -self.g*2 - 1)

    xmin = xmin - self.g*dx
    xmax = xmax + self.g*dx
    ymin = ymin - self.g*dy
    ymax = ymax + self.g*dy

    #x-dir
    for i = 1:self.nx
        self.x[i] = xmin + i*dx
    end

    for i = 1:self.nx
        self.xi[i] = self.x[i] - 0.5dx
    end

    #y-dir
    for i = 1:self.ny
        self.y[i] = ymin + i*dy
    end

    for i = 1:self.ny
        self.yi[i] = self.y[i] - 0.5dy
    end


    return self
end

#Shoctube initial data
function setup_tubey(self::data2d)
    rchange = 0.5(self.y[self.ny - self.g] - self.y[self.g + 1])

    rho1 = 1.0
    rho2 = 0.125
    press1 = 1.0
    press2 = 0.1


    for i = 1:self.ny, j = 1:self.nx
        if self.y[i] < rchange
            self.rho[i, j] = rho1
            self.press[i ,j] = press1
            self.eps[i, j] = press1/rho1/(gamma - 1.0)
            self.velx[i, j] = 0.0
            self.vely[i, j] = 0.0
        else
            self.rho[i, j] = rho2
            self.press[i, j] = press2
            self.eps[i, j] = press2/rho2/(gamma - 1.0)
            self.velx[i, j] = 0.0
            self.vely[i, j] = 0.0
        end
    end

    return self
end

#Shoctube initial data
function setup_tubex(self::data2d)
    rchange = 0.5(self.x[self.nx - self.g] - self.x[self.g + 1])

    rho1 = 1.0
    rho2 = 0.125
    press1 = 1.0
    press2 = 0.1

    for i = 1:self.ny, j = 1:self.nx
        if self.x[j] < rchange
            self.rho[i, j] = rho1
            self.press[i ,j] = press1
            self.eps[i, j] = press1/rho1/(gamma - 1.0)
            self.velx[i, j] = 0.0
            self.vely[i, j] = 0.0
        else
            self.rho[i, j] = rho2
            self.press[i, j] = press2
            self.eps[i, j] = press2/rho2/(gamma - 1.0)
            self.velx[i, j] = 0.0
            self.vely[i, j] = 0.0
        end
    end

    return self
end

#Shoctube initial data
function setup_tubexy(self::data2d)

    @assert hyd.x[1] == hyd.y[1]
    @assert hyd.x[end] == hyd.y[end]
    @assert hyd.nx == hyd.ny

    rho1 = 1.0
    rho2 = 0.125
    press1 = 1.0
    press2 = 0.1

    self.rho[:,:] = rho1
    self.press[:,:] = press1
    self.eps[:,:] = press1/rho1/(gamma - 1.0)

    for i = 1:self.ny
        for j = 1:i
                self.rho[i, j] = rho2
                self.press[i, j] = press2
                self.eps[i, j] = press2/rho2/(gamma - 1.0)
            end
        end

    return self
end


#Shoctube initial data
function setup_blast(self::data2d)

    rho1 = 1.0
    press1 = 0.1
    press2 = 10.0

    self.rho[:,:] = rho1*ones(self.ny, self.nx)
    self.press[:,:] = press1*ones(self.ny, self.nx)
    self.eps[:,:] = press1./rho1./(gamma - 1.0)
    self.velx[:,:] = zeros(self.ny, self.nx)
    self.vely[:,:] = zeros(self.ny, self.nx)

    midx = int(self.nx/2)
    midy = int(self.ny/2)


    #circle=[1,3,4,4,5,5]
    #circle=[2,4,5,6,6,6,7,7]
    circle=[3,5,6,7,8,9,9,10,10,10]
    N = length(circle)
    for j = 0:N-1
        yy = j
        xx = circle[N-j]

        self.press[midy+yy, midx-xx:midx+xx] = press2
        self.eps[midy+yy, midx-xx:midx+xx] = press2./rho1./(gamma - 1.0)
        self.press[midy-yy, midx-xx:midx+xx] = press2
        self.eps[midy-yy, midx-xx:midx+xx] = press2./rho1./(gamma - 1.0)
    end

    return self
end

#Two colliding spheres
function setup_collision(self::data2d)

    rho1 = 0.01
    rho2 = 10.0
    press1 = 0.01
    vel2 = 2.0

    self.rho[:,:] = rho1*ones(self.ny, self.nx)
    self.press[:,:] = press1*ones(self.ny, self.nx)
    self.eps[:,:] = press1./rho1./(gamma - 1.0)
    self.velx[:,:] = zeros(self.ny, self.nx)
    self.vely[:,:] = zeros(self.ny, self.nx)

    #first ball
    midx = int(self.nx/2)
    midy = int(0.25self.ny)

    #circle=[1,3,4,4,5,5]
    #circle=[2,4,5,6,6,6,7,7]
    circle=[3,5,6,7,8,9,9,10,10,10]
    N = length(circle)
    for j = 0:N-1
        yy = j
        xx = circle[N-j]

        self.rho[midy+yy, midx-xx:midx+xx] = rho2
        self.rho[midy-yy, midx-xx:midx+xx] = rho2

        self.vely[midy+yy, midx-xx:midx+xx] = vel2
        self.vely[midy-yy, midx-xx:midx+xx] = vel2

    end

    #second ball
    midx = int(self.nx/2)
    midy = int(0.75self.ny)
    #circle=[1,3,4,4,5,5]
    #circle=[2,4,5,6,6,6,7,7]
    circle=[3,5,6,7,8,9,9,10,10,10]
    N = length(circle)
    for j = 0:N-1
        yy = j
        xx = circle[N-j]

        self.rho[midy+yy, midx-xx:midx+xx] = rho2
        self.rho[midy-yy, midx-xx:midx+xx] = rho2

        self.vely[midy+yy, midx-xx:midx+xx] = -vel2
        self.vely[midy-yy, midx-xx:midx+xx] = -vel2
    end

    return self
end

#One free falling spheres
function setup_fall(self::data2d)

    rho1 = 1.0
    rho2 = 2.0
    press1 = 2.5
    press2 = 1.25

    self.rho[:,:] = rho1*ones(self.ny, self.nx)
    self.press[:,:] = press1*ones(self.ny, self.nx)
    self.eps[:,:] = press1./rho1./(gamma - 1.0)
    self.velx[:,:] = zeros(self.ny, self.nx)
    self.vely[:,:] = zeros(self.ny, self.nx)

    #first ball
    midx = int(self.nx/2)
    midy = int(0.8self.ny)

    #circle=[1,3,4,4,5,5]
    #circle=[2,4,5,6,6,6,7,7]
    circle=[3,5,6,7,8,9,9,10,10,10]
    N = length(circle)
    for j = 0:N-1
        yy = j
        xx = circle[N-j]

        self.rho[midy+yy, midx-xx:midx+xx] = rho2
        self.rho[midy-yy, midx-xx:midx+xx] = rho2

        self.press[midy+yy, midx-xx:midx+xx] = press2
        self.press[midy-yy, midx-xx:midx+xx] = press2
    end

    return self
end


#Taylor instability initial data
function setup_taylor(self::data2d)
    rchange = int(0.8self.ny)

    grav = 0.005
    rho1 = 2.0
    rho2 = 1.0
    press1 = 0.01

    self.rho[1:(rchange-1), :] = rho2*ones(rchange-1, self.nx)
    self.rho[rchange:(self.ny), :] = rho1*ones((self.ny-rchange)+1, self.nx)

    self.press[:,:] = press1*ones(self.ny, self.nx)
    self.eps[:,:] = press1./rho1./(gamma - 1.0)

    self.velx[:,:] = zeros(self.ny, self.nx)
    self.vely[:,:] = zeros(self.ny, self.nx)


    offs = 35
    offsy= 30
    dx = self.x[self.nx-offs-1] - self.x[offs]

    function siny(j)
        sin(pi*(self.y[j]-self.y[rchange-offsy])/(self.y[rchange+offsy]-self.y[rchange-offsy]))
    end


#    for j = rchange:self.ny

    #for j = (rchange-offsy):(rchange+offsy)
    for j = (rchange):(rchange+offsy)

        self.vely[j, offs:(self.nx-offs-1)] = Float64[-0.01siny(j)*sin(pi*(self.x[n]-self.x[offs])/dx) for n = offs:(self.nx-offs-1)]
#        self.vely[j, offs:(self.nx-offs-1)] = Float64[-10.0sin(pi*(self.x[n]-self.x[offs])/dx) for n = offs:(self.nx-offs-1)]
    end

    for j = 1:self.ny
        self.press[j,:] = 2.5 - self.rho[j,20]*grav*self.y[j]
    end


    return self
end

#Taylor instability initial data according to Athena code
function setup_taylor2(self::data2d)
    rchange = int(0.5self.ny)

    A = 0.01
    grav = 0.005
    rho1 = 2.0
    rho2 = 1.0

    self.rho[1:(rchange-1), :] = rho2*ones(rchange-1, self.nx)
    self.rho[rchange:(self.ny), :] = rho1*ones((self.ny-rchange)+1, self.nx)

    self.press = zeros(self.ny, self.nx)

    self.velx[:,:] = zeros(self.ny, self.nx)
    self.vely[:,:] = zeros(self.ny, self.nx)

    Lx = abs(self.x[self.nx - self.g] - self.x[self.g+1])
    Ly = abs(self.y[self.ny - self.g] - self.y[self.g+1])

    for j = (self.g+1):(self.ny - self.g), i = (self.g+1):(self.nx - self.g)
        self.vely[j,i] = A*((1.0+cos(2pi*self.x[i]/Lx)*(1.0+cos(2pi*self.y[j]/Ly))))
    end

    #set pressure to hydrostatic equiblibrium
    for j = (self.g+1):(self.ny-self.g)
        self.press[j,:] = 2.5 - self.rho[j,20]*grav*self.y[j]
    end
    self.eps[:,:] = self.press[:,:] ./ self.rho[:,:] ./ (gamma - 1.0)

    return self
end



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

    #TODO: check corner boundaries; are they correct?
    #TODO: Do we need to force boundary condities to left and right states too?

    return hyd
end

#splice 2dim data into x-axis pencils;
#y-axis data is rotated into x-axis pencils, too
#
#to save time we copy only the edge values
function splice(hyd::data2d, j::Int64, dim::Int64)

    if dim == 1 #x-dim

        @assert j < hyd.ny
        nn = hyd.nx
        xs = 1:nn
        ys = j
        hyd1 = data1d(nn)
        hyd1.x = hyd.x
        hyd1.xi = hyd.xi

        #normal velocity
        hyd1.velp[:] = vec(hyd.velxp[ys, xs, dim])
        hyd1.velm[:] = vec(hyd.velxm[ys, xs, dim])

        #transverse velocity
        hyd1.veltp[:] = vec(hyd.velyp[ys, xs, dim])
        hyd1.veltm[:] = vec(hyd.velym[ys, xs, dim])

        #normal momentum
        hyd1.qm[:, 2] = vec(hyd.qm[ys, xs, 2, dim])
        hyd1.qp[:, 2] = vec(hyd.qp[ys, xs, 2, dim])

        #transverse momentum
        hyd1.qm[:, 3] = vec(hyd.qm[ys, xs, 3, dim])
        hyd1.qp[:, 3] = vec(hyd.qp[ys, xs, 3, dim])

    elseif dim == 2 #y-dim
        
        @assert j < hyd.nx
        nn = hyd.ny
        xs = j
        ys = 1:nn
        hyd1 = data1d(nn)
        hyd1.x = hyd.y
        hyd1.xi = hyd.yi

        #normal velocity
        hyd1.velp[:] = vec(hyd.velyp[ys, xs, dim])
        hyd1.velm[:] = vec(hyd.velym[ys, xs, dim])

        #transverse velocity
        hyd1.veltp[:] = vec(hyd.velxp[ys, xs, dim])
        hyd1.veltm[:] = vec(hyd.velxm[ys, xs, dim])

        #normal momentum
        hyd1.qm[:, 2] = vec(hyd.qm[ys, xs, 3, dim])
        hyd1.qp[:, 2] = vec(hyd.qp[ys, xs, 3, dim])

        #transverse momentum
        hyd1.qm[:, 3] = vec(hyd.qm[ys, xs, 2, dim])
        hyd1.qp[:, 3] = vec(hyd.qp[ys, xs, 2, dim])

    else
        error("dim = $dim is not valid")
    end


    #hyd1.g = hyd.g
    #hyd1.rho[:] = vec(hyd.rho[i, :])

    #density
    hyd1.rhom[:] = vec(hyd.rhom[ys, xs, dim])
    hyd1.rhop[:] = vec(hyd.rhop[ys, xs, dim])

    #internal energy
    hyd1.epsm[:] = vec(hyd.epsm[ys, xs, dim])
    hyd1.epsp[:] = vec(hyd.epsp[ys, xs, dim])

    #pressure
    hyd1.pressm[:] = vec(hyd.pressm[ys, xs, dim])
    hyd1.pressp[:] = vec(hyd.pressp[ys, xs, dim])

    #conserved variables
    hyd1.qm[:, 1] = vec(hyd.qm[ys, xs, 1, dim])
    hyd1.qp[:, 1] = vec(hyd.qp[ys, xs, 1, dim])
    hyd1.qm[:, 4] = vec(hyd.qm[ys, xs, 4, dim])
    hyd1.qp[:, 4] = vec(hyd.qp[ys, xs, 4, dim])

    return hyd1
end

