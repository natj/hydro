#Grid setup

type data2d
    x::Array{Float64, 1} #cell centers for X
    y::Array{Float64, 1} #cell centers for Y
    xi::Array{Float64, 1} #cell LEFT interfaces
    yi::Array{Float64, 1} #cell LOWER interfaces

    rho::Array{Float64, 2}
    rhop::Array{Float64, 2}
    rhom::Array{Float64, 2}

    vel::Array{Float64, 2}
    velp::Array{Float64, 2}
    velm::Array{Float64, 2}

    eps::Array{Float64, 2}
    epsp::Array{Float64, 2}
    epsm::Array{Float64, 2}

    press::Array{Float64, 2}
    pressp::Array{Float64, 2}
    pressm::Array{Float64, 2}

    q::Array{Float64, 3} #conserved quantities
    qp::Array{Float64, 3}
    qm::Array{Float64, 3}

    nx::Int64
    ny::Int64
    g::Float64 #ghost cells

    function data2d(nx::Int64, ny::Int64)

        x = zeros(nx)
        y = zeros(ny)
        xi = zeros(nx)
        yi = zeros(ny)

        rho = zeros(ny, nx)
        rhop = zeros(ny, nx)
        rhom = zeros(ny, nx)

        vel = zeros(ny, nx)
        velp = zeros(ny, nx)
        velm = zeros(ny, nx)

        eps = zeros(ny, nx)
        epsp = zeros(ny, nx)
        epsm = zeros(ny, nx)

        press = zeros(ny, nx)
        pressp = zeros(ny, nx)
        pressm = zeros(ny, nx)

        q = zeros(3, ny, nx)
        qp = zeros(3, ny, nx)
        qm = zeros(3, ny, nx)
        g = 3

        new(x,
            y,
            xi,
            yi,
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
            nx,
            ny,
            g)
    end
end



#1-dim grid
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

    rho1 = 0.4
    press1 = 0.2
    press2 = 0.99

    self.rho[:,:] = rho1*ones(self.ny, self.nx)
    self.press[:,:] = press1*ones(self.ny, self.nx)
    self.eps[:,:] = press1./rho1./(gamma - 1.0)
    self.vel[:,:] = zeros(self.ny, self.nx)

    midx = int(self.nx/2)
    midy = int(self.ny/2)
    self.press[midy-2:midy+2, midx-2:midx+2] = press2
    self.eps[midy-2:midy+2, midx-2:midx+2] = press2./rho1./(gamma - 1.0)

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
    hyd.rho[:, 1:hyd.g] = hyd.rho[:, hyd.g+1]
    hyd.vel[:, 1:hyd.g] = hyd.vel[:, hyd.g+1]
    hyd.eps[:, 1:hyd.g] = hyd.eps[:, hyd.g+1]
    hyd.press[:, 1:hyd.g] = hyd.press[:, hyd.g+1]

    hyd.rho[:, (hyd.nx-hyd.g+1) : hyd.nx] = hyd.rho[:, hyd.nx-hyd.g]
    hyd.vel[:, (hyd.nx-hyd.g+1) : hyd.nx] = hyd.vel[:, hyd.nx-hyd.g]
    hyd.eps[:, (hyd.nx-hyd.g+1) : hyd.nx] = hyd.eps[:, hyd.nx-hyd.g]
    hyd.press[:, (hyd.nx-hyd.g+1) : hyd.nx] = hyd.press[:, hyd.nx-hyd.g]

    #y boundaries
    hyd.rho[1:hyd.g, :] = hyd.rho[hyd.g+1, :]
    hyd.vel[1:hyd.g, :] = hyd.vel[hyd.g+1, :]
    hyd.eps[1:hyd.g, :] = hyd.eps[hyd.g+1, :]
    hyd.press[1:hyd.g, :] = hyd.press[hyd.g+1, :]

    hyd.rho[(hyd.nx-hyd.g+1) : hyd.nx, :] = hyd.rho[hyd.nx-hyd.g, :]
    hyd.vel[(hyd.nx-hyd.g+1) : hyd.nx, :] = hyd.vel[hyd.nx-hyd.g, :]
    hyd.eps[(hyd.nx-hyd.g+1) : hyd.nx, :] = hyd.eps[hyd.nx-hyd.g, :]
    hyd.press[(hyd.nx-hyd.g+1) : hyd.nx, :] = hyd.press[hyd.nx-hyd.g, :]

    #TODO: check corner boundaries

    return hyd
end



