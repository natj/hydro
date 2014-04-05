#Grid setup

type data
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
    g::Float64 #ghost cells

    function data(nzones::Int64)

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

        q = zeros(3,nzones)
        qp = zeros(3,nzones)
        qm = zeros(3,nzones)
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



#1-dim grid
function grid_setup(self, xmin, xmax)
    dx = (xmax - xmin) / (self.n -self.g*2 - 1)
    xmin = xmin - self.g*dx
    xmax = xmax + self.g*dx

    for i = 1:self.n
        self.x[i] = xmin + i*dx
    end

    for i = 1:self.n
        self.xi[i] = self.x[i] - 0.5dx
    end

    return self
end

#Shoctube initial data
function setup_ID(self)
    rchange = 0.5(self.x[self.n - self.g] - self.x[self.g + 1])

    rho1 = 1.0
    rho2 = 0.125
    press1 = 1.0
    press2 = 0.1

    for i = 1:self.n
        if self.x[i] < rchange
            self.rho[i] = rho1
            self.press[i] = press1
            self.eps[i] = press1/rho1/(gamma - 1.0)
            self.vel[i] = 0.0
        else
            self.rho[i] = rho2
            self.press[i] = press2
            self.eps[i] = press2/rho2/(gamma - 1.0)
            self.vel[i] = 0.0
        end
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
    hyd.rho[1:hyd.g] = hyd.rho[hyd.g+1]
    hyd.vel[1:hyd.g] = hyd.vel[hyd.g+1]
    hyd.eps[1:hyd.g] = hyd.eps[hyd.g+1]
    hyd.press[1:hyd.g] = hyd.press[hyd.g+1]

    hyd.rho[(hyd.n-hyd.g+1) : hyd.n] = hyd.rho[hyd.n-hyd.g]
    hyd.vel[(hyd.n-hyd.g+1) : hyd.n] = hyd.vel[hyd.n-hyd.g]
    hyd.eps[(hyd.n-hyd.g+1) : hyd.n] = hyd.eps[hyd.n-hyd.g]
    hyd.press[(hyd.n-hyd.g+1) : hyd.n] = hyd.press[hyd.n-hyd.g]

    return hyd
end



