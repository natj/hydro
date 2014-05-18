#Grid setup

type data1d
    x::Array{Float64,1} #cell centers
    xi::Array{Float64,1} #cell LEFT interfaces

    rho::Array{Float64,1}
    rhop::Array{Float64,1}
    rhom::Array{Float64,1}

    velx::Array{Float64,1}
    velxp::Array{Float64,1}
    velxm::Array{Float64,1}

    vely::Array{Float64,1}
    velyp::Array{Float64,1}
    velym::Array{Float64,1}

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

        velx = zeros(nzones)
        velxp = zeros(nzones)
        velxm = zeros(nzones)

        vely = zeros(nzones)
        velyp = zeros(nzones)
        velym = zeros(nzones)

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
function grid_setup(self::data2d, xmin, xmax, ymin, ymax)
    dx = (xmax - xmin) / (self.nx -self.g*2 - 1)
    dy = (ymax - ymin) / (self.ny -self.g*2 - 1)

    xmin = xmin - self.g*dx
    xmax = xmax + self.g*dx
    ymin = ymin - self.g*dy
    ymax = ymax + self.g*dy

    #x setup
    for i = 0:self.nx-1
        self.x[i+1] = xmin + i*dx
    end

    for i = 0:self.nx-1
        self.xi[i+1] = self.x[i+1] - 0.5dx
    end

    #y setup
    for i = 0:self.ny-1
        self.y[i+1] = ymin + i*dy
    end

    for i = 0:self.ny-1
        self.yi[i+1] = self.y[i+1] - 0.5dy
    end


    return self
end


#splice 2dim data into x-axis pencils
#(to save time we copy only the edge values)
function splice(hyd::data2d, j::Int64, dim::Int64)

    if dim == 1 #x-dim
        #println("x j=$j $(hyd.ny)")
        @assert j <= hyd.ny
        nn = hyd.nx
        xs = 1:nn
        ys = j
        hyd1 = data1d(nn)
        hyd1.x = hyd.x
        hyd1.xi = hyd.xi

        hyd1.velxp[:] = vec(hyd.velxp[ys, xs, dim])
        hyd1.velxm[:] = vec(hyd.velxm[ys, xs, dim])

        hyd1.velyp[:] = vec(hyd.velyp[ys, xs, dim])
        hyd1.velym[:] = vec(hyd.velym[ys, xs, dim])

        hyd1.qm[:, 2] = vec(hyd.qm[ys, xs, 2, dim])
        hyd1.qp[:, 2] = vec(hyd.qp[ys, xs, 2, dim])

        hyd1.qm[:, 3] = vec(hyd.qm[ys, xs, 3, dim])
        hyd1.qp[:, 3] = vec(hyd.qp[ys, xs, 3, dim])

    elseif dim == 2 #y-dim
        #println("y j=$j $(hyd.nx)")
        @assert j <= hyd.nx
        nn = hyd.ny
        xs = j
        ys = 1:nn
        hyd1 = data1d(nn)
        hyd1.x = hyd.y
        hyd1.xi = hyd.yi

        hyd1.velxp[:] = vec(hyd.velyp[ys, xs, dim])
        hyd1.velxm[:] = vec(hyd.velym[ys, xs, dim])

        hyd1.velyp[:] = vec(hyd.velxp[ys, xs, dim])
        hyd1.velym[:] = vec(hyd.velxm[ys, xs, dim])

        hyd1.qm[:, 2] = vec(hyd.qm[ys, xs, 3, dim])
        hyd1.qp[:, 2] = vec(hyd.qp[ys, xs, 3, dim])

        hyd1.qm[:, 3] = vec(hyd.qm[ys, xs, 2, dim])
        hyd1.qp[:, 3] = vec(hyd.qp[ys, xs, 2, dim])

    else
        error("dim = $dim is not valid")
    end


    hyd1.g = hyd.g

    hyd1.rhom[:] = vec(hyd.rhom[ys, xs, dim])
    hyd1.rhop[:] = vec(hyd.rhop[ys, xs, dim])

    hyd1.epsm[:] = vec(hyd.epsm[ys, xs, dim])
    hyd1.epsp[:] = vec(hyd.epsp[ys, xs, dim])

    hyd1.pressm[:] = vec(hyd.pressm[ys, xs, dim])
    hyd1.pressp[:] = vec(hyd.pressp[ys, xs, dim])

    hyd1.qm[:, 1] = vec(hyd.qm[ys, xs, 1, dim])
    hyd1.qp[:, 1] = vec(hyd.qp[ys, xs, 1, dim])

    hyd1.qm[:, 4] = vec(hyd.qm[ys, xs, 4, dim])
    hyd1.qp[:, 4] = vec(hyd.qp[ys, xs, 4, dim])

    return hyd1
end

