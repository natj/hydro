#Gravity related stuff


#Additional source terms
function source_terms(hyd::data2d, dt)

    #self-gravity
    #gxflx, gyflx = selfgravity(hyd.rho, hyd.nx, hyd.ny, r32x, r32y)

    #constant gravity

    #x-dir
    dir = 1

    #+
    gxflx, gyflx = ygravity(hyd.rhop[:,:,dir], hyd.nx, hyd.ny)
    hyd.velxp[:,:,dir] += 0.5*gxflx
    hyd.velyp[:,:,dir] += 0.5*gyflx
    hyd.epsp[:,:,dir] += 0.5*((hyd.velxp[:,:,dir] .* gxflx) .+ (hyd.velyp[:,:,dir] .* gyflx))

    #-
    gxflx, gyflx = ygravity(hyd.rhom[:,:,dir], hyd.nx, hyd.ny)
    hyd.velxm[:,:,dir] += 0.5*gxflx
    hyd.velym[:,:,dir] += 0.5*gyflx
    hyd.epsm[:,:,dir] += 0.5*((hyd.velxm[:,:,dir] .* gxflx) .+ (hyd.velym[:,:,dir] .* gyflx))


    #y-dir
    dir = 2

    #+
    gxflx, gyflx = ygravity(hyd.rhop[:,:,dir], hyd.nx, hyd.ny)
    hyd.velxp[:,:,dir] += 0.5*gxflx
    hyd.velyp[:,:,dir] += 0.5*gyflx
    hyd.epsp[:,:,dir] += 0.5*((hyd.velxp[:,:,dir] .* gxflx) .+ (hyd.velyp[:,:,dir] .* gyflx))

    #-
    gxflx, gyflx = ygravity(hyd.rhom[:,:,dir], hyd.nx, hyd.ny)
    hyd.velxm[:,:,dir] += 0.5*gxflx
    hyd.velym[:,:,dir] += 0.5*gyflx
    hyd.epsm[:,:,dir] += 0.5*((hyd.velxm[:,:,dir] .* gxflx) .+ (hyd.velym[:,:,dir] .* gyflx))


    return hyd
end

#Constant y-dir gravity
function ygravity(rho::AbstractMatrix,
                  nx::Int64,
                  ny::Int64;
                  gx=0.0,
                  gy=-1.0)

    gxflux = gx .* rho
    gyflux = gy .* rho

    return gxflux, gyflux
end

#Calculate x(x^2 + y^2)^-3/2 and y(x^2 + y^2)^-3/2 kernel
#for convolution used in self-gravity
function r32kernel(x, y)

    Nx = length(x)
    Ny = length(y)
    dx = x[2] - x[1]
    dy = y[2] - y[1]

    r32x = zeros(2Ny-1, 2Nx-1)
    r32y = zeros(2Ny-1, 2Nx-1)

    r32q = zeros(Ny, Nx)

    #construct 1 quarter of the kernel
    for j = 0:(Ny-1), i = 0:(Nx-1)
        r32q[j+1, i+1] = ((i*dx)^2.0 + (j*dy)^2.0)^(-1.5)
    end
    r32q[1,1] = 0.0

    #x - axis component
    r32x[1:Ny, Nx:2Nx-1] = flipud(r32q)
    r32x[Ny:2Ny-1, Nx:2Nx-1] = r32q

    #multiply by x
    for j = 1:2Ny-1, i = Nx:2Nx-1
        r32x[j, i] = r32x[j, i]*dx*(i-Nx)
    end

    #flip negative values
    r32x[:, 1:Nx-1] = -1.0*fliplr(r32x[:, Nx+1:end])


    #y - axis component
    r32y[Ny:2Ny-1, 1:Nx] = fliplr(r32q)
    r32y[Ny:2Ny-1, Nx:2Nx-1] = r32q

    #multiply by y
    for j = Ny:2Ny-1, i = 1:2Nx-1
        r32y[j, i] = r32y[j, i]*dy*(j-Ny)
    end

    #flip negative values
    r32y[1:Ny-1, :] = -1.0*flipud(r32y[Ny+1:end, :])

    return r32x, r32y
end

#Self-gravity
function selfgravity(rho::AbstractMatrix,
                     nx::Int64,
                     ny::Int64,
                     r32x::AbstractMatrix,
                     r32y::AbstractMatrix)

    xr = nx:(2nx-1)
    yr = ny:(2ny-1)

    gxflux = conv2(rho, r32x)[xr, yr]
    gyflux = conv2(rho, r32y)[xr, yr]

    gxflux = gxflux .* rho
    gyflux = gyflux .* rho

    return gxflux, gyflux
end


