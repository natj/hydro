#Gravity related stuff

function ygravity(rho::AbstractMatrix,
                  nx::Int64,
                  ny::Int64;
                  gx=0.0,
                  gy=5.0)

    gxflux = gx .* rho
    gyflux = gy .* rho

    return gxflux, gyflux
end

#Calculate x(x^2 + y^2)^-3/2 and y(x^2 + y^2)^-3/2
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


