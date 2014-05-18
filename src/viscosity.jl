#Artificial viscosity


#Artificial viscosity according to Colella & Woodward 1984 (eq. 4.5)
function avisc_CW(u, v, g, nx, ny, dx, dy)

    const cvisc = 0.1

    avisco_x = zeros(ny, nx)
    avisco_y = zeros(ny, nx)

    for j = (g-1):(ny-g+1), i = (g-1):(nx-g+1)

        #x-interface
        divU_x = (u[j, i] - u[j, i-1])/dx + 0.25(v[j+1, i] + v[j+1, i-1] - v[j-1, i] - v[j-1, i-1])/dx
        avisco_x[j, i] = cvisc*max(-divU_x*dx, 0.0)

        #y-interface value
        divU_y = 0.25(u[j, i+1] + u[j-1, i+1] - u[j, i-1] - u[j-1, i-1])/dx + (v[j, i] - v[j-1, i])/dy
        avisco_y[j, i] = cvisc*max(-divU_y*dy, 0.0)
    end

    return avisco_x, avisco_y
end

#General function for artificial viscosity 
function artificial_viscosity(hyd::data2d, dt)

    dx = abs(hyd.x[2]-hyd.x[1])
    dy = abs(hyd.y[2]-hyd.y[1])
    aviscx, aviscy = avisc_CW(hyd.velx, hyd.vely, hyd.g, hyd.nx, hyd.ny, dx, dy)

    fx = zeros(hyd.ny, hyd.nx, 4)
    fy = zeros(hyd.ny, hyd.nx, 4)

    for j = (hyd.g):(hyd.ny-hyd.g+1), i = (hyd.g):(hyd.nx-hyd.g+1)
        for k = 1:4
            fx[j, i, k] = aviscx[j, i]*(hyd.q[j, i-1, k] - hyd.q[j, i, k])
            fy[j, i, k] = aviscy[j, i]*(hyd.q[j-1, i, k] - hyd.q[j, i, k])
        end
    end

    return fx, fy
end
