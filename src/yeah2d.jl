#YEt Another Hydro code
#2-dim


include("grid.jl")
include("visualize.jl")
include("eos.jl")
include("reconstr.jl")
include("solvers.jl")
include("gravity.jl")

using Winston

#2-dim
function prim2con(rho::AbstractMatrix,
                  velx::AbstractMatrix,
                  vely::AbstractMatrix,
                  eps::AbstractMatrix)

    q = zeros(size(rho, 1), size(rho, 2), 4)
    q[:, :, 1] = rho
    q[:, :, 2] = rho .* velx
    q[:, :, 3] = rho .* vely
    q[:, :, 4] = rho .* eps .+ 0.5rho .* (velx.^2.0 .+ vely.^2.0)
    #q[:, :, 4] = clamp(q[:,:,4], 1.0e-10, 1.0e10)

    return q
end


#2-dim
function con2prim(q)
    rho = q[:, :, 1]
    velx = q[:, :, 2] ./ rho
    vely = q[:, :, 3] ./ rho
    eps = q[:, :, 4] ./ rho .- 0.5(velx.^2.0 .+ vely.^2.0)
    #eps = clamp(eps, 1.0e-10, 1.0e10)

    #for i = 1:50, j = 1:100
    #    if eps[j, i] < 0.0
    #        println("eps = $(eps[j, i]) i=$i j=$j")
    #        println("$(q[j, i, 4]) $(rho[j,i]) $(0.5(velx[j,i]^2.0 + vely[j,i]^2.0)) ")
    #    end
    #end


    press = eos_press(rho, eps, gamma)

    return rho, eps, press, velx, vely
end

#time step calculation
function calc_dt(hyd, dtp)
    cs = sqrt(eos_cs2(hyd.rho, hyd.eps, gamma))
    dtnew = 1.0
    for j = (hyd.g+1):(hyd.ny-hyd.g+1), i = (hyd.g+1):(hyd.nx-hyd.g+1)
        dtnew = min(dtnew, (hyd.x[i+1] - hyd.x[i]) / max(abs(hyd.velx[j, i]+cs[j, i]), abs(hyd.velx[j, i]-cs[j, i])))
        dtnew = min(dtnew, (hyd.y[j+1] - hyd.y[j]) / max(abs(hyd.vely[j, i]+cs[j, i]), abs(hyd.vely[j, i]-cs[j, i])))
    end

    dtnew = min(cfl*dtnew, 1.05*dtp)

    return dtnew
end

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

#Colella & Woodward 1984 (eq. 4.5)
function avisc_CW(u, v, g, nx, ny, dx, dy)

    const cvisc = 0.1

    avisco_x = zeros(ny, nx)
    avisco_y = zeros(ny, nx)

    for j = g:(ny-g+2), i = g:(nx-g+2)

        #x-interface
        divU_x = (u[j, i] - u[j, i-1])/dx + 0.25(v[j+1, i] + v[j+1, i-1] - v[j-1, i] - v[j-1, i-1])/dy
        avisco_x[j, i] = cvisc*max(-divU_x*dx, 0.0)

        #y-interface value
        divU_y = 0.25(u[j, i+1] + u[j-1, i+1] - u[j, i-1] - u[j-1, i-1])/dx + (v[j, i] - v[j-1, i])/dy
        avisco_y[j, i] = cvisc*max(-divU_y*dy, 0.0)
    end

    return avisco_x, avisco_y
end

function artificial_viscosity(hyd::data2d, dt)

    dx = abs(hyd.x[2]-hyd.x[1])
    dy = abs(hyd.y[2]-hyd.y[1])
    aviscx, aviscy = avisc_CW(hyd.velx, hyd.vely, hyd.g, hyd.nx, hyd.ny, dx, dy)

    fx = zeros(hyd.ny, hyd.nx, 4)
    fy = zeros(hyd.ny, hyd.nx, 4)

    for k = 1:4
        for j = hyd.g:(hyd.ny-hyd.g+2), i = hyd.g:(hyd.nx-hyd.g+2)
            fx[j, i, k] = aviscx[j, i]*(hyd.q[j, i-1, k] - hyd.q[j, i, k])
            fy[j, i, k] = aviscy[j, i]*(hyd.q[j-1, i, k] - hyd.q[j, i, k])
        end
    end

    return fx, fy
end


function calc_rhs(hyd, dt, iter)

    #reconstruction and prim2con
    hyd = reconstruct(hyd)

    #apply source terms for the interfaces
    #hyd = source_terms(hyd, dt)

    #compute flux difference
    #fluxdiff = sflux(hyd, dt, iter)
    fx, fy = uflux(hyd, dt)

    #add artificial viscosity
    vfx, vfy = artificial_viscosity(hyd, dt)
    fx += vfx
    fy += vfy

    #hyd = source_terms(hyd, 2.0dt)

    dx = abs(hyd.x[2]-hyd.x[1])
    dy = abs(hyd.y[2]-hyd.y[1])

    fluxdiff = zeros(hyd.ny, hyd.nx, 4)
    for j = (hyd.g+1):(hyd.ny-hyd.g+1), i = (hyd.g+1):(hyd.nx-hyd.g+1)
        for k = 1:4
            fluxdiff[j, i, k] = (fx[j, i, k] - fx[j, i-1, k])/dx + (fy[j, i, k] - fy[j-1, i, k])/dy
        end
    end


    #return RHS = -fluxdiff
    return -fluxdiff
end



###############
# main program
###############

function evolve(hyd, tend, gamma, cfl, nx, ny)

    dt = 1.0e-5
    dtp = dt

    #get initial timestep
    dt = calc_dt(hyd, dt)

    #initial prim2con
    hyd.q = prim2con(hyd.rho, hyd.velx, hyd.vely, hyd.eps)

    t = 0.0
    i = 1

    visualize(hyd)

    while t < tend

        if i % 10 == 0
            println("$i $t $dt")
            visualize(hyd)
        end

        #calculate new timestep
        dt = calc_dt(hyd, dt)

        #save old state
        #hydold = hyd
        qold = copy(hyd.q)

        #calc rhs
        k1 = calc_rhs(hyd, 0.5dt, i)
        #calculate intermediate step
        hyd.q = qold + 0.5dt*k1
        #con2prim
        hyd.rho, hyd.eps, hyd.press, hyd.velx, hyd.vely = con2prim(hyd.q)

        #gxflx, gyflx = ygravity(0.5(hyd.rho[:,:]+hydold.rho[:,:]), hyd.nx, hyd.ny)
        #hyd.velx[:,:] -= 0.5dt*gxflx
        #hyd.vely[:,:] -= 0.5dt*gyflx
        #hyd.eps[:,:] -= 0.5dt*((hyd.velx[:,:] .* gxflx) .+ (hyd.vely[:,:] .* gyflx))

        #boundaries
        hyd = apply_bcs(hyd)

        #calc rhs
        k2 = calc_rhs(hyd, dt, i)
        #apply update
        hyd.q = qold + dt*(0.5k1 + 0.5k2)
        #con2prim
        hyd.rho, hyd.eps, hyd.press, hyd.velx, hyd.vely = con2prim(hyd.q)
        #apply bcs
        hyd = apply_bcs(hyd)

        #update time
        t += dt
        i += 1

        #if i > 100; break; end
    end

    return hyd
end


#basic parameters
gamma = 1.4
cfl = 0.5

nx = 100
ny = 100
tend = 5.0

#initialize
hyd = data2d(nx, ny)
#r32x, r32y = r32kernel(hyd.x, hyd.y)

#set up grid
hyd = grid_setup(hyd, 0.0, 1.0, 0.0, 1.0)
#hyd = grid_setup(hyd, 0.0, 1.0, -0.5, 0.5)

#set up initial data
#hyd = setup_taylor2(hyd)
hyd = setup_blast(hyd)
#hyd = setup_tubexy(hyd)
#hyd = setup_tubey(hyd)
#hyd = setup_collision(hyd)
#hyd = setup_fall(hyd)
#hyd = setup_kh(hyd)

visualize(hyd)

#main integration loop
hyd = evolve(hyd, tend, gamma, cfl, nx, ny)
#hyd = @profile evolve(hyd, tend, gamma, cfl, nx, ny)














