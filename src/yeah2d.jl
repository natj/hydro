#YEt Another Hydro code
#2-dim


include("grid.jl")
include("eos.jl")
include("reconstr.jl")
include("solvers.jl")

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
    q[:, :, 4] = rho .* eps .+ 0.5rho .* (velx.^2.0 + vely.^2.0)

    return q
end


#2-dim
function con2prim(q)
    rho = q[:, :, 1]
    velx = q[:, :, 2] ./ rho
    vely = q[:, :, 3] ./ rho
    eps = q[:, :, 4] ./ rho - 0.5(velx.^2.0 + vely.^2.0)
    eps = clamp(eps, 1.0e-5, 1.0e10)
    press = eos_press(rho, eps, gamma)

    return rho, eps, press, velx, vely
end

#time step calculation
function calc_dt(hyd, dtp)
    cs = sqrt(eos_cs2(hyd.rho, hyd.eps, gamma))
    dtnew = 1.0
    for j = (hyd.g+1):(hyd.ny-hyd.g+1), i = (hyd.g+1):(hyd.nx-hyd.g+1)
        dtnew = min(dtnew, (hyd.x[i+1] - hyd.x[i]) / max(abs(hyd.velx[j, i]+cs[j, i]), abs(hyd.velx[j, i]-cs[j, i])))
    end
    for j = (hyd.g+1):(hyd.ny-hyd.g+1), i = (hyd.g+1):(hyd.nx-hyd.g+1)
        dtnew = min(dtnew, (hyd.y[j+1] - hyd.y[j]) / max(abs(hyd.vely[j, i]+cs[j, i]), abs(hyd.vely[j, i]-cs[j, i])))
    end

    dtnew = min(cfl*dtnew, 1.05*dtp)

    return dtnew
end


function calc_rhs(hyd, dt, iter)
    #reconstruction and prim2con
    hyd = reconstruct(hyd)

    #compute flux difference
    fluxdiff = sflux(hyd, dt, iter)

    #return RHS = -fluxdiff
    return -fluxdiff
end

function visualize(hyd)

    cm = Uint32[Color.convert(Color.RGB24,c) for c in flipud(Color.colormap("RdBu"))]

    hdata = hyd.rho[(hyd.g+1):(hyd.ny-hyd.g), (hyd.g+1):(hyd.nx-hyd.g)]
    p1=FramedPlot()
    #clims = (minimum(hdata), maximum(hdata))
    clims = (0.0, 1.0)
    img = Winston.data2rgb(hdata, clims, cm)
    add(p1, Image((hyd.x[hyd.g+1], hyd.x[hyd.nx-hyd.g]), (hyd.y[hyd.g+1], hyd.y[hyd.ny-hyd.g]), img;))
    setattr(p1, xrange=(hyd.x[hyd.g+1], hyd.x[hyd.nx-hyd.g]))
    setattr(p1, yrange=(hyd.y[hyd.g+1], hyd.y[hyd.ny-hyd.g]))
    setattr(p1, title="rho")

    #pressure
    hdata = hyd.press[(hyd.g+1):(hyd.ny-hyd.g), (hyd.g+1):(hyd.nx-hyd.g)]
    p2=FramedPlot()
    #clims = (minimum(hdata), maximum(hdata))
    clims = (0.0, 1.0)
    img = Winston.data2rgb(hdata, clims, cm)
    add(p2, Image((hyd.x[hyd.g+1], hyd.x[hyd.nx-hyd.g]), (hyd.y[hyd.g+1], hyd.y[hyd.ny-hyd.g]), img;))
    setattr(p2, xrange=(hyd.x[hyd.g+1], hyd.x[hyd.nx-hyd.g]))
    setattr(p2, yrange=(hyd.y[hyd.g+1], hyd.y[hyd.ny-hyd.g]))
    setattr(p2, title="press")

    #vel
    hdata = sqrt(hyd.velx.^2.0 .+ hyd.vely.^2.0)[(hyd.g+1):(hyd.ny-hyd.g), (hyd.g+1):(hyd.nx-hyd.g)]
    p3=FramedPlot()
    #clims = (minimum(hdata), maximum(hdata))
    clims = (0.0, 1.0)
    img = Winston.data2rgb(hdata, clims, cm)
    add(p3, Image((hyd.x[hyd.g+1], hyd.x[hyd.nx-hyd.g]), (hyd.y[hyd.g+1], hyd.y[hyd.ny-hyd.g]), img;))
    setattr(p3, xrange=(hyd.x[hyd.g+1], hyd.x[hyd.nx-hyd.g]))
    setattr(p3, yrange=(hyd.y[hyd.g+1], hyd.y[hyd.ny-hyd.g]))
    setattr(p3, title="vel")

    #eps
    hdata = hyd.eps[(hyd.g+1):(hyd.ny-hyd.g), (hyd.g+1):(hyd.nx-hyd.g)]
    p4=FramedPlot()
    clims = (minimum(hdata), maximum(hdata))
    #clims = (0.0, 1.0)
    img = Winston.data2rgb(hdata, clims, cm)
    add(p4, Image((hyd.x[hyd.g+1], hyd.x[hyd.nx-hyd.g]), (hyd.y[hyd.g+1], hyd.y[hyd.ny-hyd.g]), img;))
    setattr(p4, xrange=(hyd.x[hyd.g+1], hyd.x[hyd.nx-hyd.g]))
    setattr(p4, yrange=(hyd.y[hyd.g+1], hyd.y[hyd.ny-hyd.g]))
    setattr(p4, title="eps")

    t = Table(2,2)
    t[1,1] = p1
    t[1,2] = p2
    t[2,1] = p3
    t[2,2] = p4
    display(t)

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
        hydold = hyd
        qold = hyd.q

        #calc rhs
        k1 = calc_rhs(hyd, 0.5dt, i)
        #calculate intermediate step
        hyd.q = qold + 0.5dt*k1
        #con2prim
        hyd.rho, hyd.eps, hyd.press, hyd.velx, hyd.vely = con2prim(hyd.q)
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

nx = 128
ny = 128
tend = 0.2

#initialize
hyd = data2d(nx, ny)

#set up grid
hyd = grid_setup(hyd, 0.0, 0.5, 0.0, 0.5)

#set up initial data
#hyd = setup_taylor(hyd)
#hyd = setup_blast(hyd)
hyd = setup_tubexy(hyd)

#main integration loop
hyd = evolve(hyd, tend, gamma, cfl, nx, ny)
#hyd = @profile evolve(hyd, tend, gamma, cfl, nx, ny)














