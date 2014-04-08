#YEt Another Hydro code
#2-dim


include("grid2d.jl")
include("eos.jl")
include("reconstr.jl")
include("solvers.jl")

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

    #dtnew = max(dtnew, 1.0e-5)
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

    #println(hyd.rho)
    cm = Uint32[Color.convert(Color.RGB24,c) for c in Color.colormap("RdBu")]

    #rho
    hdata = hyd.rho
    p1=FramedPlot()
    #clims = (minimum(hdata), maximum(hdata))
    clims = (0.0, 1.0)
    img = Winston.data2rgb(hdata, clims, cm)
    add(p1, Image((hyd.x[1], hyd.x[end]), (hyd.y[1], hyd.y[end]), img;))

    #pressure
    hdata = hyd.press
    p2=FramedPlot()
    #clims = (minimum(hdata), maximum(hdata))
    clims = (0.0, 1.0)
    img = Winston.data2rgb(hdata, clims, cm)
    add(p2, Image((hyd.x[1], hyd.x[end]), (hyd.y[1], hyd.y[end]), img;))

    #vel
    hdata = sqrt(hyd.velx.^2.0 .+ hyd.vely.^2.0)
    p3=FramedPlot()
    #clims = (minimum(hdata), maximum(hdata))
    clims = (0.0, 1.0)
    img = Winston.data2rgb(hdata, clims, cm)
    add(p3, Image((hyd.x[1], hyd.x[end]), (hyd.y[1], hyd.y[end]), img;))

    #eps
    hdata = hyd.eps
    p4=FramedPlot()
    #clims = (minimum(hdata), maximum(hdata))
    clims = (0.0, 1.0)
    img = Winston.data2rgb(hdata, clims, cm)
    add(p4, Image((hyd.x[1], hyd.x[end]), (hyd.y[1], hyd.y[end]), img;))

    t = Table(1,4)
    t[1,1] = p1
    t[1,2] = p2
    t[1,3] = p3
    t[1,4] = p4
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

    #display
    #plot(hyd.x, hyd.rho, "r-")
    visualize(hyd)


    while t < tend
        if i % 10 == 0
            #sleep(1.0)
            #output
            println("$i $t $dt")
    #        p=plot(hyd.x, hyd.rho, "r-")
    #        p=oplot(hyd.x, hyd.vel, "b-")
    #        p=oplot(hyd.x, hyd.press, "g-")
    #        display(p)
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

        #if i > 10; break; end
    end

    return hyd
end


#basic parameters
gamma = 1.4
cfl = 0.5

nx = 50
ny = 200
tend = 10.0

#initialize
hyd = data2d(nx, ny)

#set up grid
hyd = grid_setup(hyd, 0.0, 0.5, 0.0, 5.0)

#set up initial data
hyd = setup_taylor(hyd)
#hyd = setup_blast(hyd)
#hyd = setup_tubey(hyd)

#main integration loop
hyd = evolve(hyd, tend, gamma, cfl, nx, ny)
#hyd = @profile evolve(hyd, tend, gamma, cfl, nx, ny)














