#2-dim version

include("grid.jl")
include("boundaries.jl")
include("visualize.jl")
include("eos.jl")
include("reconstruct.jl")
include("rsolvers.jl")
include("gravity.jl")
include("viscosity.jl")
include("integrators.jl")

using Winston

#Primary variables to conserved variables
function prim2con(rho::AbstractMatrix,
                  velx::AbstractMatrix,
                  vely::AbstractMatrix,
                  eps::AbstractMatrix)

    q = zeros(size(rho, 1), size(rho, 2), 4)
    q[:, :, 1] = rho
    q[:, :, 2] = rho .* velx
    q[:, :, 3] = rho .* vely
    q[:, :, 4] = rho .* eps .+ 0.5rho .* (velx.^2.0 .+ vely.^2.0)

    return q
end


#Conserved variables to primary variables
function con2prim(q)
    rho = q[:, :, 1]
    velx = q[:, :, 2] ./ rho
    vely = q[:, :, 3] ./ rho
    eps = q[:, :, 4] ./ rho .- 0.5(velx.^2.0 .+ vely.^2.0)

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


#snapshot of the simulation
function snapshot(hyd, it)

            cm = Uint32[Color.convert(Color.RGB24,c) for c in flipud(Color.colormap("RdBu"))]
            xs = hyd.g+1
            xe = hyd.nx-hyd.g
            ye = hyd.ny-hyd.g

            hdata = hyd.rho[xs:ye, xs:xe]
            pf=FramedPlot()
            clims = (0.8, 2.1)
            img = Winston.data2rgb(hdata, clims, cm)
            add(pf, Image((hyd.x[xs], hyd.x[xe]), (hyd.y[xs], hyd.y[ye]), img;))
            setattr(pf, xrange=(hyd.x[xs], hyd.x[xe]))
            setattr(pf, yrange=(hyd.y[xs], hyd.y[ye]))
            #setattr(pf, title="rho")
            Winston.file(pf, "film_$(it).png")

    return nothing
end

#Construct the right hand side of the group of equations
function calc_rhs(hyd, dt)

    #reconstruction and prim2con
    hyd = reconstruct(hyd)

    #apply source terms for the interfaces
    #hyd = source_terms(hyd, dt)

    #unsplit flux
    fx, fy = uflux(hyd, dt)

    #add artificial viscosity
    #vfx, vfy = artificial_viscosity(hyd, dt)
    #fx += vfx
    #fy += vfy

    #flux difference
    fluxdiff = zeros(hyd.ny, hyd.nx, 4)
    for j = (hyd.g-1):(hyd.ny-hyd.g+1), i = (hyd.g-1):(hyd.nx-hyd.g+1)
        dx = abs(hyd.x[i+1]-hyd.x[i])
        dy = abs(hyd.y[j+1]-hyd.y[j])

        for k = 1:4
            fluxdiff[j, i, k] = (fx[j, i, k] - fx[j, i-1, k])/dx + (fy[j, i, k] - fy[j-1, i, k])/dy
        end
    end

    return -fluxdiff
end



###############
# main program
###############

function evolve(hyd, tend, gamma, cfl, nx, ny)

    time = linspace(0.0, tend, int(tend*100))
    it = 1

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

        if time[it] <= t
            #snapshot(hyd, it)
            it += 1
        end

        #calculate new timestep
        dt = calc_dt(hyd, dt)

        #time integration
        hyd = RK2(hyd, dt)

        #update time
        t += dt
        i += 1

    end

    return hyd
end


#basic parameters
gamma = 1.4
cfl = 0.6

nx = 8
ny = 500
tend = 5.0


#initialize
hyd = data2d(nx, ny)
#r32x, r32y = r32kernel(hyd.x, hyd.y)

#set up grid

include("../tests/1dim_hd.jl")
include("../tests/2dim_hd.jl")

hyd = grid_setup(hyd, 0.0, 1.0, 0.0, 1.0)
#hyd = grid_setup(hyd, 0.0, 1.0, 0.0, 0.5)
#hyd = grid_setup(hyd, -0.25, 0.25, -0.75, 0.75)
#hyd = grid_setup(hyd, 0., 1., -0.5, 0.5)

#set up initial data
#hyd = setup_taylor2(hyd)
#hyd = setup_blast(hyd)
#hyd = setup_tubexy(hyd)
hyd = setup_tubey(hyd)
#hyd = setup_collision(hyd)
#hyd = setup_fall(hyd)
#hyd = setup_kh(hyd)

visualize(hyd)

#main integration loop
hyd = evolve(hyd, tend, gamma, cfl, nx, ny)
#hyd = @profile evolve(hyd, tend, gamma, cfl, nx, ny)














