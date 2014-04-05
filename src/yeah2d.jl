#YEt Another Hydro code
#2-dim

#basic parameters
gamma = 1.4
cfl = 0.5
dt = 1.0e-5
dtp = dt

nx = 200
ny = 300
tend = 0.2

include("grid2d.jl")
include("eos.jl")
include("reconstr.jl")
include("solvers.jl")


#2-dim
function prim2con(rho::AbstractMatrix,
                  vel::AbstractMatrix,
                  eps::AbstractMatrix)

    q = zeros(3, size(rho, 1), size(rho, 2))
    q[1,:,:] = rho
    q[2,:,:] = rho .* vel
    q[3,:,:] = rho .* eps .+ 0.5rho .* vel.^2.0

    return q
end


#2-dim
function con2prim(q)
    rho = vec(q[1,:,:])
    vel = vec(q[2,:,:]' ./ rho)
    eps = vec(q[3,:,:]' ./ rho - 0.5vel.^2.0)
    press = eos_press(rho, eps, gamma)

    return rho, eps, press, vel
end

#time step calculation
function calc_dt(hyd, dtp)
    cs = sqrt(eos_cs2(hyd.rho, hyd.eps, gamma))
    dtnew = 1.0
    for j = (hyd.g+1):(hyd.ny-hyd.g+1), i = (hyd.g+1):(hyd.nx-hyd.g+1)
        dtnew = min(dtnew, min((hyd.y[i+1] - hyd.y[i]), hyd.x[i+1] - hyd.x[i+1]) / max(abs(hyd.vel[j, i]+cs[j, i]), abs(hyd.vel[j, i]-cs[j, i])))
    end
    dtnew = min(cfl*dtnew, 1.05*dtp)

    return dtnew
end


function calc_rhs(hyd)
    #reconstruction and prim2con
    hyd = reconstruct(hyd)
    #compute flux difference
    fluxdiff = hlle(hyd)
    #return RHS = -fluxdiff
    return -fluxdiff
end



###############
# main program
###############

#initialize
hyd = data2d(nx, ny)

#set up grid
hyd = grid_setup(hyd, 0.0, 1.0, 0.0, 1.5)

#set up initial data
hyd = setup_blast(hyd)

#get initial timestep
dt = calc_dt(hyd, dt)

#initial prim2con
hyd.q = prim2con(hyd.rho, hyd.vel, hyd.eps)

t = 0.0
i = 1

#display
#plot(hyd.x, hyd.rho, "r-")

while t < tend

    if i % 10 == 0
        sleep(1.0)
        #output
        println("$i $t $dt")
        p=plot(hyd.x, hyd.rho, "r-")
        p=oplot(hyd.x, hyd.vel, "b-")
        p=oplot(hyd.x, hyd.press, "g-")
        display(p)
    end

    #calculate new timestep
    dt = calc_dt(hyd, dt)

    #save old state
    hydold = hyd
    qold = hyd.q

    #calc rhs
    k1 = calc_rhs(hyd)
    #calculate intermediate step
    hyd.q = qold + 0.5dt*k1
    #con2prim
    hyd.rho, hyd.eps, hyd.press, hyd.vel = con2prim(hyd.q)
    #boundaries
    hyd = apply_bcs(hyd)

    #calc rhs
    k2 = calc_rhs(hyd)
    #apply update
    hyd.q = qold + dt*(0.5k1 + 0.5k2)
    #con2prim
    hyd.rho, hyd.eps, hyd.press, hyd.vel = con2prim(hyd.q)
    #apply bcs
    hyd = apply_bcs(hyd)

    #update time
    t += dt
    i += 1

end


















