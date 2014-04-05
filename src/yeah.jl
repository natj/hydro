#YEt Another Hydro code

#basic parameters
gamma = 1.4
cfl = 0.5
dt = 1.0e-5
dtp = dt

nzones = 200
tend = 0.2

include("grid.jl")


function prim2con(rho, vel, eps)
    q = zeros(3, size(rho, 1))
    q[1,:] = rho
    q[2,:] = rho .* vel
    q[3,:] = rho .* eps .+ 0.5rho .* vel.^2.0

    return q
end

function con2prim(q)
    rho = vec(q[1,:])
    vel = vec(q[2,:]' ./ rho)
    eps = vec(q[3,:]' ./ rho - 0.5vel.^2.0)
    press = eos_press(rho, eps, gamma)

    return rho, eps, press, vel
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

function minmod(a,b)
    if a*b < 0.0
        return 0.0
    elseif abs(a) < abs(b)
        return a
    else
        return b
    end

    nothing
end

signum(x,y) = y >= 0.0 ? abs(x) : -abs(x)
#function signum(x,y)
#    if y >= 0.0
#        return abs(x)
#    else
#        return -abs(x)
#    end
#end


function tvd_mc_reconstruct(n, g, f, x, xi)
    fp = zeros(n)
    fm = zeros(n)

    for i = g:(n-g+2)
        dx_up = x[i] - x[i-1]
        dx_down = x[i+1] - x[i]
        dx_m = x[i] -xi[i]
        dx_p = xi[i+1] - x[i]
        df_up = (f[i]-f[i-1]) / dx_up
        df_down = (f[i+1]-f[i]) / dx_down

        if df_up*df_down < 0.0
            delta = 0.0
        else
            delta = signum(min(2.0abs(df_up), 2.0abs(df_down), 0.5(abs(df_up)+abs(df_down))), df_up + df_down)
        end

        fp[i] = f[i] + delta*dx_p
        fm[i] = f[i] - delta*dx_m
    end

    return fp, fm
end


function reconstruct(hyd)

    hyd.rhop, hyd.rhom = tvd_mc_reconstruct(hyd.n,
                                            hyd.g,
                                            hyd.rho,
                                            hyd.x,
                                            hyd.xi)
    hyd.epsp, hyd.epsm = tvd_mc_reconstruct(hyd.n,
                                            hyd.g,
                                            hyd.eps,
                                            hyd.x,
                                            hyd.xi)
    hyd.velp, hyd.velm = tvd_mc_reconstruct(hyd.n,
                                            hyd.g,
                                            hyd.vel,
                                            hyd.x,
                                            hyd.xi)

    hyd.pressp = eos_press(hyd.rhop, hyd.epsp, gamma)
    hyd.pressm = eos_press(hyd.rhom, hyd.epsm, gamma)

    hyd.qp = prim2con(hyd.rhop, hyd.velp, hyd.epsp)
    hyd.qm = prim2con(hyd.rhom, hyd.velm, hyd.epsm)

    return hyd
end

#equation of state
function eos_press(rho, eps, gamma)
    press = (gamma- 1.0) .* rho .* eps
    return press
end

function eos_cs2(rho, eps, gamma)
    prs = (gamma - 1.0) .* rho .* eps
    dpde = (gamma - 1.0) .* rho
    dpdrho = (gamma - 1.0) .* eps
    cs2 = dpdrho .+ dpde .* prs ./ (rho + 1.0e-30).^2.0
    return cs2
end

#time step calculation
function calc_dt(hyd, dtp)
    cs = sqrt(eos_cs2(hyd.rho, hyd.eps, gamma))
    dtnew = 1.0
    for i = (hyd.g+1):(hyd.n-hyd.g+1)
        dtnew = min(dtnew, (hyd.x[i+1] - hyd.x[i]) / max(abs(hyd.vel[i]+cs[i]), abs(hyd.vel[i]-cs[i])))
    end

    dtnew = min(cfl*dtnew, 1.05*dtp)

    return dtnew
end

#HLLE solver
function hlle(hyd)
    fluxdiff = zeros(3, hyd.n)

    #compute eigenvalues
    evl = zeros(3, hyd.n)
    evr = zeros(3, hyd.n)
    smin = zeros(hyd.n)
    smax = zeros(hyd.n)
    csp = sqrt(eos_cs2(hyd.rhop, hyd.epsp, gamma))
    csm = sqrt(eos_cs2(hyd.rhom, hyd.epsm, gamma))

    for i = 2:(hyd.n-1)
        evl[1,i] = hyd.velp[i]
        evr[1,i] = hyd.velm[i+1]
        evl[2,i] = hyd.velp[i] - csp[i]
        evr[2,i] = hyd.velm[i+1] - csm[i+1]
        evl[3,i] = hyd.velp[i] + csp[i]
        evr[3,i] = hyd.velm[i+1] +csm[i+1]

        #min and max eigenvalues
        smin[i] = min(evl[1,i], evl[2,i], evl[3,i],
                      evr[1,i], evr[2,i], evr[3,i], 0.0)
        smax[i] = max(evl[1,i], evl[2,i], evl[3,i],
                      evr[1,i], evr[2,i], evr[3,i], 0.0)
    end

    #set up flux left L and right R of the interface
    #at i+1/2
    fluxl = zeros(3, hyd.n)
    fluxr = zeros(3, hyd.n)

    for i = 2:(hyd.n-1)
        fluxl[1,i] = hyd.qp[1,i] * hyd.velp[i]
        fluxl[2,i] = hyd.qp[2,i] * hyd.velp[i] + hyd.pressp[i]
        fluxl[3,i] = (hyd.qp[3,i] + hyd.pressp[i]) * hyd.velp[i]

        fluxr[1,i] = hyd.qm[1,i+1] * hyd.velm[i+1]
        fluxr[2,i] = hyd.qm[2,i+1] * hyd.velm[i+1] + hyd.pressm[i+1]
        fluxr[3,i] = (hyd.qm[3,i+1] + hyd.pressm[i+1]) * hyd.velm[i+1]
    end

    #solve the Riemann problem for the i+1/2 interface
    ds = smax .- smin
    flux = zeros(3, hyd.n)
    for i = hyd.g:(hyd.n-hyd.g+1)
        flux[:,i] = (smax[i]*fluxl[:,i] .- smin[i]*fluxr[:,i] .+ smax[i]*smin[i]*(hyd.qm[:,i+1] - hyd.qp[:,i])) / ds[i]
    end

    #flux difference
    for i = (hyd.g+1):(hyd.n-hyd.g+1)
        rm = hyd.xi[i]
        rp = hyd.xi[i+1]
        dxi = 1.0/(rp - rm)
        fluxdiff[:,i] = dxi * (flux[:,i]  .- flux[:,i-1])
    end

    return fluxdiff
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
hyd = data(nzones)

#set up grid
hyd = grid_setup(hyd, 0.0, 1.0)

#set up initial data
hyd = setup_ID(hyd)

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


















