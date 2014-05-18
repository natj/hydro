#Time integrators


#Second order Runge-Kutta
function RK2(hyd::data2d, dt)

    #save old state
    old_rho = copy(hyd.rho[:,:])
    qold = copy(hyd.q)

    #calc rhs
    k1 = calc_rhs(hyd, 0.5dt)
    #calculate intermediate step
    hyd.q = qold + 0.5dt*k1
    #con2prim
    hyd.rho, hyd.eps, hyd.press, hyd.velx, hyd.vely = con2prim(hyd.q)

    #gxflx, gyflx = ygravity(0.5(hyd.rho[:,:]+old_rho[:,:]), hyd.nx, hyd.ny)
    #hyd.velx[:,:] -= 0.5dt*gxflx
    #hyd.vely[:,:] -= 0.5dt*gyflx
    #hyd.eps[:,:] -= 0.5dt*((hyd.velx[:,:] .* gxflx) .+ (hyd.vely[:,:] .* gyflx))

    #boundaries
    hyd = apply_refl_bcs(hyd)

    #calc rhs
    k2 = calc_rhs(hyd, dt)
    #apply update
    hyd.q = qold + dt*(0.5k1 + 0.5k2)
    #con2prim
    hyd.rho, hyd.eps, hyd.press, hyd.velx, hyd.vely = con2prim(hyd.q)

    #gxflx, gyflx = ygravity(0.5(hyd.rho[:,:]+old_rho[:,:]), hyd.nx, hyd.ny)
    #hyd.velx[:,:] -= dt*gxflx
    #hyd.vely[:,:] -= dt*gyflx
    #hyd.eps[:,:] -= dt*((hyd.velx[:,:] .* gxflx) .+ (hyd.vely[:,:] .* gyflx))

    #apply bcs
    hyd = apply_refl_bcs(hyd)

    return hyd
end
