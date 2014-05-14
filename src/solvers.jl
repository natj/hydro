#Rieman problem solvers


#HLLE solver
function hlle(hyd::data1d)
    fluxdiff = zeros(hyd.n, 3)

    #compute eigenvalues
    evl = zeros(hyd.n, 3)
    evr = zeros(hyd.n, 3)
    smin = zeros(hyd.n)
    smax = zeros(hyd.n)
    csp = sqrt(eos_cs2(hyd.rhop, hyd.epsp, gamma))
    csm = sqrt(eos_cs2(hyd.rhom, hyd.epsm, gamma))

    for i = 2:(hyd.n-1)
        evl[i,1] = hyd.velp[i]
        evr[i,1] = hyd.velm[i+1]
        evl[i,2] = hyd.velp[i] - csp[i]
        evr[i,2] = hyd.velm[i+1] - csm[i+1]
        evl[i,3] = hyd.velp[i] + csp[i]
        evr[i,3] = hyd.velm[i+1] +csm[i+1]

        #min and max eigenvalues
        smin[i] = min(evl[i,1], evl[i,2], evl[i,3],
                      evr[i,1], evr[i,2], evr[i,3], 0.0)
        smax[i] = max(evl[i,1], evl[i,2], evl[i,3],
                      evr[i,1], evr[i,2], evr[i,3], 0.0)
    end

    #set up flux left L and right R of the interface
    #at i+1/2
    fluxl = zeros(hyd.n, 3)
    fluxr = zeros(hyd.n, 3)

    for i = 2:(hyd.n-1)
        fluxl[i,1] = hyd.qp[i,1] * hyd.velp[i]
        fluxl[i,2] = hyd.qp[i,2] * hyd.velp[i] + hyd.pressp[i]
        fluxl[i,3] = (hyd.qp[i,3] + hyd.pressp[i]) * hyd.velp[i]

        fluxr[i,1] = hyd.qm[i+1,1] * hyd.velm[i+1]
        fluxr[i,2] = hyd.qm[i+1,2] * hyd.velm[i+1] + hyd.pressm[i+1]
        fluxr[i,3] = (hyd.qm[i+1,3] + hyd.pressm[i+1]) * hyd.velm[i+1]
    end

    #solve the Riemann problem for the i+1/2 interface
    ds = smax .- smin
    flux = zeros(hyd.n, 3)

    for j = 1:3
        for i = hyd.g:(hyd.n-hyd.g+1)
            flux[i,j] = (smax[i]*fluxl[i,j] - smin[i]*fluxr[i,j] + smax[i]*smin[i]*(hyd.qm[i+1,j] - hyd.qp[i,j])) / ds[i]
        end
    end

    #flux difference
    for j = 1:3
        for i = (hyd.g+1):(hyd.n-hyd.g+1)
            rm = hyd.xi[i]
            rp = hyd.xi[i+1]
            dxi = 1.0/(rp - rm)
            fluxdiff[i, j] = dxi * (flux[i, j]  - flux[i-1, j])
        end
    end

    return fluxdiff
end



#HLLC solver
function hllc(hyd::data1d)
    
    #compute eigenvalues
    smin = zeros(hyd.n)
    smax = zeros(hyd.n)
    smid = zeros(hyd.n)
    csp = sqrt(eos_cs2(hyd.rhop, hyd.epsp, gamma))
    csm = sqrt(eos_cs2(hyd.rhom, hyd.epsm, gamma))

    for i = 2:(hyd.n-1)
        #min and max eigenvalues

        #minus = right
        #plus = left
        #min = left
        #max = right

        smin[i] = min(hyd.velp[i],   hyd.velp[i] - csp[i],     hyd.velp[i] + csp[i],
                      hyd.velm[i+1], hyd.velm[i+1] - csm[i+1], hyd.velm[i+1] +csm[i+1])
        smax[i] = max(hyd.velp[i],   hyd.velp[i] - csp[i],     hyd.velp[i] + csp[i],
                      hyd.velm[i+1], hyd.velm[i+1] - csm[i+1], hyd.velm[i+1] +csm[i+1])
        smid[i] = (hyd.pressm[i+1] - hyd.pressp[i] + hyd.rhop[i]*hyd.velp[i]*(smin[i] - hyd.velp[i]) - hyd.rhom[i+1]*hyd.velm[i+1]*(smax[i] - hyd.velm[i+1]))/(hyd.rhop[i]*(smin[i]-hyd.velp[i]) - hyd.rhom[i+1]*(smax[i]-hyd.velm[i+1]))
    end

    #set up flux left L and right R of the interface
    #at i+1/2
    fluxl = zeros(hyd.n, 3)
    fluxr = zeros(hyd.n, 3)

    for i = 2:(hyd.n-1)
        fluxl[i,1] = hyd.qp[i,1] * hyd.velp[i]
        fluxl[i,2] = hyd.qp[i,2] * hyd.velp[i] + hyd.pressp[i]
        fluxl[i,3] = (hyd.qp[i,3] + hyd.pressp[i]) * hyd.velp[i]

        fluxr[i, 1] = hyd.qm[i+1, 1] * hyd.velm[i+1]
        fluxr[i, 2] = hyd.qm[i+1, 2] * hyd.velm[i+1] + hyd.pressm[i+1]
        fluxr[i, 3] = (hyd.qm[i+1, 3] + hyd.pressm[i+1]) * hyd.velm[i+1]
    end


    #solve the Riemann problem for the i+1/2 interface
    ds = smax .- smin
    flux = zeros(hyd.n, 3)

    for i = (hyd.g-1):(hyd.n-hyd.g+1)
    #for i = 1:hyd.n
        if smin[i] >= 0.0
            flux[i,1] = fluxl[i,1]
            flux[i,2] = fluxl[i,2]
            flux[i,3] = fluxl[i,3]
        elseif smin[i] <= 0.0 && smid[i] >= 0.0
            #f_*,left

            #original hllc (p_*,l = p_*,r) (not working)
            #qsl0 = hyd.rhop[i]*((smin[i]-hyd.velp[i])/(smin[i]-smid[i]))
            #qsl3 = hyd.epsp[i]/hyd.rhop[i]+(smid[i]-hyd.velp[i])*(smid[i]+hyd.pressp[i]/(hyd.rhop[i]*(smin[i]-hyd.velp[i])))
            #flux[i,1] = fluxl[i,1] + smin[i]*(qsl0 - hyd.qp[i, 1])
            #flux[i,2] = fluxl[i,2] + smin[i]*(qsl0*smid[i] - hyd.qp[i, 2])
            #flux[i,3] = fluxl[i,3] + smin[i]*(qsl0*qsl3 - hyd.qp[i, 3])

            #hllc v1.1
            #flux[i,1] = smid[i]*(smin[i]*hyd.qp[i,1]-fluxl[i,1])/(smin[i]-smid[i])
            #flux[i,2] = (smid[i]*(smin[i]*hyd.qp[i,2]-fluxl[i,2])+smin[i]*(hyd.pressp[i]+hyd.rhop[i]*(smin[i]-hyd.velp[i])*(smid[i]-hyd.velp[i])))/(smin[i]-smid[i])
            #flux[i,3] = (smid[i]*(smin[i]*hyd.qp[i,3]-fluxl[i,3])+smid[i]*smin[i]*(hyd.pressp[i]+hyd.rhop[i]*(smin[i]-hyd.velp[i])*(smid[i]-hyd.velp[i])))/(smin[i]-smid[i])

            #hllc v1.2 (arithmetic mean of star-area pressures)
            plr = 0.5*(hyd.pressp[i]+hyd.pressm[i+1]+hyd.rhop[i]*(smin[i]-hyd.velp[i])*(smid[i]-hyd.velp[i])+hyd.rhom[i+1]*(smax[i]-hyd.velm[i+1])*(smid[i]-hyd.velm[i+1]))
            flux[i,1] = smid[i]*(smin[i]*hyd.qp[i,1]-fluxl[i,1])/(smin[i]-smid[i])
            flux[i,2] = (smid[i]*(smin[i]*hyd.qp[i,2]-fluxl[i,2]) + smin[i]*plr)/(smin[i]-smid[i])
            flux[i,3] = (smid[i]*(smin[i]*hyd.qp[i,3]-fluxl[i,3]) + smid[i]*smin[i]*plr)/(smin[i]-smid[i])


        elseif smid[i] <= 0.0 && smax[i] > 0.0
            #f_*,right

            #original hllc (p_*,l = p_*,r) (not working)
            #qsl0 = hyd.rhom[i+1]*((smax[i]-hyd.velm[i+1])/(smax[i]-smid[i]))
            #qsl3 = (hyd.epsm[i+1]/hyd.rhom[i+1])+(smid[i]-hyd.velm[i+1])*(smid[i]+hyd.pressm[i+1]/(hyd.rhom[i+1]*(smax[i]-hyd.velm[i+1])))
            #flux[i,1] = fluxr[i,1] + smax[i]*(qsl0 - hyd.qm[i+1, 1])
            #flux[i,2] = fluxr[i,2] + smax[i]*(qsl0*smid[i] - hyd.qm[i+1, 2])
            #flux[i,3] = fluxr[i,3] + smax[i]*(qsl0*qsl3 - hyd.qm[i+1, 3])

            #hllc v1.1
            #flux[i,1] = smid[i]*(smax[i]*hyd.qm[i+1,1]-fluxr[i,1])/(smax[i]-smid[i])
            #flux[i,2] = (smid[i]*(smax[i]*hyd.qm[i+1,2]-fluxr[i,2])+smax[i]*(hyd.pressm[i+1]+hyd.rhom[i+1]*(smax[i]-hyd.velm[i+1])*(smid[i]-hyd.velm[i+1])))/(smax[i]-smid[i])
            #flux[i,3] = (smid[i]*(smax[i]*hyd.qm[i+1,3]-fluxr[i,3])+smid[i]*smax[i]*(hyd.pressm[i+1]+hyd.rhom[i+1]*(smax[i]-hyd.velm[i+1])*(smid[i]-hyd.velm[i+1])))/(smax[i]-smid[i])

            #hllc v1.2 (arithmetic mean of star-area pressures)
            plr = 0.5*(hyd.pressp[i]+hyd.pressm[i+1]+hyd.rhop[i]*(smin[i]-hyd.velp[i])*(smid[i]-hyd.velp[i])+hyd.rhom[i+1]*(smax[i]-hyd.velm[i+1])*(smid[i]-hyd.velm[i+1]))
            flux[i,1] = smid[i]*(smax[i]*hyd.qm[i+1,1]-fluxr[i,1])/(smax[i]-smid[i])
            flux[i,2] = (smid[i]*(smax[i]*hyd.qm[i+1,2]-fluxr[i,2]) + smax[i]*plr)/(smax[i]-smid[i])
            flux[i,3] = (smid[i]*(smax[i]*hyd.qm[i+1,3]-fluxr[i,3]) + smid[i]*smax[i]*plr)/(smax[i]-smid[i])


        elseif smax[i] <= 0.0
            flux[i,1] = fluxr[i,1]
            flux[i,2] = fluxr[i,2]
            flux[i,3] = fluxr[i,3]
        end
    end

    #flux difference
    #fluxdiff = zeros(hyd.n, 3)
    #for j = 1:3
    #    for i = (hyd.g+1):(hyd.n-hyd.g+1)
    #        rm = hyd.xi[i]
    #        rp = hyd.xi[i+1]
    #        dxi = 1.0/(rp - rm)
    #        fluxdiff[i, j] = dxi * (flux[i, j]  - flux[i-1, j])
    #    end
    #end

    return flux
end


#x-sweep
function xsweep(hyd::data2d)
    flux = zeros(hyd.ny, hyd.nx, 4)
    for j = (hyd.g-1):(hyd.ny-hyd.g+1)
        hyd1 = splice(hyd, j, 1)
        fluxi = hllc(hyd1)

        for i = (hyd.g-1):(hyd.nx-hyd.g+1)
            flux[j, i, 1] = fluxi[i, 1]
            flux[j, i, 2] = fluxi[i, 2]
            flux[j, i, 4] = fluxi[i, 3]
        end
    end

    return flux
end

#y-sweep
function ysweep(hyd::data2d)
    flux = zeros(hyd.ny, hyd.nx, 4)
    for j = (hyd.g-1):(hyd.nx-hyd.g+1)
        hyd1 = splice(hyd, j, 2)
        fluxi = hllc(hyd1)

        for i = (hyd.g-1):(hyd.ny-hyd.g+1)
            flux[i, j, 1] = fluxi[i, 1]
            flux[i, j, 3] = fluxi[i, 2]
            flux[i, j, 4] = fluxi[i, 3]
        end
    end

    return flux
end


#Split 2-dimensional flux
function ssflux(hyd::data2d, dt, iter)
    fluxdiff = mod(iter, 2) == 0 ? xsweep(hyd) : ysweep(hyd)
    return fluxdiff
end

function ssyflux(hyd::data2d, dt, iter)
    fluxdiff = ysweep(hyd)
    return fluxdiff
end


#Split 2-dimensional flux
function sflux(hyd::data2d, dt, iter)

    hyd2 = hyd

    qold = copy(hyd2.q)
    fluxdiff = mod(iter, 2) == 0 ? xsweep(hyd2) : ysweep(hyd2)
    hyd2.q = qold - dt*fluxdiff
    hyd2.rho, hyd2.eps, hyd2.press, hyd2.velx, hyd2.vely = con2prim(hyd2.q)
    hyd2 = apply_bcs(hyd2)
    hyd2 = reconstruct(hyd2)
    fluxdiff = mod(iter, 2) == 0 ? ysweep(hyd2) : xsweep(hyd2)

    return fluxdiff
end


function uflux(hyd::data2d, dt)

    #transverse fluxes
    fx = xsweep(hyd)
    fy = ysweep(hyd)

    # U_xl[i,j,:] = U_xl[i,j,:] - 0.5*dt/dy * (F_y[i-1,j+1,:] - F_y[i-1,j,:])
    # U_xr[i,j,:] = U_xr[i,j,:] - 0.5*dt/dy * (F_y[i,j+1,:] - F_y[i,j,:])

    # U_yl[i,j,:] = U_yl[i,j,:] - 0.5*dt/dx * (F_x[i+1,j-1,:] - F_x[i,j-1,:])
    # U_yr[i,j,:] = U_yr[i,j,:] - 0.5*dt/dx * (F_x[i+1,j,:] - F_x[i,j,:])

    #TODO: relax dx & dy calculations
    for j = (hyd.g-1):(hyd.ny-hyd.g+1), i = (hyd.g-1):(hyd.nx-hyd.g+1)
        for k = 1:4
            hyd.qp[j,i,k,1] -= 0.5dt*(fy[j+1, i-1, k] - fy[j, i-1, k])/(hyd.yi[j+1]-hyd.yi[j])
            hyd.qm[j,i,k,1] -= 0.5dt*(fy[j+1, i, k] - fy[j, i, k])/(hyd.yi[j+1]-hyd.yi[j])

            hyd.qp[j,i,k,2] -= 0.5dt*(fx[j-1, i+1, k] - fx[j-1, i, k])/(hyd.xi[i+1]-hyd.xi[i])
            hyd.qm[j,i,k,2] -= 0.5dt*(fx[j, i+1, k] - fx[j, i, k])/(hyd.xi[i+1]-hyd.xi[i])
        end
    end

    #update conserved variables
    hyd.rhom[:,:,1], hyd.epsm[:,:,1], hyd.pressm[:,:,1], hyd.velxm[:,:,1], hyd.velym[:,:,1] = con2prim(hyd.qm[:,:,:,1])
    hyd.rhop[:,:,1], hyd.epsp[:,:,1], hyd.pressp[:,:,1], hyd.velxp[:,:,1], hyd.velyp[:,:,1] = con2prim(hyd.qp[:,:,:,1])

    hyd.rhom[:,:,2], hyd.epsm[:,:,2], hyd.pressm[:,:,2], hyd.velxm[:,:,2], hyd.velym[:,:,2] = con2prim(hyd.qm[:,:,:,2])
    hyd.rhop[:,:,2], hyd.epsp[:,:,2], hyd.pressp[:,:,2], hyd.velxp[:,:,2], hyd.velyp[:,:,2] = con2prim(hyd.qp[:,:,:,2])
  
    #fluxes normal to the surface
    fx = xsweep(hyd)
    fy = ysweep(hyd)

    return fx, fy
end










