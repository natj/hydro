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
    for i = hyd.g:(hyd.n-hyd.g+1)
        flux[i,:] = (smax[i]*fluxl[i,:] .- smin[i]*fluxr[i,:] .+ smax[i]*smin[i]*(hyd.qm[i+1,:] - hyd.qp[i,:])) / ds[i]
    end

    #flux difference
    for i = (hyd.g+1):(hyd.n-hyd.g+1)
        rm = hyd.xi[i]
        rp = hyd.xi[i+1]
        dxi = 1.0/(rp - rm)
        fluxdiff[i,:] = dxi * (flux[i,:]  .- flux[i-1,:])
    end

    return fluxdiff
end

#Split 2-dimensional flux
#Strang splitting
function sflux(hyd::data2d, iter)

    #x-sweep
    function xsweep(hyd::data2d)
        fluxdiff = zeros(hyd.ny, hyd.nx, 4)
        for j = (hyd.g+1):(hyd.ny-hyd.g+1)
            #println("x  j=$j")
            hyd1 = splice(hyd, j, 1)
            fluxdiffi = hlle(hyd1)

            fluxdiff[j, :, 1] = fluxdiffi[:, 1]
            fluxdiff[j, :, 2] = fluxdiffi[:, 2]
            fluxdiff[j, :, 4] = fluxdiffi[:, 3]
        end
        return fluxdiff
    end

    #y-sweep
    function ysweep(hyd::data2d)
        fluxdiff = zeros(hyd.ny, hyd.nx, 4)
        for j = (hyd.g+1):(hyd.nx-hyd.g+1)
            #println("y  j=$j")
            hyd1 = splice(hyd, j, 2)
            fluxdiffi = hlle(hyd1)

            fluxdiff[:, j, 1] = fluxdiffi[:, 1]'
            fluxdiff[:, j, 3] = fluxdiffi[:, 2]'
            fluxdiff[:, j, 4] = fluxdiffi[:, 3]'
        end
        return fluxdiff
    end

    qold = hyd.q
    fluxdiff = mod(iter, 2) == 0 ? xsweep(hyd) : ysweep(hyd)
    hyd.q = qold - dt*fluxdiff
    hyd.rho, hyd.eps, hyd.press, hyd.velx, hyd.vely = con2prim(hyd.q)
    hyd = apply_bcs(hyd)
    hyd = reconstruct(hyd)
    fluxdiff = mod(iter, 2) == 0 ? ysweep(hyd) : xsweep(hyd)

    return fluxdiff
end













