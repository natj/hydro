#Rieman problem solvers


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