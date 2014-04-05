#Reconstruction

function minmod(a,b)
    if a*b < 0.0
        return 0.0
    elseif abs(a) < abs(b)
        return a
    else
        return b
    end
end

signum(x,y) = y >= 0.0 ? abs(x) : -abs(x)

function tvd_mc_reconstruction(n, g, f, x, xi)
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

    hyd.rhop, hyd.rhom = tvd_mc_reconstruction(hyd.n,
                                               hyd.g,
                                               hyd.rho,
                                               hyd.x,
                                               hyd.xi)
    hyd.epsp, hyd.epsm = tvd_mc_reconstruction(hyd.n,
                                               hyd.g,
                                               hyd.eps,
                                               hyd.x,
                                               hyd.xi)
    hyd.velp, hyd.velm = tvd_mc_reconstruction(hyd.n,
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
