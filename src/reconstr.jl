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

#2-dim reconstruction
function reconstruct(hyd::data2d)

  #x-pencils
  for j = hyd.g:(hyd.ny-hyd.g+2)
    rhop, rhom = tvd_mc_reconstruction(hyd.nx,
                                       hyd.g,
                                       vec(hyd.rho[j,:]),
                                       hyd.x,
                                       hyd.xi)

    hyd.rhop[j,:,1] = rhop'
    hyd.rhom[j,:,1] = rhom'

    epsp, epsm = tvd_mc_reconstruction(hyd.nx,
                                       hyd.g,
                                       vec(hyd.eps[j,:]),
                                       hyd.x,
                                       hyd.xi)

    hyd.epsp[j,:,1] = epsp'
    hyd.epsm[j,:,1] = epsm'

    velxp, velxm = tvd_mc_reconstruction(hyd.nx,
                                         hyd.g,
                                         vec(hyd.velx[j,:]),
                                         hyd.x,
                                         hyd.xi)

    hyd.velxp[j,:,1] = velxp'
    hyd.velxm[j,:,1] = velxm'

    velyp, velym = tvd_mc_reconstruction(hyd.nx,
                                         hyd.g,
                                         vec(hyd.vely[j,:]),
                                         hyd.x,
                                         hyd.xi)

    hyd.velyp[j,:,1] = velyp'
    hyd.velym[j,:,1] = velym'
  end

  #y-pencils
  for i = hyd.g:(hyd.nx-hyd.g+2)
    hyd.rhop[:,i,2], hyd.rhom[:,i,2] = tvd_mc_reconstruction(hyd.ny,
                                                             hyd.g,
                                                             hyd.rho[:,i],
                                                             hyd.y,
                                                             hyd.yi)

    hyd.epsp[:,i,2], hyd.epsm[:,i,2] = tvd_mc_reconstruction(hyd.ny,
                                                             hyd.g,
                                                             hyd.eps[:,i,1],
                                                             hyd.y,
                                                             hyd.yi)

    hyd.velxp[:,i,2], hyd.velxm[:,i,2] = tvd_mc_reconstruction(hyd.ny,
                                                               hyd.g,
                                                               hyd.velx[:,i],
                                                               hyd.y,
                                                               hyd.yi)

    hyd.velyp[:,i,2], hyd.velym[:,i,2] = tvd_mc_reconstruction(hyd.ny,
                                                               hyd.g,
                                                               hyd.vely[:,i],
                                                               hyd.y,
                                                               hyd.yi)
  end


  hyd.pressp[:,:,1] = eos_press(hyd.rhop[:,:,1], hyd.epsp[:,:,1], gamma)
  hyd.pressm[:,:,1] = eos_press(hyd.rhom[:,:,1], hyd.epsm[:,:,1], gamma)
  hyd.pressp[:,:,2] = eos_press(hyd.rhop[:,:,2], hyd.epsp[:,:,2], gamma)
  hyd.pressm[:,:,2] = eos_press(hyd.rhom[:,:,2], hyd.epsm[:,:,2], gamma)

  hyd.qp[:,:,:,1] = prim2con(hyd.rhop[:,:,1], hyd.velxp[:,:,1], hyd.velyp[:,:,1], hyd.epsp[:,:,1])
  hyd.qm[:,:,:,1] = prim2con(hyd.rhom[:,:,1], hyd.velxm[:,:,1], hyd.velym[:,:,1], hyd.epsm[:,:,1])
  hyd.qp[:,:,:,2] = prim2con(hyd.rhop[:,:,2], hyd.velxp[:,:,2], hyd.velyp[:,:,2], hyd.epsp[:,:,2])
  hyd.qm[:,:,:,2] = prim2con(hyd.rhom[:,:,2], hyd.velxm[:,:,2], hyd.velym[:,:,2], hyd.epsm[:,:,2])

  return hyd
end
