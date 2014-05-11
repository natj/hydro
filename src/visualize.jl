function visualize(data::AbstractMatrix, xmin, ymin, xmax, ymax; g=1)

    cm = Uint32[Color.convert(Color.RGB24,c) for c in flipud(Color.colormap("RdBu"))]

    ny, nx = size(data)
    xs = g
    xe = nx - g
    ye = ny - g

    hdata = data[xs:ye, xs:xe]
    p = FramedPlot()
    clims = (minimum(hdata), maximum(hdata))
    img = Winston.data2rgb(hdata, clims, cm)
    add(p, Image((xmin, xmax), (ymin, ymax), img;))
    setattr(p, xrange=(xmin, xmax))
    setattr(p, yrange=(ymin, ymax))
    display(p)
    return p
end

function visualize(hyd::data2d)

    cm = Uint32[Color.convert(Color.RGB24,c) for c in flipud(Color.colormap("RdBu"))]
    #cm = Uint32[Color.convert(Color.RGB24,c) for c in Color.colormap("RdBu")]


    xs = hyd.g
    xe = hyd.nx-hyd.g-1
    ye = hyd.ny-hyd.g-1


    hdata = hyd.rho[xs:ye, xs:xe]
    p1=FramedPlot()
    #clims = (minimum(hdata), maximum(hdata))
    clims = (0.0, 2.1)
    #clims = (0.0, maximum(hdata))
    img = Winston.data2rgb(hdata, clims, cm)
    add(p1, Image((hyd.x[xs], hyd.x[xe]), (hyd.y[xs], hyd.y[ye]), img;))
    setattr(p1, xrange=(hyd.x[xs], hyd.x[xe]))
    setattr(p1, yrange=(hyd.y[xs], hyd.y[ye]))
    setattr(p1, title="rho")

    #additional window for 1-dim figs
    xy = Float64[sqrt(hyd.x[i]^2+hyd.y[i]^2.0) for i = 1:hyd.nx]
    rhoxy = Float64[hyd.rho[hyd.nx-i+1,i] for i = hyd.nx:-1:1]
    velxy = Float64[sqrt(hyd.velx[hyd.nx-i+1,i]^2.0 + hyd.velx[hyd.nx-i+1,i]^2.0) for i = hyd.nx:-1:1]

    velxxy = Float64[hyd.velx[hyd.nx-i+1,i] for i = hyd.nx:-1:1]
    velyxy = Float64[hyd.vely[hyd.nx-i+1,i] for i = hyd.nx:-1:1]

    pressxy = Float64[hyd.press[hyd.nx-i+1,i] for i = hyd.nx:-1:1]
    epsxy = Float64[hyd.eps[hyd.nx-i+1,i] for i = hyd.nx:-1:1]

    #p11 = plot(hyd.y, rhoo)#, yrange=[0.0, 3.0])
    #p11 = plot(hyd.y, hyd.rho[:,50])#, yrange=[0.0, 3.0])
    #p11 = oplot(hyd.y, hyd.rho[:,3], "r--")

    p11 =  plot(xy, rhoxy, "r;")
    #p11 = oplot(xy, velxy, "b-")

    p11 = oplot(xy, abs(velxxy), "b;")
    p11 = oplot(xy, abs(velyxy), "c;")

    p11 = oplot(xy, pressxy, "g;")
    #p11 = oplot(hyd.y, epsxy, "k")

    #mid = int(hyd.ny/2)
    #p11 = oplot(hyd.y, hyd.rho[:,mid], "r")
    #p11 = oplot(hyd.y, abs(hyd.vely[:,mid]), "b")
    #p11 = oplot(hyd.y, hyd.press[:,mid], "g")
    #p11 = oplot(hyd.y, hyd.eps[:,mid], "k")

    #pressure
    hdata = hyd.press[xs:ye, xs:xe]
    p2=FramedPlot()
    #clims = (minimum(hdata), maximum(hdata))
    #clims = (0.0, 1.0)
    clims = (0.0, maximum(hdata))
    img = Winston.data2rgb(hdata, clims, cm)
    add(p2, Image((hyd.x[xs], hyd.x[xe]), (hyd.y[xs], hyd.y[ye]), img;))
    setattr(p2, xrange=(hyd.x[xs], hyd.x[xe]))
    setattr(p2, yrange=(hyd.y[xs], hyd.y[ye]))
    setattr(p2, title="press")

    #vel
    hdata = sqrt(hyd.velx.^2.0 .+ hyd.vely.^2.0)[xs:ye, xs:xe]
    p3=FramedPlot()
    #clims = (minimum(hdata), maximum(hdata))
    clims = (0.0, maximum(hdata))
    #clims = (0.0, 1.0)
    img = Winston.data2rgb(hdata, clims, cm)
    add(p3, Image((hyd.x[xs], hyd.x[xe]), (hyd.y[xs], hyd.y[ye]), img;))
    setattr(p3, xrange=(hyd.x[xs], hyd.x[xe]))
    setattr(p3, yrange=(hyd.y[xs], hyd.y[ye]))
    setattr(p3, title="vel")

    #eps
    hdata = hyd.eps[xs:ye, xs:xe]
    p4=FramedPlot()
    #clims = (minimum(hdata), maximum(hdata))
    #clims = (0.0, 1.0)
    clims = (0.0, maximum(hdata))
    img = Winston.data2rgb(hdata, clims, cm)
    add(p4, Image((hyd.x[xs], hyd.x[xe]), (hyd.y[xs], hyd.y[ye]), img;))
    setattr(p4, xrange=(hyd.x[xs], hyd.x[xe]))
    setattr(p4, yrange=(hyd.y[xs], hyd.y[ye]))
    setattr(p4, title="eps")

    t = Table(2,2)
    t[1,1] = p1
    #t[1,2] = p2
    t[1,2] = p11
    t[2,1] = p3
    t[2,2] = p4
    display(t)

end
