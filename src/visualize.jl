#Visualization

#Colormap
cm = Uint32[Color.convert(Color.RGB24,c) for c in flipud(Color.colormap("RdBu"))]
#cm = Uint32[Color.convert(Color.RGB24,c) for c in Color.colormap("RdBu")]

#Visualize common matrix
function visualize_M(data::AbstractMatrix, 
                     x, y;
                     g=1, 
                     title="", 
                     dmin=0.0, 
                     dmax=:find)

    ny, nx = size(data)
    xs = g +1 
    xe = nx - g + 1
    ye = ny - g + 1

    xmin = x[xs]
    xmax = x[xe]
    ymin =  y[xs]
    ymax =  y[ye]    

    hdata = data[xs:ye, xs:xe]
    p = FramedPlot()

    #look or set min and max values
    if dmin == :find
        hmin = minimum(hdata)
    elseif isa(dmin, Real)
        hmin = dmin
    end

    if dmax == :find
        hmax = maximum(hdata)
    elseif isa(dmax, Real)
        hmax = dmax
    end

    clims = (minimum(hdata), maximum(hdata))

    #make image
    img = Winston.data2rgb(hdata, clims, cm)
    add(p, Image((xmin, xmax), (ymin, ymax), img;))
    setattr(p, xrange=(xmin, xmax))
    setattr(p, yrange=(ymin, ymax))
    setattr(p, title=title)

    return p
end

function visualize(hyd::data2d)

    #density
    p1 = visualize_M(hyd.rho,
                     hyd.x, hyd.y,
                     g=3,
                     title="density <i>\\rho</i>",
                     dmin = :find,
                     dmax = :find)

    #pressure
    p2 = visualize_M(hyd.press,
                     hyd.x, hyd.y,
                     g=3,
                     title="pressure <i>p</i>",
                     dmin = 0.0,
                     dmax = :find)


    #velocity
    hdata = sqrt(hyd.velx.^2.0 .+ hyd.vely.^2.0)
    #hdata = hyd.velx[xs:ye, xs:xe]

    p3 = visualize_M(hdata,
                     hyd.x, hyd.y,
                     g=3,
                     title="velocity <i>(u^2 + v^2)^{1/2}</i>",
                     dmin = :find,
                     dmax = :find)

    #internal energy
    p4 = visualize_M(hyd.press,
                     hyd.x, hyd.y,
                     g=3,
                     title="internal energy <i>ε</i>",
                     dmin = 0.0,
                     dmax = :find)


    ######
    #additional window for 1-dim figs

    #diagonal slice
    xy = Float64[sqrt(hyd.x[i]^2+hyd.y[i]^2.0) for i = 1:hyd.nx]
    rhoxy = Float64[hyd.rho[hyd.nx-i+1,i] for i = hyd.nx:-1:1]
    velxy = Float64[sqrt(hyd.velx[hyd.nx-i+1,i]^2.0 + hyd.velx[hyd.nx-i+1,i]^2.0) for i = hyd.nx:-1:1]

    velxxy = Float64[hyd.velx[hyd.nx-i+1,i] for i = hyd.nx:-1:1]
    velyxy = Float64[hyd.vely[hyd.nx-i+1,i] for i = hyd.nx:-1:1]

    pressxy = Float64[hyd.press[hyd.nx-i+1,i] for i = hyd.nx:-1:1]
    epsxy = Float64[hyd.eps[hyd.nx-i+1,i] for i = hyd.nx:-1:1]

    p11 =  plot(xy, rhoxy, "r",title="diagonal slice")
    p11 = oplot(xy, abs(velxxy), "b")
    p11 = oplot(xy, abs(velyxy), "c")
    p11 = oplot(xy, pressxy, "g")
    #p11 = oplot(hyd.y, epsxy, "k")


    #middle slice
    mid = int(hyd.nx/2)
    p12 = plot(hyd.y, hyd.rho[:,mid], "r",title="middle slice")
    p12 = oplot(hyd.y, hyd.vely[:,mid], "b")
    p12 = oplot(hyd.y, hyd.press[:,mid], "g")
    #p12 = oplot(hyd.y, hyd.eps[:,mid], "k")


    t = Table(2,3)
    t[1,1] = p1
    t[1,2] = p2
    t[2,1] = p3
    t[2,2] = p4
    
    t[1,3] = p11
    t[2,3] = p12

    display(t)

end
