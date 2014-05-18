#Diagonal shocktube
function setup_tubexy(self::data2d)

    @assert hyd.x[1] == hyd.y[1]
    @assert hyd.x[end] == hyd.y[end]
    @assert hyd.nx == hyd.ny

    rho1 = 1.0
    rho2 = 0.125
    press1 = 1.0
    press2 = 0.14

    self.rho[:,:] = rho1
    self.press[:,:] = press1
    self.eps[:,:] = press1/rho1/(gamma - 1.0)

    for i = 1:self.ny
        for j = 1:i
                self.rho[i, j] = rho2
                self.press[i, j] = press2
                self.eps[i, j] = press2/rho2/(gamma - 1.0)
            end
        end

    return self
end


#Blast wave
function setup_blast(self::data2d)

    rho1 = 1.0
    press1 = 0.1
    press2 = 10.0

    self.rho[:,:] = rho1*ones(self.ny, self.nx)
    self.press[:,:] = press1*ones(self.ny, self.nx)
    self.eps[:,:] = press1./rho1./(gamma - 1.0)
    self.velx[:,:] = zeros(self.ny, self.nx)
    self.vely[:,:] = zeros(self.ny, self.nx)

    midx = int(self.nx/2)
    midy = int(self.ny/2)


    #circle=[1,3,4,4,5,5]
    #circle=[2,4,5,6,6,6,7,7]
    circle=[3,5,6,7,8,9,9,10,10,10]
    N = length(circle)
    for j = 0:N-1
        yy = j
        xx = circle[N-j]

        self.press[midy+yy, midx-xx:midx+xx] = press2
        self.eps[midy+yy, midx-xx:midx+xx] = press2./rho1./(gamma - 1.0)
        self.press[midy-yy, midx-xx:midx+xx] = press2
        self.eps[midy-yy, midx-xx:midx+xx] = press2./rho1./(gamma - 1.0)
    end

    return self
end

#Two colliding spheres
function setup_collision(self::data2d)

    rho1 = 0.01
    rho2 = 10.0
    press1 = 0.01
    vel2 = 2.0

    self.rho[:,:] = rho1*ones(self.ny, self.nx)
    self.press[:,:] = press1*ones(self.ny, self.nx)
    self.eps[:,:] = press1./rho1./(gamma - 1.0)
    self.velx[:,:] = zeros(self.ny, self.nx)
    self.vely[:,:] = zeros(self.ny, self.nx)

    #first ball
    midx = int(self.nx/2)
    midy = int(0.25self.ny)

    #circle=[1,3,4,4,5,5]
    #circle=[2,4,5,6,6,6,7,7]
    circle=[3,5,6,7,8,9,9,10,10,10]
    N = length(circle)
    for j = 0:N-1
        yy = j
        xx = circle[N-j]

        self.rho[midy+yy, midx-xx:midx+xx] = rho2
        self.rho[midy-yy, midx-xx:midx+xx] = rho2

        self.vely[midy+yy, midx-xx:midx+xx] = vel2
        self.vely[midy-yy, midx-xx:midx+xx] = vel2

    end

    #second ball
    midx = int(self.nx/2)
    midy = int(0.75self.ny)
    #circle=[1,3,4,4,5,5]
    #circle=[2,4,5,6,6,6,7,7]
    circle=[3,5,6,7,8,9,9,10,10,10]
    N = length(circle)
    for j = 0:N-1
        yy = j
        xx = circle[N-j]

        self.rho[midy+yy, midx-xx:midx+xx] = rho2
        self.rho[midy-yy, midx-xx:midx+xx] = rho2

        self.vely[midy+yy, midx-xx:midx+xx] = -vel2
        self.vely[midy-yy, midx-xx:midx+xx] = -vel2
    end

    return self
end

#Free falling spheres
function setup_fall(self::data2d)

    rho1 = 1.0
    rho2 = 2.0
    press1 = 2.5
    press2 = 1.25

    self.rho[:,:] = rho1*ones(self.ny, self.nx)
    self.press[:,:] = press1*ones(self.ny, self.nx)
    self.eps[:,:] = press1./rho1./(gamma - 1.0)
    self.velx[:,:] = zeros(self.ny, self.nx)
    self.vely[:,:] = zeros(self.ny, self.nx)

    #first ball
    midx = int(self.nx/2)
    midy = int(0.8self.ny)

    #circle=[1,3,4,4,5,5]
    #circle=[2,4,5,6,6,6,7,7]
    circle=[3,5,6,7,8,9,9,10,10,10]
    N = length(circle)
    for j = 0:N-1
        yy = j
        xx = circle[N-j]

        self.rho[midy+yy, midx-xx:midx+xx] = rho2
        self.rho[midy-yy, midx-xx:midx+xx] = rho2

        self.press[midy+yy, midx-xx:midx+xx] = press2
        self.press[midy-yy, midx-xx:midx+xx] = press2
    end

    return self
end


#Taylor instability initial data
function setup_taylor(self::data2d)
    rchange = int(0.8self.ny)

    grav = 0.005
    rho1 = 2.0
    rho2 = 1.0
    press1 = 0.01

    self.rho[1:(rchange-1), :] = rho2*ones(rchange-1, self.nx)
    self.rho[rchange:(self.ny), :] = rho1*ones((self.ny-rchange)+1, self.nx)

    self.press[:,:] = press1*ones(self.ny, self.nx)
    self.eps[:,:] = press1./rho1./(gamma - 1.0)

    self.velx[:,:] = zeros(self.ny, self.nx)
    self.vely[:,:] = zeros(self.ny, self.nx)


    offs = 35
    offsy= 30
    dx = self.x[self.nx-offs-1] - self.x[offs]

    function siny(j)
        sin(pi*(self.y[j]-self.y[rchange-offsy])/(self.y[rchange+offsy]-self.y[rchange-offsy]))
    end


#    for j = rchange:self.ny

    #for j = (rchange-offsy):(rchange+offsy)
    for j = (rchange):(rchange+offsy)

        self.vely[j, offs:(self.nx-offs-1)] = Float64[-0.01siny(j)*sin(pi*(self.x[n]-self.x[offs])/dx) for n = offs:(self.nx-offs-1)]
#        self.vely[j, offs:(self.nx-offs-1)] = Float64[-10.0sin(pi*(self.x[n]-self.x[offs])/dx) for n = offs:(self.nx-offs-1)]
    end

    for j = 1:self.ny
        self.press[j,:] = 2.5 - self.rho[j,20]*grav*self.y[j]
    end


    return self
end

#Taylor instability initial data according to Athena code
function setup_taylor2(self::data2d)
    rchange = int(0.5self.ny)

    A = 0.01
    grav = 0.005
    rho1 = 2.0
    rho2 = 1.0

    self.rho[1:(rchange-1), :] = rho2*ones(rchange-1, self.nx)
    self.rho[rchange:(self.ny), :] = rho1*ones((self.ny-rchange)+1, self.nx)

    self.press = zeros(self.ny, self.nx)

    self.velx[:,:] = zeros(self.ny, self.nx)
    self.vely[:,:] = zeros(self.ny, self.nx)

    Lx = abs(self.x[self.nx - self.g] - self.x[self.g+1])
    Ly = abs(self.y[self.ny - self.g] - self.y[self.g+1])

    for j = (self.g+1):(self.ny - self.g), i = (self.g+1):(self.nx - self.g)
        self.vely[j,i] = A*((1.0+cos(2pi*self.x[i]/Lx)*(1.0+cos(2pi*self.y[j]/Ly))))
    end

    #set pressure to hydrostatic equiblibrium
    for j = (self.g+1):(self.ny-self.g)
        self.press[j,:] = 2.5 - self.rho[j,20]*grav
    end
    self.eps[:,:] = self.press[:,:] ./ self.rho[:,:] ./ (gamma - 1.0)

    return self
end


#Kelvin-Helmholtz -instability
function setup_kh(self::data2d)
    
    A = 0.01
    rho1 = 1.0
    rho2 = 2.0
    vel = 0.5

    for j = 1:self.ny, i = 1:self.nx
        if abs(self.y[j]) >= 0.25
            self.rho[j, i] = rho1
            self.velx[j,i] = vel
        else
            self.rho[j, i] = rho2
            self.velx[j,i] = -vel
        end
    end

    Lx = 0.5*abs(self.x[self.nx - self.g] - self.x[self.g+1])
    Ly = 0.5*abs(self.y[self.ny - self.g] - self.y[self.g+1])

    for j = (self.g+1):(self.ny - self.g), i = (self.g+1):(self.nx - self.g)
        self.velx[j,i] += A*((1.0+cos(2pi*self.x[i]/Lx)*(1.0+cos(2pi*self.y[j]/Ly))))
        #self.vely[j,i] += A*((1.0+sin(2pi*self.x[i]/Lx)*(1.0+sin(2pi*self.y[j]/Ly))))
    end

    self.press = ones(self.ny, self.nx)*2.5
    self.vely[:,:] = zeros(self.ny, self.nx)
    self.eps[:,:] = self.press[:,:] ./ self.rho[:,:] ./ (gamma - 1.0)

    return self
end

