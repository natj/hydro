#Shoctube initial data
function setup_tubey(self::data2d)
    rchange = 0.5(self.y[self.ny - self.g] - self.y[self.g + 1])

    rho1 = 1.0
    rho2 = 0.125
    press1 = 1.0
    press2 = 0.1


    for i = 1:self.ny, j = 1:self.nx
        if self.y[i] < rchange
            self.rho[i, j] = rho1
            self.press[i ,j] = press1
            self.eps[i, j] = press1/rho1/(gamma - 1.0)
            self.velx[i, j] = 0.0
            self.vely[i, j] = 0.0
        else
            self.rho[i, j] = rho2
            self.press[i, j] = press2
            self.eps[i, j] = press2/rho2/(gamma - 1.0)
            self.velx[i, j] = 0.0
            self.vely[i, j] = 0.0
        end
    end

    return self
end

#Shoctube initial data
function setup_tubex(self::data2d)
    rchange = 0.5(self.x[self.nx - self.g] - self.x[self.g + 1])

    rho1 = 1.0
    rho2 = 0.125
    press1 = 1.0
    press2 = 0.1

    for i = 1:self.ny, j = 1:self.nx
        if self.x[j] < rchange
            self.rho[i, j] = rho1
            self.press[i ,j] = press1
            self.eps[i, j] = press1/rho1/(gamma - 1.0)
            self.velx[i, j] = 0.0
            self.vely[i, j] = 0.0
        else
            self.rho[i, j] = rho2
            self.press[i, j] = press2
            self.eps[i, j] = press2/rho2/(gamma - 1.0)
            self.velx[i, j] = 0.0
            self.vely[i, j] = 0.0
        end
    end

    return self
end

