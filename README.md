# Hydro #
#### Fairly simple & very modular 1- and 2-dimensional hydrodynamic code

Hydro is a 2 dimensional hydrodynamic code written completely in julia.
It aims to be modular and flexible foundation for more complex computational needs.

## Tests & Showcase

### Kelvin-Helmholtz instability
[![kh](http://natj.github.io/hydro/img/kh.png)](https://vimeo.com/95607699)

## Technical specs
Like mentioned, hydro is widely modular and currently sports these under its hood:


### rsolvers.jl (Riemann solvers)
* HLLE (Harten, Lax, van Leer, Einfeldt)
* unsplitted HLLC (Harten, Lax, van Leer, Contact)


### integrators.jl (time stepping)
* Second order Runge-Kutta (RK2)


### reconstruction.jl
* linear piecewise interpolation between cell centered values to cell edges

#### limiters
* minmod [Roe 1986]
* MC (monotoniced central) [van Leer 1977]


### boundaries.jl (boundary conditions)
* outflow
* reflective
* periodic (toroidal)


### eos.jl (equation of states)
* basic gamma law equation of state


### gravity.jl
* constant y-direction
* self-gravity


### grid.jl
* rectangular 


### viscosity.jl
* artificial viscosity



## Inspiration, influences & references

* Initial code structure: Caltech course [ay190](http://www.tapir.caltech.edu/~cott/ay190/) - Computational astrophysics (by Christian Ott)
* Modifications & additions similar to [pyro2](http://bender.astro.sunysb.edu/hydro_by_example/index.html) code (unsplit method, artificial viscosity, gravity)
* [Tests](http://www.astro.virginia.edu/VITA/ATHENA/athena_testsuite.html) inspired by [Athena](https://trac.princeton.edu/Athena/) code
