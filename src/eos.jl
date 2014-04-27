#equation of states

#gamma law
function eos_press(rho, eps, gamma)
    const pmin = 1.0e-10
    const pmax = 1.0e20

    press = (gamma - 1.0) .* rho .* eps
    press = clamp(press, pmin, pmax)
    return press
end

#soundspeed
function eos_cs2(rho, eps, gamma)
    prs = (gamma - 1.0) .* rho .* eps
    dpde = (gamma - 1.0) .* rho
    dpdrho = (gamma - 1.0) .* eps
    cs2 = dpdrho .+ dpde .* prs ./ (rho + 1.0e-30).^2.0


    for i = 1:length(cs2)
        if cs2[i] < 0.0
            println("cs2 = $(cs2[i]) i=$i")
            println(eps[i])
            #println(dpdrho[i])
            #println(dpde[i])
            #println(prs[i])
        end
    end


    return cs2
end