#equation of states

#ideal gas
function eos_press(rho, eps, gamma)
    press = (gamma - 1.0) .* rho .* eps
    return press
end

#soundspeed
function eos_cs2(rho, eps, gamma)
    prs = (gamma - 1.0) .* rho .* eps
    dpde = (gamma - 1.0) .* rho
    dpdrho = (gamma - 1.0) .* eps
    cs2 = dpdrho .+ dpde .* prs ./ (rho + 1.0e-30).^2.0
    return cs2
end