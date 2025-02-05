using Pkg
Pkg.add(["Catalyst", "DifferentialEquations", "Plots"])

using Catalyst
using DifferentialEquations
using Plots

# Define the species
#=
Define the species:
Cp is the inactivated species
ACp is the activated species
RNA is the RNA strand  
A is 2-aminoimidazole
CppC is the catalyzed dimer (1,3-di(Cystine-5′-phosphoro)-2-aminoimidazolium)
=#
@variables t
@species Cp(t) RNA(t) CppC(t) ACp(t) A(t)
@parameters K_m v_max Keff Ki V(t)

# Define V_func(t)
V_func(t) = 5.05 + 4.95 * sin(2 * pi * t / 10 - 0.96)

# Define the reactions
rxns = @reaction_network begin
    k1/V_func^2, 2ACp → CppC + A                        # Formation of the catalyst (dimer of Cp)
    mm(CppC, v_max, Keff), CppC → Cp + RNA         # Catalyzed polymer growth
    k3/V_func^2, Cp + A → ACp                               # Formation of Activated Cystine
    k4/V_func, CppC → ACp + Cp                            # Degradation of the Catalyzed dimer
    k5/V_func, ACp → Cp + A                               # Degradation of the Activated Cystine
    k6/V_func, CppC + A → 2ACp                            # Degradation of the catalyzed dimer with 2-aminoimidazole
end

# Initial conditions for the species
u0 = [ACp => 100.0, CppC => 0.0, A => 0.0, RNA => 1.0, Cp => 0.0]

# Parameters for the reactions (rate constants)
k1 = 4.55e3  # h^-1 mM^-1
k3 = 0.1  # h^-1
k4 = 0.167  # h^-1
k5 = 3.33e-3  # h^-1
k6 = 0.183  # h^-1 mM^-1

# Additional parameters
K_m = 1.06  # mM
v_max = 19.5  # h^-1
Ki = 24.7  # mM
Keff = K_m * (1 + (Cp + ACp) / Ki)  # Effective Michaelis constant

p = [k1, k3, k4, k5, k6, K_m, v_max, Ki]

tspan = (0.0, 50.0)

# Create the ODE problem
prob = ODEProblem(rxns, u0, tspan, p)

# Solve the ODE problem
sol = solve(prob, Tsit5())

# Calculate the absolute number of RNA molecules over time
rna_count = sol[RNA]

# Desired polymer size in number of molecules
desired_polymer_size = 50  # This is the target size (number of RNA molecules)

# Find the time it takes to reach the desired RNA size
time_to_reach_size = findfirst(>(desired_polymer_size), rna_count)

# Plot the results
plot(sol, vars=[Cp, RNA, CppC, ACp, A], xlabel="Time", ylabel="Concentration", legend=true)
title!("Polymer Growth with Catalyst")

# Plot the size of the polymer over time
plot(sol.t, rna_count, xlabel="Time", ylabel="Number of RNA Molecules", legend=false)
title!("Number of RNA Molecules Over Time")

# Display the time to reach the desired RNA size
println("Time to reach the desired polymer size ($desired_polymer_size units): $time_to_reach_size")