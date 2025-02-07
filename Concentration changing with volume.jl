# Import necessary packages
using Catalyst
using DifferentialEquations
using Plots

# Define the reaction network
rn = @reaction_network begin
    k1, A + B --> C
    k2, C --> A + B
end

# Define the volume function
V_func(t) = 5.05 + 4.95 * sin(2 * pi * t / 10 - 0.96)

# Define the initial conditions and parameters
u0 = [1.0, 1.0, 0.0]  # Initial concentrations for A, B, and C
k1 = 1.0
k2 = 0.5
p = [k1, k2]

# Define the initial time span for the simulation
tspan = (0.0, 10.0)

# Define the ODE function to account for volume changes
function ode_func(dy, y, p, t)
    A, B, C = y
    k1, k2 = p

    # Get the current volume
    V = V_func(t)

    # Define the reaction rates
    r1 = k1 * (A / V) * (B / V) * V
    r2 = k2 * (C / V) * V

    # Define the ODEs for concentrations
    dy[1] = -r1 + r2  # dA/dt
    dy[2] = -r1 + r2  # dB/dt
    dy[3] = r1 - r2   # dC/dt
end

# Create the ODE problem
prob = ODEProblem(ode_func, u0, tspan, p)

# Solve the ODE problem
sol = solve(prob, Tsit5())

# Plot the results
plot(sol, vars=(0, 1), label="A", xlabel="Time", ylabel="Concentration", title="Concentrations over Time")
plot!(sol, vars=(0, 2), label="B")
plot!(sol, vars=(0, 3), label="C")

# Plot the volume function
t = range(0, stop=10, length=100)
V = V_func.(t)
plot(t, V, label="Volume", ylabel="Volume", secondary=true)

# Show the plot
display(plot(sol, vars=(0, 1), label="A", xlabel="Time", ylabel="Concentration", title="Concentrations over Time"))
plot!(sol, vars=(0, 2), label="B")
plot!(sol, vars=(0, 3), label="C")
plot!(t, V, label="Volume", ylabel="Volume", secondary=true)