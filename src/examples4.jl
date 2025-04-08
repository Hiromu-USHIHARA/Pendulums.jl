include("./ODESolvers.jl")
include("./Pendulum.jl")
using .ODESolvers
using Plots
# %%
const G, L = 9.81, 1.0
u0_ = [0.001, 0., 0.001, 0.]
tspan_ = (0., 50.)
dt_ = 0.005
p0_ = G, L
# %%
ts_, us_ = solve(double_pendulum_force, u0_, tspan_, dt_, p0_, method=:symplectic_euler)
θ₁_num = [u[1] for u in us_]
θ₂_num = [u[3] for u in us_]
# %%
Ω₊, Ω₋ = sqrt(2 + sqrt(2)) * sqrt(G / L), sqrt(2 - sqrt(2)) * sqrt(G / L)
C₊, C₋ = (u0_[1] - u0_[3]/sqrt(2))/2, (u0_[1] + u0_[3]/sqrt(2))/2
# %%
θ₁_analytical = [C₊ * cos(Ω₊ * t) + C₋ * cos(Ω₋ * t) for t in ts_]
θ₂_analytical = [-sqrt(2)*C₊ * cos(Ω₊ * t) + sqrt(2)*C₋ * cos(Ω₋ * t) for t in ts_]
# %%
fig = plot(
    xlabel="Time", ylabel="Angle"
)
plot!(ts_, θ₁_num, label="θ₁ (Numerical)", lw=2)
plot!(ts_, θ₁_analytical, label="θ₁ (Analytical)", lw=2, ls=:dash)
savefig(fig, "./double_pendulum.linear_approximation1.png")
fig = plot(
    xlabel="Time", ylabel="Angle"
)
plot!(ts_, θ₂_num, label="θ₂ (numerical)", lw=2)
plot!(ts_, θ₂_analytical, label="θ₂ (theoretical)", ls=:dash, lw=2)
savefig(fig, "./double_pendulum.linear_approximation2.png")