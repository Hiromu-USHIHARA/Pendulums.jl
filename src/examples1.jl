include("./ODESolvers.jl")
include("./Pendulum.jl")
using .ODESolvers
using Plots
# %%
const G, L = 9.81, 1.0
u0_ = [π/4, 0.]
tspan_ = (0., 50.)
dt_ = 0.01
# %%
# simple
p0_ = G, L
simulate_pendulum(simple_pendulum_force, single_pendulum_position, u0_, tspan_, dt_, p0_, gifname="./simple_pendulum.gif")
# %%
# damped
p1_ = G, L, 0.2, 0.0, 0.0
simulate_pendulum(damped_driven_pendulum_force, single_pendulum_position, u0_, tspan_, dt_, p1_, gifname="damped_pendulum.gif")
# %%
# damped + driven
p2_ = G, L, 0.2, 1.5, 2.0
simulate_pendulum(damped_driven_pendulum_force, single_pendulum_position, u0_, tspan_, dt_, p2_, gifname="damped_driven_pendulum.gif")
# %%
# varying γ
γs_ = [0.0, 0.1, 0.3, 0.5]
A_, Ω_ = 1.5, 2.0
# %%
fig = plot()
for (i, γ) in enumerate(γs_)
    p = G, L, γ, 0.,0.
    ts, us = solve(damped_driven_pendulum_force, u0_, tspan_, dt_, p)
    plot_amplitude!(ts, us, label="γ=$(γ)", color=i)
end
savefig(fig, "./damped_amplitude.png")
# %%
fig = plot()
for (i, γ) in enumerate(γs_)
    p = G, L, γ, A_, Ω_
    ts, us = solve(damped_driven_pendulum_force, u0_, tspan_, dt_, p)
    plot_amplitude!(ts, us, label="γ=$(γ)", color=i)
end
savefig(fig, "./damped_driven_amplitude.png")