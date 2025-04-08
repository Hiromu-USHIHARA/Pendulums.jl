include("./ODESolvers.jl")
include("./Pendulum.jl")
using .ODESolvers
using Plots
# %%
const G, L = 9.81, 1.0
u0_ = [π/4, 0.]
tspan_ = (0., 100.)
dt_ = 0.01
# %%
# Orbital in phase space
p2_ = G, L, 0.2, 1.5, 2.0
ts, us = solve(damped_driven_pendulum_force, u0_, tspan_, dt_, p2_)
fig = plot_phase(us)
savefig(fig, "./damped_driven_phase_space.png")
# %%
γs_ = [0.0, 0.1, 0.3, 0.5]
fig = plot()
for (i, γ) in enumerate(γs_)
    p = G, L, γ, 0.,0.
    ts, us = solve(damped_driven_pendulum_force, u0_, tspan_, dt_, p, method=:symplectic_euler)
    plot_phase!(us, label="γ=$(γ)", color=i)
end
savefig(fig, "./damped_phase_space.png")
# %%
# Spectrum
Ωs_ = [2., 10., 20., 50.]
fig = plot()
for (i, Ω) in enumerate(Ωs_)
    p = G, L, .2, 2. ,Ω
    ts, us = solve(damped_driven_pendulum_force, u0_, tspan_, dt_, p)
    plot_spectrum!(ts, us, label="Ω=$(Ω)", color=i)
end
lens!([0., 10.], [0., 30.], inset = (1, bbox(0.2, 0.25, 0.8, 0.6)))
savefig(fig, "./damped_driven_spectrum.png")

