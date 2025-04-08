include("./ODESolvers.jl")
include("./Pendulum.jl")
using .ODESolvers
using Plots
# %%
const G, L = 9.81, 1.0
u0_ = [π/2, 0., π, 0.]
tspan_ = (0., 50.)
dt_ = 0.01
# %%
# double pendulum
p0_ = G, L
simulate_pendulum(double_pendulum_force, double_pendulum_position, u0_, tspan_, dt_, p0_, gifname="./double_pendulum.gif", is_single=false)
#
simulate_pendulum_multiple_w_trace(double_pendulum_force, double_pendulum_position, u0_, tspan_, dt_, p0_, gifname="./double_pendulum_trace.gif")
# %%
# initial condition sensitivity
dt_ = 0.01
tspan_ = (0., 50.)
u01_ = [π/2, 0., π + 1.0e-5, 0.]
u02_ = [π/2, 0., π - 1.0e-5, 0.]
ts_, us1_ = solve(double_pendulum_force, u01_, tspan_, dt_, p0_, method=:symplectic_euler)
_, us2_ = solve(double_pendulum_force, u02_, tspan_, dt_, p0_, method=:symplectic_euler)
positions1_ = [double_pendulum_position(u, p0_) for u in us1_]
positions2_ = [double_pendulum_position(u, p0_) for u in us2_]
xss1_ = [[pos[i][1] for pos in positions1_] for i in 1:2]
yss1_ = [[pos[i][2] for pos in positions1_] for i in 1:2]
xss2_ = [[pos[i][1] for pos in positions2_] for i in 1:2]
yss2_ = [[pos[i][2] for pos in positions2_] for i in 1:2]
anim=Animation()
for i in 1:10:length(ts_)
    xvals1_=[0.; [xss1_[j][i] for j in 1:length(xss1_)]]
    yvals1_=[0.; [yss1_[j][i] for j in 1:length(yss1_)]]
    xvals2_=[0.; [xss2_[j][i] for j in 1:length(xss2_)]]
    yvals2_=[0.; [yss2_[j][i] for j in 1:length(yss2_)]]
    plot(
        legend=nothing, aspect_ratio=1., title="t=$(round(ts_[i], digits=2))",
        xlims=(-2.4*p0_[2], 2.4*p0_[2]),
        ylims=(-2.4*p0_[2], 2.4*p0_[2])
    )
    plot!(
        xvals1_, yvals1_, lw=3, color=1
    )
    scatter!([xvals1_[end]], [yvals1_[end]], markersize=8, color=1)
    plot!(
        xvals2_, yvals2_, lw=3, color=2
    )
    scatter!([xvals2_[end]], [yvals2_[end]], markersize=8, color=2)
    #
    frame(anim)
end
gif(anim, "./double_pendulum_initial_condition.gif", fps=1000)