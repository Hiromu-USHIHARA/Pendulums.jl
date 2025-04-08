include("./ODESolvers.jl")
using .ODESolvers
using Plots
# %%
"""
    simple_pendulum_force(t, u, p)

Force on single pendulum (simple harmonic oscillator)
"""
function simple_pendulum_force(t, u, p)
    θ, ω = u
    g, l = p
    return ω, -g/l * sin(θ)
end
# %%
"""
    single_pendulum_position(u)

Position of single pendulum
"""
function single_pendulum_position(u, p)
    θ, ω = u
    g, l = p
    return l*sin(θ), -l*cos(θ)
end
# %%
"""
    simulate_pendulum(force, to_position, u0, tspan, dt, p; gifname="anim.gif", is_single=true)

- `force(t,u,p)`        : Force on position `u` at time `t` with parameters `p`
- `to_position(u,p)`    : Conversion to (x,y)
- `u0`                  : Initial condition
- `tspan`               : Time range (`t0`, `tmax`)
- `dt`                  : Time step
- `p`                   : Parameters (g, L)

Returns: Animation GIF
"""
function simulate_pendulum(force, to_position, u0, tspan, dt, p; gifname="anim.gif", is_single=true)
    if is_single
        return simulate_pendulum_single(force, to_position, u0, tspan, dt, p; gifname=gifname)
    else
        return simulate_pendulum_multiple(force, to_position, u0, tspan, dt, p; gifname=gifname)
    end
end
# %%
function simulate_pendulum_single(force, to_position, u0, tspan, dt, p; gifname="anim.gif")
    ts, us = solve(force, u0, tspan, dt, p)
    positions=[to_position(u,p) for u in us]
    xs=[pos[1] for pos in positions]
    ys=[pos[2] for pos in positions]
    #
    l = p[2]
    #
    anim = Animation()
    for i in 1:10:length(xs)
        plot(
            [0.0, xs[i]], [0.0, ys[i]], lw=3,
            xlims=(-1.2*l, 1.2*l),
            ylims=(-1.2*l, 0.2),
            legend=nothing, aspect_ratio=1., title="t=$(ts[i])"
        )
        scatter!([xs[i]], [ys[i]], markersize=8, color=:orange)
        frame(anim)
    end
    #
    gif(anim, gifname, fps=100)
end
# %%
function simulate_pendulum_multiple(force, to_position, u0, tspan, dt, p; gifname="anim.gif")
    ts, us = solve(force, u0, tspan, dt, p)
    positions=[to_position(u,p) for u in us] # [((x1, y1), (x2, y2), ...)(t)]
    Nparticles=length(positions[1])
    #
    xss = [[pos[i][1] for pos in positions] for i in 1:Nparticles]
    yss = [[pos[i][2] for pos in positions] for i in 1:Nparticles]
    #
    l = p[2]
    Nl = Nparticles*l
    #
    anim=Animation()
    for i in 1:10:length(ts)
        xvals = [0.; [xss[j][i] for j in 1:length(xss)]]
        yvals = [0.; [yss[j][i] for j in 1:length(yss)]]
        #
        plot(
            xvals, yvals, lw=3,
            xlims=(-1.2*Nl, 1.2*Nl),
            ylims=(-1.2*Nl, 1.2*Nl),
            legend=nothing, aspect_ratio=1., title="t=$(round(ts[i], digits=2))"
        )
        scatter!([xvals[end]], [yvals[end]], markersize=8, color=:orange)
        frame(anim)
    end
    #
    gif(anim, gifname, fps=100)
end
# %%
"""
    damped_driven_pendulum_force(t, u, p)

Force on single pendulum with damping and driving forces
p = (g, l, γ, A, Ω)
"""
function damped_driven_pendulum_force(t, u, p)
    θ, ω = u
    g, l, γ, A, Ω = p
    return ω, -γ*ω - g/l *sin(θ) + A*cos(Ω*t)
end
# %%
"""
    plot_amplitude(ts, us; label="", color=1)

Time dependence of amplitude
"""
function plot_amplitude(ts, us; label="", color=1)
    θs = [u[1] for u in us]
    return plot(ts, θs, label=label, xlabel="Time", ylabel="Angle", lw=2, color=color)
end
# %%
function plot_amplitude!(ts, us; label="", color=1)
    θs = [u[1] for u in us]
    plot!(ts, θs, label=label, xlabel="Time", ylabel="Angle", lw=2, color=color)
end
# %%
"""
    plot_phase(us; xlabel="θ", ylabel="ω", label="", color=1)
"""
function plot_phase(us; xlabel="θ", ylabel="ω", label="", color=1)
    θs = [u[1] for u in us]
    ωs = [u[2] for u in us]
    return plot(θs, ωs, xlabel=xlabel, ylabel=ylabel, label=label, color=color)
end
# %%
function plot_phase!(us; xlabel="θ", ylabel="ω", label="", color=1)
    θs = [u[1] for u in us]
    ωs = [u[2] for u in us]
    plot!(θs, ωs, xlabel=xlabel, ylabel=ylabel, label=label, color=color)
end
# %%
using FFTW
# %%
"""
    plot_spectrum(ts, us; dt=nothing, label="", color=1)
"""
function plot_spectrum(ts, us; dt=nothing, label="", color=1)
    θs=[u[1] for u in us]
    N=length(θs)
    #
    if isnothing(dt)
        dt = ts[2] - ts[1]
    end
    #
    fs = fft(θs .- sum(θs)/N)
    freqs = (0:N-1) ./ (N*dt)
    #
    return plot(
        freqs[1:div(N,2)], abs.(fs[1:div(N,2)]),
        xlabel="Frequency", ylabel="Amplitude", lw=2, label=label, color=color
    )
end
# %%
function plot_spectrum!(ts, us; dt=nothing, label="", color=1)
    θs=[u[1] for u in us]
    N=length(θs)
    #
    if isnothing(dt)
        dt = ts[2] - ts[1]
    end
    #
    fs = fft(θs .- sum(θs)/N)
    freqs = (0:N-1) ./ (N*dt)
    #
    plot!(
        freqs[1:div(N,2)], abs.(fs[1:div(N,2)]),
        xlabel="Frequency", ylabel="Amplitude", lw=2, label=label, color=color
    )
end
# %%
"""
    double_pendulum_force(t, u, p)

Force on double pendulum
- `u`   : [θ₁, ω₁, θ₂, ω₂]
- `p`   : (g, l)
"""
function double_pendulum_force(t, u, p)
    θ₁, ω₁, θ₂, ω₂ = u
    # θ₁, θ₂ = θ₁%2π, θ₂%2π
    g, l = p
    #
    Δ = θ₁ - θ₂
    denominator = 2. - cos(Δ)
    #
    dθ₁, dθ₂ = ω₁, ω₂
    #
    dω₁ = (-ω₁^2*sin(Δ)*cos(Δ)-ω₂^2*sin(Δ)+g/l*(sin(θ₂)*cos(Δ)-2*sin(θ₁))) / denominator
    dω₂ = (ω₂^2*cos(Δ)*sin(Δ)+2*ω₁^2*sin(Δ) + 2*g/l*(sin(θ₁)*cos(Δ)-sin(θ₂))) / denominator
    #
    return [dθ₁, dω₁, dθ₂, dω₂]
end
# %%
"""
    double_pendulum_position(u, p)
"""
function double_pendulum_position(u, p)
    θ₁, _, θ₂, _ = u
    _, l = p
    #
    x₁ = l*sin(θ₁)
    y₁ = -l*cos(θ₁)
    x₂ = x₁ + l*sin(θ₂)
    y₂ = y₁ - l*cos(θ₂)
    #
    return (x₁, y₁), (x₂, y₂)
end
# %%
"""
    simulate_pendulum_multiple_w_trace(force, to_position, u0, tspan, dt, p; gifname="anim.gif")
"""
function simulate_pendulum_multiple_w_trace(force, to_position, u0, tspan, dt, p; gifname="anim.gif")
    ts, us = solve(force, u0, tspan, dt, p, method=:symplectic_euler)
    positions=[to_position(u,p) for u in us] # [((x1, y1), (x2, y2), ...)(t)]
    Nparticles=length(positions[1])
    #
    xss = [[pos[i][1] for pos in positions] for i in 1:Nparticles]
    yss = [[pos[i][2] for pos in positions] for i in 1:Nparticles]
    #
    l = p[2]
    Nl = Nparticles*l
    trace_x, trace_y = Float64[], Float64[]
    #
    anim=Animation()
    for i in 1:10:length(ts)
        xvals = [0.; [xss[j][i] for j in 1:length(xss)]]
        yvals = [0.; [yss[j][i] for j in 1:length(yss)]]
        push!(trace_x, xvals[end])
        push!(trace_y, yvals[end])
        #
        plot(
            xvals, yvals, lw=3,
            xlims=(-1.2*Nl, 1.2*Nl),
            ylims=(-1.2*Nl, 1.2*Nl),
            legend=nothing, aspect_ratio=1., title="t=$(round(ts[i], digits=2))"
        )
        plot!(trace_x, trace_y, lc=:gray, lw=1, label="")
        scatter!([xvals[end]], [yvals[end]], markersize=8, color=:orange)
        frame(anim)
    end
    #
    gif(anim, gifname, fps=100)
end
# %%
