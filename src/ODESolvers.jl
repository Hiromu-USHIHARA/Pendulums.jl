module ODESolvers
export solve
"""
    solve(f, u0, tspan, dt, p; method=:rk4)

ODE solver
"""
function solve(f, u0, tspan, dt, p; method=:rk4)
    if method==:rk4
        return solve_rk4(f, u0, tspan, dt, p)
    elseif method==:symplectic_euler
        return solve_symplectic_euler(f, u0, tspan, dt, p)
    else
        error("Unsupported method: $method")
    end
end
# %%
"""
    solve_rk4(f, u0, tspan, dt, p)

ODE solver with 4th order Runge--Kutta algorithm
"""
function solve_rk4(f, u0, tspan, dt, p)
    t0, tmax=tspan
    ts=t0:dt:tmax
    us=Vector{typeof(u0)}(undef, length(ts))
    us[begin]=copy(u0)
    #
    for i in 1:length(ts)-1
        t=ts[i]
        u=us[i]
        #
        k1=dt .*f(t,u,p)
        k2=dt .*f(t+dt/2,u .+ k1 ./2, p)
        k3=dt .*f(t+dt/2,u .+ k2 ./2, p)
        k4=dt .*f(t+dt, u .+ k3, p)
        #
        us[i+1]=u .+ (k1 .+ 2.0 .*k2 .+  2.0 .*k3 .+ k4) ./ 6.
    end
    return ts, us
end 
# %%
"""
    solve_symplectic_euler(f, u0, tspan, dt, p)

ODE solver with the sympelctic Euler method for u=[x1, v1, x2, v2, ...]
"""
function solve_symplectic_euler(f, u0, tspan, dt, p)
    t0, tmax=tspan
    ts=t0:dt:tmax
    us=Vector{typeof(u0)}(undef, length(ts))
    us[begin]=copy(u0)
    N=length(u0)÷2
    #
    for i in 1:length(ts)-1
        t=ts[i]
        u=copy(us[i])
        #
        du=f(t,u,p)
        u_new=copy(u)
        for j in 1:N
            a_j = du[2j]
            u_new[2j] += dt*a_j
        end
        for j in 1:N
            u_new[2j-1] += dt*u_new[2j]
        end
        #
        us[i+1]=u_new
    end
    return ts, us
end 
end