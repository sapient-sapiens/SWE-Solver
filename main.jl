module MainModule 

include("init.jl")

using .ShallowWaterEquation

export run_simulation, central_diff_x, central_diff_y, rhs, rk4_step

function central_diff_x(f, dx, nx, ny)
    f_im, f_ip = circshift(f, (1, 0)), circshift(f, (-1, 0))
    return (f_ip - f_im) / (2 * dx)
end

function central_diff_y(f, dy, nx, ny)
    f_jm, f_jp = circshift(f, (0, 1)), circshift(f, (0, -1))
    return (f_jp - f_jm) / (2 * dy)
end

function rhs(h, u, v, dx, dy, nx, ny, g, f)
    dh_dx, dh_dy = central_diff_x(h, dx, nx, ny), central_diff_y(h, dy, nx, ny)
    du_dx, du_dy = central_diff_x(u, dx, nx, ny), central_diff_y(u, dy, nx, ny)
    dv_dx, dv_dy = central_diff_x(v, dx, nx, ny), central_diff_y(v, dy, nx, ny)
    
    hu, hv = h .* u, h .* v
    d_hu_dx, d_hv_dy = central_diff_x(hu, dx, nx, ny), central_diff_y(hv, dy, nx, ny)
    
    dh_dt = -(d_hu_dx + d_hv_dy)
    du_dt = -u .* du_dx - v .* du_dy + f .* v - g .* dh_dx
    dv_dt = -u .* dv_dx - v .* dv_dy - f .* u - g .* dh_dy
    
    return dh_dt, du_dt, dv_dt
end

function rk4_step(h, u, v, dx, dy, nx, ny, g, f, dt)
    k1_h, k1_u, k1_v = rhs(h, u, v, dx, dy, nx, ny, g, f)
    
    h2, u2, v2 = h + 0.5*dt*k1_h, u + 0.5*dt*k1_u, v + 0.5*dt*k1_v
    k2_h, k2_u, k2_v = rhs(h2, u2, v2, dx, dy, nx, ny, g, f)
    
    h3, u3, v3 = h + 0.5*dt*k2_h, u + 0.5*dt*k2_u, v + 0.5*dt*k2_v
    k3_h, k3_u, k3_v = rhs(h3, u3, v3, dx, dy, nx, ny, g, f)
    
    h4, u4, v4 = h + dt*k3_h, u + dt*k3_u, v + dt*k3_v
    k4_h, k4_u, k4_v = rhs(h4, u4, v4, dx, dy, nx, ny, g, f)
    
    h_new = h + (dt/6.0)*(k1_h + 2.0*k2_h + 2.0*k3_h + k4_h)
    u_new = u + (dt/6.0)*(k1_u + 2.0*k2_u + 2.0*k3_u + k4_u)
    v_new = v + (dt/6.0)*(k1_v + 2.0*k2_v + 2.0*k3_v + k4_v)
    
    return h_new, u_new, v_new
end

function run_simulation(nx, ny, Lx, Ly, T_Final, dt; g = 1.0, f = 1.0)
    x, y, dx, dy = get_grid(nx, ny, Lx, Ly) 
    h, u, v = rossby_wave_initialization(x, y, nx, ny)
    steps = Int(T_Final / dt)
    h_series, u_series, v_series = zeros(nx, ny, steps), zeros(nx, ny, steps), zeros(nx, ny, steps)
    t = 0.0 
    t_series = zeros(steps)
    for step in 1:steps
        h, u, v = rk4_step(h, u, v, dx, dy, nx, ny, g, f, dt) 
        h_series[:, :, step] = h 
        u_series[:, :, step] = u 
        v_series[:, :, step] = v 
        t += dt 
        t_series[step] = t
    end 
    return x, y, t_series, h_series, u_series, v_series
end 

end 