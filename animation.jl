include("main.jl")

using .MainModule
using Plots 

Lx, Ly = 1.0, 1.0
g = 1.0
f = 1.0
dt = 0.001
nx, ny = 50, 50
T_final = 0.5

x, y, t_series, h_series, u_series, v_series = 
    run_simulation(nx, ny, Lx, Ly, T_final, dt; g=g, f=f)

# Create animation
anim = Animation()
for i in 1:max(1, div(size(h_series, 3), 20)):size(h_series, 3)
    frame(anim, heatmap(x, y, h_series[:,:,i]', title="Height (t = $(round(t_series[i], digits=3)))",
                       xlabel="x", ylabel="y", aspect_ratio=:equal, color=:viridis))
end

gif(anim, "shallow_water_animations.gif", fps=10) 