module ShallowWaterEquation 

export rossby_wave_initialization, random_initialization, get_grid

function get_grid(nx::Int, ny::Int, Lx::Float64, Ly::Float64) 
    dx = Lx/nx 
    dy = Ly/ny
    x = range(dx/2, Lx-dx/2, length = nx) 
    y = range(dy/2, Ly -dy/2, length = ny) 
    return x, y, dx, dy  
end

function rossby_wave_initialization(x, y, nx, ny)
    h = zeros(nx, ny) 
    u = zeros(nx, ny) 
    v = zeros(nx, ny) 
    k = 2*pi/last(x) 
    h0 = 10.0 
    
    for j in 1:ny
        for i in 1:nx
            h[i, j] = h0 + 1.0*sin(k*x[i])*cos(k*y[j])
            u[i, j] = -1.0*k*sin(k*x[i])*sin(k*y[j])
            v[i, j] = 1.0*k*cos(k*x[i])*cos(k*y[j])
        end 
    end 
    return h, u, v 
end 

function random_initialization(x, y, nx, ny) 
    h = zeros(nx, ny) 
    u = zeros(nx, ny) 
    v = zeros(nx, ny) 
    amplitude = 5 
    for i in 1:nx 
        for j in 1:ny 
            h[i, j] = 10.0 + amplitude*randn() 
            u[i, j] = amplitude*randn() 
            v[i, j] = amplitude*randn()
        end
    end
    return h, u, v 
end

end
