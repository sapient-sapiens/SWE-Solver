Lx=4*pi 
Nx=100
dx=Lx/Nx
x = 0:dx:Lx

Ui = randn(length(x)) 
plot(x, U[i]) 
a=2
D=zeros(Nx, Nx) 
for i = 2:Nx
    D[i, i] = 1
     D[i, i-1] = -1
end 
D[1, 1] = 1
D[1, end] = -1 
D = 1/dx * D 
#u_t = -au_x
T = 5
dt = 0.01
Nt = Int64(T/dt) 
U.=Ui
Ulist = zeros(Nt, Nx) 
for i=1:Nt
    dUdt = -a*D*U
    U. = U + dt*dUdt
    Ulist[i,:].= U
end 


