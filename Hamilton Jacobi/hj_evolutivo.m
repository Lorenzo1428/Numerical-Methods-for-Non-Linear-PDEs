clc
clear 
close all


id = input("1. Burgers 2. H(p) = |p| 3. eikonale evolutiva 4. burgers mod \n");
id = max(mod(floor(id),5),1);

T = 1*(id ~= 4) + 2*(id == 4);
dx = 0.004;

b = 1*(id >= 3) + 2*(id < 3);
x = -b:dx:b;
Nx = length(x);
I = 2:Nx-1;
L = I - 1;
R = I + 1;

u0 = @(x) zeros(1,length(x)).*(id >= 3) + abs(x).*(id < 3);
dt = 0.5*dx;
Nt = T/dt;

U = u0(x);

switch id
    case 1
        f = @(u) (0.5*u.*u);
        f_sol = @(x,t) (abs(x) - t*0.5).*(abs(x) > t) + (0.5*abs(x).^2/t).*(abs(x) <= t);
    case 2
        f = @(u) abs(u);
        f_sol = @(x,t) (abs(x) - t).*(abs(x) > t) + (zeros(1,length(x))).*(abs(x) <= t);
    case 3
        f = @(u) abs(u) - 1;
        f_sol = @(x,t) 1 - abs(x);
    case 4
        f = @(u) 0.5*u.^2 - 0.5;
        f_sol = @(x,t) 1 - abs(x);
end

flux = input("1. Lax-F 2. Upwind correct \n");
flux = max(mod(floor(flux),3),1);
F = @(U,V) ( 0.5*(f(U) + f(V)) + 0.5*(dx/dt)*(U - V) ).*(flux == 1) + ( f(max(U,0)) + f(min(V,0)) - f(0)).*(flux ~= 1);  

p = plot(x,U,x,f_sol(x,T));
ylim([-0.2,1.5]);
legend("u approx","u exact");
title("T = 0");

n = 0;
while n < Nt
    n = n+1;
    dUl = (U(I) - U(L))/dx;
    dUr = (U(R) - U(I))/dx;
    U(I) = U(I) - dt*( F(dUl,dUr) );
    U(1) = ( 2*U(2)   - U(3) ).*(id < 3);
    U(end) = ( 2*U(end-1) - U(end-2) ).*(id < 3);
    p(1).YData = U;
    title("T = "+ n*dt);
    pause(0.001);
    drawnow;
end







