clc
clear 
close all

%iconale stazionale con upwind corretto
T = 2;
dx = 0.01;
x = -1:dx:1;
Nx = length(x);
I = 2:Nx-1;
L = I - 1;
R = I + 1;

u0 = @(x) zeros(length(x),1);
dt = dx/2;
Nt = T/dt;

V = ones(Nx,1);
Vs = zeros(Nx,1);
err = 1;
p = plot(x,V);
ylim([-0.5,1.5])
it = 0;
while err > 1e-6
    it = it + 1;
    dVl = (V(I) - V(L));
    dVr = (V(R) - V(I));
    A = dVl.*(dVl > 0);
    B = dVr.*(dVr < 0);

    Vs(I) = V(I).*(dVl <= 0 & dVr >= 0) +...
        (V(L) + dx).*( dVl > 0 & dVr >= 0) +  (V(R) + dx).*(dVl <= 0 & dVr < 0) +...
        (0.5*(V(R) + V(L) + sqrt(2*V(R).*V(L) - V(R).^2 - V(L).^2 + 2*dx*dx) ) ).*(dVl > 0 & dVr < 0);
    err = norm(Vs - V,inf);
    V = Vs;
    p.YData = Vs;
    pause(0.01);
    drawnow;
end