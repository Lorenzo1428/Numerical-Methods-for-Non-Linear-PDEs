clc 
clear 
close all

f = figure;
T = 5;
dx = 0.01;
x = -2:dx:2;
Nx = length(x);

u0 = @(x) 0.5*exp(-20*x.^2);
fb = @(u) 0.5*u.*u;
fv = @(v) (2/3)*v.^(3/2);

U = u0(x);
V = u0(x).^2;
dt = dx;%/max(abs(uL),abs(uR));
Nt = floor(T/dt);
lambda = dt/dx;

I = 1:Nx;
L = I - 1;
R = I + 1;
L(1) = 1;
R(end) = Nx - 1;

Su = zeros(Nx,1);
Sv = zeros(Nx,1);

p = plot(x,U,x,sqrt(V),x,Su,"--",x,Sv,"--",x, Sv - Su, ":");
xlim([-2,2]);
ylim([-1,2]);
legend("u","v","Su","Sv","diff")

n = 0;
%gudonov
while n < Nt
    n = n+1;
    U(I) = U(I) - lambda*(F(U(I),U(R),fb) - F(U(L),U(I),fb));
    V(I) =  V(I) - lambda*(F(V(I),V(R),fv) - F(V(L),V(I),fv));
    Su = 0.5*(U(R) + U(I));
    Sv = (2/3)*(V(R).^3 + V(I).^3)./(V(R) - V(I));
    p(1).YData = U;
    p(2).YData = sqrt(V);
    p(3).YData = Su;
    p(4).YData = Sv;
    p(5).YData = Sv - Su;
    drawnow;
end

function flux = F(U1,U2,f)
    flux = max(f(U1),f(U2)).*(U1 >= U2) + f(U2).*(U1 <= U2 & U2 < 0) + f(U1).*(U1 < U2 & U1 > 0);  
end