clc
clear 
close all

T = 1;
dx = 0.001;
x = -1:dx:1;
Nx = length(x);

id_flux = input("1. Burgers 2. Traffic \n");
id_flux = max(mod(floor(id_flux),3),1);
if id_flux == 1 %burgers
    id = input("1. Godunov 2. Lax-F 3. Upwind correct \n");
    id = max(mod(floor(id),4),1);
    uL = 0.5;
    uR = 0.2;
    fb = @(u) 0.5*u.*u;
    csi = @(u) (1/3)*u.^3;
    dt = dx/max(abs(uL),abs(uR));

else %traffic
    id = input("1. Lax-F 2. Upwind correct \n");
    id = max(mod(floor(id),3),1);
    uL = 0.2;
    uR = 0.5;
    fb = @(u) u.*(1-u);
    csi = @(u) 0.5*u.*u - (2/3)*u.^3;
    dt = dx;
end

u0 = @(x) uL*(x < 0) + uR*(x >= 0);
mu = @(u) 0.5*u.*u;

U = u0(x);
Nt = floor(T/dt);
lambda = dt/dx;

I = 1:Nx;
L = I - 1;
R = I + 1;
L(1) = 1;
R(end) = Nx - 1;

n = 0;

Us = zeros(1,Nx);
M = zeros(1,Nx);

p = plot(x,U,x,M);
xlim([-1,1]);
ylim([-0.5,1]);
legend("u","mu")
title("T= "+0);

if id_flux == 1

    switch id
        case 1
            %godunov
            while n < Nt
                n = n+1;
                Us(I) = U(I) - lambda*(F_god(U(I),U(R),fb) - F_god(U(L),U(I),fb));
                M(I) = mu(Us(I)) - mu(U(I)) + lambda*(F_god(U(I),U(R),csi) - F_god(U(L),U(I),csi));
                U = Us;
                p(1).YData = U;
                p(2).YData = 5*M;
                title("T= "+n*dt);
                drawnow;
                pause(0.05)
            end
        case 2
            %lax-f
            while n < Nt
                n = n+1;
                Us(I) = U(I) - lambda*(F_laxf(U(I),U(R),dx,dt,fb) - F_laxf(U(L),U(I),dx,dt,fb));
                M(I) = mu(Us(I)) - mu(U(I)) + lambda*(F_laxf(U(I),U(R),dx,dt,csi) - F_laxf(U(L),U(I),dx,dt,csi));
                U = Us;
                p(1).YData = U;
                p(2).YData = 2*M;
                title("T= "+n*dt);
                drawnow;
                pause(0.05)
            end
        case 3
                %upwind correct convex
                while n < Nt
                    n = n+1;
                    Us(I) = U(I) - lambda*(F_up_convex(U(I),U(R),fb) - F_up_convex(U(L),U(I),fb));
                    M(I) = mu(Us(I)) - mu(U(I)) + lambda*(F_up_convex(U(I),U(R),csi) - F_up_convex(U(L),U(I),csi));
                    U = Us;
                    p(1).YData = U;
                    p(2).YData = 5*M;
                    title("T= "+n*dt);
                    drawnow;
                    pause(0.05)
                end
    end
else 
    switch id
        case 1
            %lax-f
            while n < Nt
                n = n+1;
                Us(I) = U(I) - lambda*(F_laxf(U(I),U(R),dx,dt,fb) - F_laxf(U(L),U(I),dx,dt,fb));
                M(I) = mu(Us(I)) - mu(U(I)) + lambda*(F_laxf(U(I),U(R),dx,dt,csi) - F_laxf(U(L),U(I),dx,dt,csi));
                U = Us;
                p(1).YData = U;
                p(2).YData = 5*M;
                title("T= "+n*dt);
                drawnow;
                pause(0.05)
            end
        case 2
                %upwind correct concave
                while n < Nt
                    n = n+1;
                    Us(I) = U(I) - lambda*(F_up_concave(U(I),U(R),fb) - F_up_concave(U(L),U(I),fb));
                    M(I) = mu(Us(I)) - mu(U(I)) + lambda*(F_up_concave(U(I),U(R),csi) - F_up_concave(U(L),U(I),csi));
                    U = Us;
                    p(1).YData = U;
                    p(2).YData = 5*M;
                    title("T= "+n*dt);
                    drawnow;
                    pause(0.05)
                end
    end

end

function flux = F_laxf(U,V,dx,dt,f) %flux lax-f
    flux = 0.5*(f(U) + f(V)) + 0.5*(dt/dx)*(U - V);
end

function flux = F_god(U1,U2,f) %flux godunov
    flux = max(f(U1),f(U2)).*(U1 >= U2) + f(U2).*(U1 <= U2 & U2 < 0) + f(U1).*(U1 < U2 & U1 > 0);  
end

function flux = F_up_convex(U1,U2,f)
    flux = f(max(U1,0)) + f(min(U2,0));
end

function flux = F_up_concave(U1,U2,f)
    flux = f(min(U1,0.5)) + f(max(U2,0.5)) - 1/4;
end