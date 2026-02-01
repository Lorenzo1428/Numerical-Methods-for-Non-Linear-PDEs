clc
clear 
close all

flux = input("Flusso 1. Burger's 2. Traffic \n");
flux = max(mod(floor(flux), 3), 1);
id = input("1.Upwind 2.Lax-f 3.Mac cormack 4.Lax-w 5. Upwind correct convex 6. Upwind correct concave ... " + ...
    "7.Godunov convex 8. Godunov concave \n" );
id = max(mod(floor(id),9),1);

T = 1;
dx = 0.001;
x = -1:dx:1;
Nx = length(x);

uL = 0;
uR = 0.5;

u0 = @(x) uL*(x < 0) + uR*(x >= 0);
if  flux == 1
    fb = @(u) 0.5*u.*u;
    fb_der = @(u) u;
    dt = dx/max(abs(uL),abs(uR));
else
    fb = @(u) u.*(1-u);
    fb_der = @(u) 1-2*u;
    dt = dx;
end

F_convex = @(U1,U2,f) max(f(U1),f(U2)).*(U1 >= U2) + f(U2).*(U1 <= U2 & U2 < 0) + f(U1).*(U1 < U2 & U1 > 0); 
F_concave = @(U1,U2,f) min(f(U1),f(U2)).*(U1 <= U2) + f(U2).*(U1 > U2 & U2 >= 0.5) + f(U1).*(U1 > U2 & U1 <= 0.5) +...
    0.25*(U1 > U2 & U1 > 0.5 & U2 < 0.5);

U = u0(x);
Nt = floor(T/dt);
lambda = dt/dx;

I = 1:Nx;
L = I - 1;
R = I + 1;
L(1) = 1;
R(end) = Nx - 1;

n = 0;

p = plot(x,U);
%xlim([-1,1]);
ylim([-1,2]);
title("T")
legend("soluzione")

switch id
    case 1
            %upwind
            while n < Nt
                n = n+1;
                U(I) = U(I) - lambda*( (fb(U(R)) - fb(U(I))).*(fb_der(U(I)) <= 0 ) + (fb(U(I)) - fb(U(L))).*(fb_der(U(I)) > 0 ) );
                p.YData = U;
                drawnow;
            end
    case 2   
            %lax f
            while n < Nt
                n = n+1;
                U(I) = 0.5*(U(R) + U(L)) - 0.5*lambda*(fb(U(R)) - fb(U(L)));
                p.YData = U;
                drawnow;
            end
    case 3
            %mac cormack
            Us = zeros(1,Nx);
            while n < Nt+1
                n = n+1;
                Us(I) = U(I) - lambda*(fb(U(R)) - fb(U(I)));
                U(I) = 0.5*(U(I) + Us(I)) - 0.5*lambda*(fb(Us(I)) - fb(Us(L)));
                p.YData = U;
                title("T=" + n*dt);
                drawnow;
            end
    case 4
            %lax w
            Up = zeros(1,Nx);
            Um = zeros(1,Nx);
            while n < Nt
                n = n+1;
                Up(I) = 0.5*(U(R) + U(I)) - 0.5*lambda*(fb(U(R)) - fb(U(I)));
                Um(I) = 0.5*(U(L) + U(I)) - 0.5*lambda*(fb(U(I)) - fb(U(L)));
                U(I) = U(I) - lambda*(fb(Up(I)) - fb(Um(I)));
                p.YData = U;
                drawnow;
            end
    case 5
            %upwind correction convex
            while n < Nt
                n = n+1;
                F_up = fb(max(U(I),0)) - fb(min(U(I),0)) + fb(min(U(R),0)) - fb(max(U(L),0));
                U(I) = U(I) - lambda*F_up;
                p.YData = U;
                drawnow;
            end
    case 6 
            %upwind correction concave
            while n < Nt
                n = n+1;
                F_up = fb(min(U(I),0.5)) - fb(max(U(I),0.5)) + fb(max(U(R),0.5)) - fb(min(U(L),0.5));
                U(I) = U(I) - lambda*F_up;
                p.YData = U;
                drawnow;
            end
    case 7
            %gudonov convex
            while n < Nt
                n = n+1;
                U(I) = U(I) - lambda*(F_convex(U(I),U(R),fb) - F_convex(U(L),U(I),fb));
                p.YData = U;
                drawnow;
            end
    case 8
            %godunov concave
            while n < Nt
                n = n+1;
                U(I) = U(I) - lambda*(F_concave(U(I),U(R),fb) - F_concave(U(L),U(I),fb));
                p.YData = U;
                drawnow;
            end
end

