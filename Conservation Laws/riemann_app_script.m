function [Us,x,dt] = riemann_app_script(T,dx,uL,uR,fb_value,id)

    x = -1:dx:1;
    Nx = length(x);
    
    u0 = @(x) uL*(x < 0) + uR*(x >= 0);
    if strcmp(fb_value,'Burgers')
        fb = @(u) 0.5*u.*u;
        fb_der = @(u) u;
        dt = dx/max(abs(uL),abs(uR));
    else
        fb = @(u) u.*(1-u);
        fb_der = @(u) 1-2*u;
        dt = dx;
    end
    F_convex = @(U1,U2,f) max(f(U1),f(U2)).*(U1 >= U2) + f(U2).*(U1 < U2 & U2 <= 0) + f(U1).*(U1 < U2 & U1 >= 0); 
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
        
    Us(:,1) = U;
    switch id
        case 1
                %upwind
                while n < Nt
                    n = n+1;                  
                    U(I) = U(I) - lambda*( (fb(U(R)) - fb(U(I))).*(fb_der(U(I)) <= 0 ) + (fb(U(I)) - fb(U(L))).*(fb_der(U(I)) > 0 ) );
                    Us(:,n+1) = U;
                end
        case 2   
                %lax f
                while n < Nt
                    n = n+1;
                    U(I) = 0.5*(U(R) + U(L)) - 0.5*lambda*(fb(U(R)) - fb(U(L)));
                    Us(:,n+1) = U;
                end
        case 3
                %macCormack
                U_s = zeros(1,Nx);
                while n < Nt+1
                    n = n+1;
                    U_s(I) = U(I) - lambda*(fb(U(R)) - fb(U(I)));
                    U(I) = 0.5*(U(I) + U_s(I)) - 0.5*lambda*(fb(U_s(I)) - fb(U_s(L)));
                    Us(:,n+1) = U;
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
                    Us(:,n+1) = U;
                end
        case 5
               %upwind correction convex
                while n < Nt
                    n = n+1;
                    F_up = fb(max(U(I),0)) - fb(min(U(I),0)) + fb(min(U(R),0)) - fb(max(U(L),0));
                    U(I) = U(I) - lambda*F_up;
                    Us(:,n+1) = U;
                end
        case 6 
                %upwind correction concave
                while n < Nt
                    n = n+1;
                    F_up = fb(min(U(I),0.5)) - fb(max(U(I),0.5)) + fb(max(U(R),0.5)) - fb(min(U(L),0.5));
                    U(I) = U(I) - lambda*F_up;
                    Us(:,n+1) = U;
                end
        case 7
                %godunov convex
                while n < Nt
                    n = n+1;
                    U(I) = U(I) - lambda*(F_convex(U(I),U(R),fb) - F_convex(U(L),U(I),fb));
                    Us(:,n+1) = U;
                end
        case 8
                %godunov concave
                while n < Nt
                    n = n+1;
                    U(I) = U(I) - lambda*(F_concave(U(I),U(R),fb) - F_concave(U(L),U(I),fb));
                    Us(:,n+1) = U;
                end
    end

end