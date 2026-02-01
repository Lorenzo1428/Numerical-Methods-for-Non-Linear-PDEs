clc
clear 
close all

T = 1;
dt = 0.1;
dx = 0.001;
Nt = floor(T/dt);
h = 0.01;
a = -1:h:1;

id = input("1. Burgers 2. H(p) = |p| 3. eikonale evolutiva \n");
id = max(mod(floor(id),4),1);

L = @(a) 0.5*abs(a)^2*(id == 1) + 0*(id == 2) + 1*(id == 3);
b = 1*(id == 3) + 2*(id ~= 3);
x = -b:dx:b;
U0 = abs(x);
U = U0;
U1(:,1) = U0;
for  n = 1:Nt    
    U_new = inf(size(U));
    for i = 1:length(a)
        xs = x - a(i)*dt;
        bound_index = xs >= x(1) & xs <= x(end);
        xs_bound = xs(bound_index);
        pos = (xs_bound-x(1))/dx +1;
        k = floor(pos);
        k = max(1, min(length(x)-1, k));

        w = pos - k;
        Us = (1-w).*U(k) + w.*U(k+1);
        Ui = dt*L(a(i)) + Us;
        U_new(bound_index) = min(U_new(bound_index),Ui);
    end 
    mask_inf = isinf(U_new);
    if any(mask_inf)
        U_new(mask_inf) = U(mask_inf); 
    end
    
    if id == 3 %condizione bordo eikonale
        U_new(1) = 0;
        U_new(end) = 0;
    end
    U = U_new;
    U1(:,n+1) = U;
end

switch id
    case 1
        f = @(x,t) (abs(x) - t*0.5).*(abs(x) > t) + (0.5*abs(x).^2/t).*(abs(x) <= t);
    case 2
        f = @(x,t) (abs(x) - t).*(abs(x) > t) + (zeros(1,length(x))).*(abs(x) <= t);
    case 3
        f = @(x,t) 2 - abs(x);
end

err = norm(U1(:,end)' - f(x,T),inf);

p = plot(x,U1(:,1),x,f(x,T),"--");
legend("U approx","U esatta");
title("T = 0");
xlim([-b b]);
ylim([0 b]);
for n = 1:Nt
    p(1).YData = U1(:,n+1);
    title("T = " + n*dt);
    drawnow
    pause(0.1)
end