clc
clear 
close all

%risolve iconale stazionaria
h = 0.01;
dx = 0.01;
da= 0.01;
x = -1:dx:1;
a = -1:da:1;
L = @(a) 1;

U0 = x.^2;
U = U0;
U1(:,1) = U0;

tol = 1e-12;
itmax = 1e3;
it = 0;
err = 1;
while err > tol && it <= itmax
    U_new = inf(size(U));
    for i = 1:length(a)
        xs = x - a(i)*h;
        bound_index = xs >= x(1) & xs <= x(end);
        xs_bound = xs(bound_index);
        pos = (xs_bound-x(1))/dx +1;
        k = floor(pos);
        k = max(1, min(length(x)-1, k));

        w = pos - k;
        Us = (1-w).*U(k) + w.*U(k+1);
        Ui = h*L(a(i)) + Us;
        U_new(bound_index) = min(U_new(bound_index),Ui);
    end 
    mask_inf = isinf(U_new);
    if any(mask_inf)
        U_new(mask_inf) = U(mask_inf); 
    end

    U_new(1) = 0;
    U_new(end) = 0;

    it = it + 1;
    err = norm(U - U_new,inf);

    U = U_new;
    U1(:,it+1) = U;
end

f = @(x,t) 1 - abs(x);

err_sol = norm(U1(:,end)' - f(x,1),inf);

p = plot(x,U1(:,1),x,f(x,1),"--");
legend("U approx","U esatta");
xlim([x(1) x(end)]);
ylim([0 1]);
for n = 2:size(U1,2)
    p(1).YData = U1(:,n);
    drawnow
    pause(0.1)
end