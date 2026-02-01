clc
clear 
close all

lambda = 0; %discount factor
gamma = 1; %peso del controllo

h = 0.05;
dx = 0.05;
epsi = 1e-6;
f = @(x,y,a) [0,0] + [ y/(sqrt(x^2+y^2) + epsi), -x/(sqrt(x^2+y^2) + epsi) ]  + a;
rho = linspace(0,1,4);
rho = rho(2:end);
theta = linspace(0,2*pi,16);
[R,T] = ndgrid(rho,theta);
A = [0, 0 ; R(:).*cos(T(:)) , R(:).*sin(T(:))];
x = -1:dx:1;
N = length(x);
Na = size(A,1);
[X,Y] = ndgrid(x);

figure;
tile_size = 0.25; 
Checker = mod(floor(X/tile_size) + floor(Y/tile_size), 2);
pcolor(X, Y, Checker,'EdgeColor', 'none');
shading flat;
hold on;
colormap([0.90 0.92 0.95; 0.95 0.97 1.0]);
axis equal;
grid on;
set(gca, 'Layer', 'top', 'GridColor', [0.6 0.6 0.6], 'GridAlpha', 0.5);
xlim([-1,1])
ylim([-1 1])

[x_click,y_click] = ginput(1);
[~, idx] = min(abs(x- x_click));
[~, idy] = min(abs(x - y_click));
xc = x(idx);
yc =x(idy);
plot(xc,yc,'bo',"Color",'r','LineWidth',3);
drawnow
hold on 

[x_click,y_click] = ginput(1);
[~, idx] = min(abs(x- x_click));
[~, idy] = min(abs(x - y_click));
xc1 = x(idx);
yc1 =x(idy);
plot(xc1,yc1,'bo',"Color",'r','LineWidth',3);
drawnow
hold on 
pause(0.5)
close all
disp("Centro 1 = [" + xc + " " + yc + "]" + newline + "Centro 2 = [" + xc1 + " " + yc1 + "]")

interp = @(t,s,Vc) (1-t)*(1-s)*Vc(1) + t*(1-s)*Vc(2) + (1-t)*s*Vc(3) + s*t*Vc(4);

tol = 1e-8;
itmax = 500;
radius = 0.1;
l = @(x,y,a) (1 + gamma*(a(1)^2 + a(2)^2))*( (x - xc)^2 + (y - yc)^2 >= radius^2 && (x- xc1)^2 + (y - yc1)^2 >= (radius + 0.1)^2);
V = 100*((X-xc).^2 + (Y-yc).^2 >= radius^2 & (X - xc1).^2 + (Y - yc1).^2 >= (radius + 0.1)^2);
Vm = zeros(N);
Vk = Vm;
err = 1;
it = 0;

while err > tol && it < itmax
    for i = 1:N
        for j = 1:N
            Vm(i,j) = inf;

            if (X(i) - xc)^2 + (X(j) - yc)^2 < radius^2 || (X(i)- xc1)^2 + (X(j) - yc1)^2 < (radius + 0.1)^2
                Vm(i,j) = 0;
                continue;
            end

            for m = 1:Na
                xs = [X(i,j), Y(i,j)] + h*f(X(i,j),Y(i,j),A(m,:));
                if xs(1) >= 1 || xs(1) <= -1 || xs(2) >= 1 || xs(2) <= -1
                    continue;
                else
                    [xi,yi] = findcell(xs,dx);
                    Ci = round((xi+1)/dx) +1;
                    Cj = round((yi+1)/dx) +1;
                    Vc = [V(Ci,Cj) V(Ci+1,Cj) V(Ci,Cj+1) V(Ci+1,Cj+1)];
                    Vij = (1/(1-h*lambda))*(interp((xs(1) - xi)/dx,(xs(2) - yi)/dx,Vc)+ h*l(X(i,j),Y(i,j),A(m,:)));
                    Vm(i,j) = min(Vm(i,j),Vij);
                end
            end
            Vk(i,j) = Vm(i,j);
        end
    end
    err = norm(Vk(:)-V(:),inf);
    V = Vk;
    it = it+1;
end
save("InfiniteHorizonData.mat");

function [xi,yi] = findcell(x,dx)
        k = (x(1)+1)/dx;
        l = (x(2)+1)/dx;
        xi = floor(k)*dx -1;
        yi = floor(l)*dx -1;
end