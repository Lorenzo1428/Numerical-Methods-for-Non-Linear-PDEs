clc
clear 
close all

%problema di controllo a orizzonte finito in cui il target si muove con velocita' costante in
%una direzione

gamma = 0; %peso del controllo

Tf = 1;
dt = 0.05;
dx = 0.025;
epsi = 1e-6;
tol = 1e-10;
f = @(x,y,a) [-1,1] + 0*[ y/(sqrt(x^2+y^2) + epsi), -x/(sqrt(x^2+y^2) + epsi) ]  + a;
rho = linspace(0,1.5,4);
rho = rho(2:end);
theta = linspace(0,2*pi,16);
[R,T] = ndgrid(rho,theta);
A = [0, 0 ; R(:).*cos(T(:)) , R(:).*sin(T(:))];
x = -1:dx:1;
N = length(x);
Na = size(A,1);
[X,Y] = ndgrid(x);

%scelta del centro del target
figure;
tile_size = 0.25; 
Checker = mod(floor(X/tile_size) + floor(Y/tile_size), 2);
h = pcolor(X, Y, Checker);
set(h, 'EdgeColor', 'none');
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
pause(0.5)
close all
disp("Centro target = [" + xc + " " + yc + "]")
close all

%velocita' del target e raggio
xvel = -0.7;
yvel = 0.3;
radius = 0.1;
%running cost
l = @(s,x,y,a) (gamma*(a(1)^2 + a(2)^2) + 10)*( (x - (xc+xvel*s))^2 + (y - (yc+yvel*s))^2 >= radius^2);
interp = @(t,s,Vc) (1-t)*(1-s)*Vc(1) + t*(1-s)*Vc(2) + (1-t)*s*Vc(3) + s*t*Vc(4);

s = 0:dt:Tf;
Nt = length(s);
s = s(Nt-1:-1:1);

xc_f = xc + xvel*Tf;
yc_f = yc + yvel*Tf;

%final cost
if abs(xc_f) > 1 || abs(yc_f) > 1
    V = ones(N,N) * 1e12;
else
    dist_sq = (X - xc_f).^2 + (Y - yc_f).^2;
    V = (10+0.5*dist_sq)*(dist_sq > radius^2) + 0*(dist_sq <= radius^2);
end

Vm = zeros(N);
Vk = Vm;

Vt = zeros(N,N,Nt);
Vt(:,:,Nt) = V;
for n = 1:Nt-1
    xc_vel = xc+xvel*s(n);
    yc_vel = yc+yvel*s(n);
    if abs(xc_vel) > 1 || abs(yc_vel) > 1
        V = inf*ones(N,N);
        V(:,:,n) = V;
        continue;
    end

    for i = 1:N
        for j = 1:N
            Vm(i,j) = 1e12;

            if (X(i,j) - (xc_vel))^2 + (Y(i,j) - (yc_vel))^2 <= radius^2
                Vm(i,j) = 0;
                continue;
            end

            for m = 1:Na
                xs = [X(i,j), Y(i,j)] + dt*f(X(i,j),Y(i,j),A(m,:));
                if abs(xs(1)) > 1 || abs(xs(2)) > 1 
                    continue;
                else
                    [xi,yi] = findcell(xs,dx);
                    Ci = round((xi+1)/dx) +1;
                    Cj = round((yi+1)/dx) +1;
                    if Ci < 1 || Ci >= N || Cj < 1 || Cj >= N
                        continue;
                    end
                    Vc = [V(Ci,Cj) V(Ci+1,Cj) V(Ci,Cj+1) V(Ci+1,Cj+1)];
                    Vij = (interp((xs(1) - xi)/dx,(xs(2) - yi)/dx,Vc)+ dt*l(s(n),X(i,j),Y(i,j),A(m,:)));
                    Vm(i,j) = min(Vm(i,j),Vij);
                end
            end
            Vk(i,j) = Vm(i,j);
        end
    end
    V = Vk;
    Vt(:,:,Nt-n) = V;
end

function [xi,yi] = findcell(x,dx)
        k = (x(1)+1)/dx;
        l = (x(2)+1)/dx;
        xi = floor(k)*dx -1;
        yi = floor(l)*dx -1;
end

clearvars i j m n
save("MobileTargetData.mat");