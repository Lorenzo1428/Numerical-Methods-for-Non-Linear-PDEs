clc
clear
close all

load("InfiniteHorizonData.mat");

theta = linspace(0,2*pi,150);
xp = xc + radius*cos(theta);
yp = yc + radius*sin(theta);
xp1 = xc1 + (radius + 0.1)*cos(theta);
yp1 = yc1 + (radius + 0.1)*sin(theta);
f1 = figure;
contourf(X,Y,V,20); 
hold on;
fill(xp,yp,'r', 'FaceAlpha', 0.4, 'EdgeColor', 'r'); 
legend("target 1")
hold on;
fill(xp1,yp1,'r', 'FaceAlpha', 0.4, 'EdgeColor', 'r'); 
legend("target 2")
hold on;

axis square
axis xy
[x_click,y_click] = ginput(1);
[~, idx] = min(abs(x - x_click));
[~, idy] = min(abs(x - y_click));
xn = [x(idx) x(idy)];
plot(xn(1),xn(2),'bo',"Color",'r','LineWidth',3)
drawnow
hold on
disp("Punto iniziale = [" + xn(1) + " " + xn(2) + "]")

p = plot(xn(1),xn(2),"k.-",Color='r',LineWidth=2.5,MarkerSize=15);
title("Optimal Control");
legend("V function","Target","Trajectory")
xlim([-1,1]);
ylim([-1,1]);

it1max = 500;
it1 = 0;
arg = 1e3*ones(Na,1);
while (xn(1) - xc)^2 + (xn(2) - yc)^2 >= radius^2 && (xn(1) - xc1)^2 + (xn(2) - yc1)^2 >= (radius+0.1)^2 && it1 < it1max
    for m = 1:Na
        xn1 = xn + h*f(xn(1),xn(2),A(m,:));
        if xn1(1) >= 1 || xn1(1) <= -1 || xn1(2) >= 1 || xn1(2) <= -1
            continue;
        else
            [xi,yi] = findcell(xn1,dx);
            Ci = round((xi+1)/dx) +1;
            Cj = round((yi+1)/dx) +1;
            Vc = [V(Ci,Cj) V(Ci+1,Cj) V(Ci,Cj+1) V(Ci+1,Cj+1)];
            Vm = interp((xn1(1) - xi)/dx,(xn1(2) - yi)/dx,Vc);
            arg(m) = (1/(1-h*lambda))*(Vm + h*l(xn(1),xn(2),A(m,:)));
        end   
    end
    [~, idx] = min(arg);
    xn = xn + h*f(xn(1), xn(2), A(idx,:));
    arg = 1e3*ones(Na,1);
    p.XData = [p.XData xn(1)];
    p.YData = [p.YData xn(2)];
    drawnow;
    pause(0.05);
    it1 = it1+1; 
end

function [xi,yi] = findcell(x,dx)
        k = (x(1)+1)/dx;
        l = (x(2)+1)/dx;
        xi = floor(k)*dx -1;
        yi = floor(l)*dx -1;
end