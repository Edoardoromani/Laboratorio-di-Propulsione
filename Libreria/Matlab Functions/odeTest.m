clc
clear
close all
global a
a = 1;
tspan = linspace(0,10,50);
y0 = [a,0];
S = solveOde(@fun,"RK45_A",tspan,y0,"breakevent",@(t,Y) event(t,Y));

y1 = S.y(:,1);
y2 = S.y(:,2);

plot(S.t,y1,'o-')
hold on
plot(S.t,y2,'o-')

plot(S.t,cos(S.t*a)*a,'k--')
plot(S.t,-sin(S.t*a)*a,'k--')
axis tight

legend("y_1","y_2","analytical"+newline+"solutions","location","best")

function dY = fun(t,Y)
    global a
    dY = zeros(size(Y));
    dY(1) = Y(2)*a;
    dY(2) = -Y(1)*a;
%     try
%         if Y(2) > 0.5*a
%             error("Y(2) grater than 0.5")
%         else
%             dY(1) = Y(2)*a;
%             dY(2) = -Y(1)*a;
%         end
%     catch
%         dY(2) = 0;
%     end
end

function cont = event(t,Y)
    global a
    cont = Y(2) <= (0.5 - 1e-6)*a;
end