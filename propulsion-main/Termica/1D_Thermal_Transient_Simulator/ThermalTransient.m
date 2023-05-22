%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  _____  __  __  ____  _   _  ____  _____                %
%                 |_   _||  ||  || __ \| | | |/ ___||_   _|               %
%                   | |  |  __  ||    /| |_| |\___ \  | |                 %
%                   |_|  |__||__||_|\_\ \___/ |____/  |_|                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Project: 1D Thermal Transient Simulator

Property of THRUST, unauthorized distribution is not allowed

Description:
  This code simulates a 1D thermal transient using Finite Volume Method
  (FVM). Ablation, multiple sheets of material and several boundary 
  conditions are implemented.

Changelog:
  > version: 1.0 - 13/11/2022 - Alessandro Rampazzo
    - created
%}

clc
clear
close all

format compact

%% macro parameters
% plot verification using the analytical solution (ablation not supported)
plotVerification = false;

% save plot
saveImage = false;
fileName = "nozzle_throat";

% available analysis types:
% - "minthickness" or "mint": find out the minimum required thickness of
%     a material to ensure a temperature less than T_lim. To make it work
%     choose a longer thickness than the minimum expected
% - "equilibrium" or "eq": find the equilibrium temperature. To make it
%     work you need to propagate more than the expected time required
analysisType = "none";

% limit temperature of AA 6061 and AA 6082: 
% https://www.makeitfrom.com/material-properties/6061-T6-Aluminum
T_lim = 273.15 + 170; % [K]

% colormap
cMap = hsv;

%% mesh generation

% material ID List:
%    material                 ID
% ________________________________
%   Graphite                  0
%   Steel                     1
%   Phenolic resin            2
%   AA 7075 - T6              3
%   Paraffin sasolwax 0907    4
%   HDPE                      5
%   AA 6082 - T6              6
%   AA 6061 - T6              7

% thickness [m]     ID     N nodes
meshParams = {
    5e-3,           2,       50;
    1e-3,           5,       20;
    7e-3,           1,       50};

% for analytical verification with steel
% meshParams = {20e-3,   1,     50};

% initialization temperature
T0 = 293.15; % [K]

% mesh generation
mesh = MeshTT(meshParams,T0);


%% first period propagation

tspan1 = 7; % [s]
CFLFactor1 = 0.5;
% enable ablation (only available for left fixed temperature boundary
% condition)
enableAbl = true;

% boundary condition
% BC1.sx.type = "hf";
% BC1.sx.value = @(T) 6000 * (3000 - T);
BC1.dx.type = "ad";

% for analytical verification with steel
BC1.sx.type = "temp";
BC1.sx.value = 1000;
% BC1.dx.type = "temp";
% BC1.dx.value = 293.15;

mesh = mesh.setBoundaryConditions(BC1);

% propagation
mesh = mesh.propagate(tspan1,CFLFactor1,enableAbl);

fprintf("\nfinished integrating first period, starting second period:\n\n")

%% second period propagation
tspan2 = 10;
CFLFactor2 = 0.5;

if tspan2 > 0

    % boundary conditions
    BC2.sx.type = "ad";
    BC2.dx.type = "ad";

    % propagation
    mesh = mesh.setBoundaryConditions(BC2);
    mesh = mesh.propagate(tspan2,CFLFactor2);
    
end
fprintf("\nfinished integrating second period\n\n")


%% plots

% interpolation
tPlot = linspace(0,max(mesh.t),20);
x_ablPlot = interp1(mesh.t,mesh.x_abl,tPlot);
xPlot = mesh.x;
[xx,tt] = meshgrid(xPlot,tPlot);
[xMesh,tMesh] = meshgrid(mesh.x,mesh.t);
TPlot = interp2(xMesh,tMesh,mesh.T,xx,tt,"linear");

clear xMesh
clear tMesh

figure(1)
hold on

% error initialization for verification
err = 0;
for i = 1:length(tPlot)
    % cropping plots to account for regression rate
    xplot2 = xPlot;
    for j = 1:length(xplot2)
        if xplot2(j) < x_ablPlot(i)
            xplot2(j) = x_ablPlot(i);
        end
    end
    
    % verification
    if plotVerification
        T_hot = 1000;
        p1 = plot(mesh.x*1e3,TPlot(i,:) - 273.15,'b');
        matID = 1; % Steel
        p2 = plot(mesh.x*1e3,Tanalitical(mesh.x,tPlot(i),T_hot,T0,mesh.x(end),matID)-273.15,'r');
        % error calculation
        if i > 5
            err = err + norm(TPlot(i,:) - Tanalitical(xplot2,tPlot(i),T_hot,T0,mesh.x(end),matID));
        end
    else
        % plotting the temperature profiles with different colors
        color = cMap(round(255*tPlot(i)/max(tPlot))+1,:);
        plot(xplot2*1e3,TPlot(i,:)-273.15,'Color',color,"LineWidth",1)
    end
end

% resizing and personalizing plot
ylim([-inf, max(TPlot,[],'all')*1.2])
title("Temperature - " + sprintf("%1.1f",tspan1+tspan2) + "s")
xlabel("x [mm]")
ylabel("T [°C]")
box on

% limit temperature line
yline(T_lim-273.15,'k--',"T = "+(T_lim-273.15)+" °C")

% material lines
xline(0,'k--',getName(mesh.mat(1)),'LabelHorizontalAlignment','right','LabelOrientation','aligned')
for i = 2:length(mesh.mat)
    if mesh.mat(i) ~= mesh.mat(i-1)
        xline((mesh.x(i-1) + mesh.x(i))*1e3/2,'k--',getName(mesh.mat(i)),'LabelHorizontalAlignment','right','LabelOrientation','aligned')
    end
end

% add analysis type
switch lower(analysisType)   
    case {"equilibrium","eq"}
        % find equilibrium temperature and plot it
        T_eq = mesh.T(end,1);
        if max(abs(mesh.T(end,:) - T_eq)./T_eq) < 1e-3
            lim = ylim;
            yticks = [T0-273.15,100:100:lim(2)];
            Ieq = find((yticks - (T_eq-273.15)) < 0);
            Ieq = Ieq(end);
            yticks = [yticks(1:Ieq),round(T_eq-273.15,2),yticks(Ieq+1:end)];
            fprintf("Equilibrium temperature is: %1.3f °C\n",T_eq-273.15)
            set(gca,"YTick",yticks)
        else
            warning("it appears equilibrium has not been reached, try increasing the time span")
        end
    case {"minthickness","mint"}
        % find minumum thickness required to have a gievn temperature
        % behind a wall
        try
            f = @(xx) interp1(mesh.x,mesh.T(end,:),xx,'pchip');
            min_thickness = fzero(@(xx) f(xx) - T_lim,[min(mesh.x),max(mesh.x)]);
            fprintf("minimum thickness to ensure T < %1.1f °C on outer face is: %1.3f mm\n",T_lim-273.15,min_thickness*1000)
        catch ME
            if ME.message == "Function values at the interval endpoints must differ in sign."
                warning("cannot find minimum thickness to ensure T < %1.1f °C on outer face",T_lim-273.15)
            end
        end
end

% if verification is activated add legend
if plotVerification
    fprintf("err : %1.3f\n",err)
    legend([p1,p2],["Numerical simulation","Analitical solution"])
end

% save image
if saveImage
    imwrite(frame2im(getframe(figure(1))),fileName + "_" + sprintf("%1.1f",tspan1+tspan2) + "s.jpg")
end




%% Functions

function T = Tanalitical(x,t,T0,Ta,l,matID)
    % only for fixed temperature
    T = Ta + (T0-Ta)*(1-x/l);
    n = 50;
    k = getLambda(matID,Ta)/(getCp(matID,Ta)*getRho(matID,Ta));
    for i = 1:n
        T = T - 2*(T0-Ta)/(i*pi) * exp(-(i*pi/l)^2*t*k) * sin(i*pi/l*x);
    end
end