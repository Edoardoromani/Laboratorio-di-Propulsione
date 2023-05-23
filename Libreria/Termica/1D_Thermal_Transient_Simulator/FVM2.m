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
  This code implements the a 1D Finite Volume Method for thermal analysis.
  It simply solves the following equation:
                  rho*Cp*dT/dt = div(lambda * grad(T))
  Ablation is implemented

Changelog:
  > version: 1.0 - 13/11/2022 - Alessandro Rampazzo
    - created
%}

function dY = FVM2(t,Y,x,mat,BC,enableAbl)
    %% Y = [T1,T2,T3,...,Tn,x_burned]

    % initialize dY
    dY = zeros(size(Y));
    
    % extract x_abl
    x_abl = Y(end);
    
    % check if ablation is active but boundary conditions are not set on
    % fixed temperature
    if enableAbl && lower(BC.sx.type) ~= "temperature" && lower(BC.sx.type) ~= "temp"
        error("ablation is active but left boundary condition is not set on temperature")
    end
    
    % find regression index
    I = find(x < x_abl);
    if isempty(I)
        regIdx = 0;
    else
        regIdx = I(end);
    end
    
    % find regression rate of the cell
    [~,II] = min(abs(x - x_abl));
    if enableAbl
        dY(end) = getRegressionRate(mat(II),t,Y(II));
    else
        dY(end) = 0;
    end

    lambdas = zeros(size(x));
    rhos = lambdas;
    cps = lambdas;

    for i = 1:length(lambdas)
        lambdas(i) = getLambda(mat(i),Y(i));
        rhos(i) = getRho(mat(i),Y(i));
        cps(i) = getCp(mat(i),Y(i));
    end
    

    %% sx boundary
    switch lower(BC.sx.type)
        case {"adiabatic","ad"}
            % correct temperature for ablated material
            Y(1:regIdx) = Y(1);

            % faces coordinates
            xf_dx = (x(regIdx+2)+x(regIdx+1))/2;
            xf_sx = x_abl;

            % cell width
            dx = xf_dx - xf_sx;   

            % thermal resistances
            Ri_dx = abs(x(regIdx+1)-xf_dx)/getLambda(mat(regIdx+1),Y(regIdx+1));
            Rdx_dx = abs(x(regIdx+2)-xf_dx)/getLambda(mat(regIdx+2),Y(regIdx+2));

            % heat
            Q_dx = (Y(regIdx+2) - Y(regIdx+1))/(Ri_dx + Rdx_dx)/dx;

            % rho * Cp
            rhoCp = getRho(mat(regIdx+1),Y(regIdx+1))*getCp(mat(regIdx+1),Y(regIdx+1));

            % dT/dt
            dY(regIdx+1) = Q_dx/rhoCp;

            % correct temperature derivative for ablated material
            dY(1:regIdx) = dY(regIdx+1);

        case {"temperature","temp"}
            Y(1:max(regIdx,1)) = BC.sx.value;
            regIdx = max(regIdx,1);
            xf_sx = max(x_abl,(x(regIdx) + x(regIdx+1))/2);
            xf_dx = (x(regIdx+2)+x(regIdx+1))/2;
            dx = xf_dx - xf_sx;
            Ri_dx = abs(x(regIdx+1)-xf_dx)/getLambda(mat(regIdx+1),Y(regIdx+1));
            Rdx_dx = abs(x(regIdx+2)-xf_dx)/getLambda(mat(regIdx+2),Y(regIdx+2));
            Ri_sx = abs(x(regIdx+1)-xf_sx)/getLambda(mat(regIdx+1),Y(regIdx+1));
            Rsx_sx = abs(x_abl-xf_sx)/getLambda(mat(regIdx),Y(regIdx));
            Q_dx = (Y(regIdx+2) - Y(regIdx+1))/(Ri_dx + Rdx_dx)/dx;
            Q_sx = (Y(regIdx) - Y(regIdx+1))/(Ri_sx + Rsx_sx)/dx;
            % To compensate numerical errors when the cell becomes too
            % small due to ablation
            if enableAbl && Q_sx + Q_dx < 0
                Q_sx = -Q_dx;
            end
            rhoCp = getRho(mat(regIdx+1),Y(regIdx+1))*getCp(mat(regIdx+1),Y(regIdx+1));
            dY(regIdx+1) = 1/rhoCp * (Q_dx + Q_sx);
        case {"heat flux","hf"}
            Y(1:regIdx) = Y(1);
            xf_dx = (x(regIdx+2)+x(regIdx+1))/2;
            xf_sx = x_abl;
            dx = xf_dx - xf_sx;   
            Ri_dx = abs(x(regIdx+1)-xf_dx)/getLambda(mat(regIdx+1),Y(regIdx+1));
            Rdx_dx = abs(x(regIdx+2)-xf_dx)/getLambda(mat(regIdx+2),Y(regIdx+2));
            Q_dx = (Y(regIdx+2) - Y(regIdx+1))/(Ri_dx + Rdx_dx)/dx;
            rhoCp = getRho(mat(regIdx+1),Y(regIdx+1))*getCp(mat(regIdx+1),Y(regIdx+1));
            dY(regIdx+1) = (Q_dx + BC.sx.value(Y(regIdx+1))/dx) / rhoCp;
            dY(1:regIdx) = dY(regIdx+1);
        otherwise
            error("boundary condition not recognized")
    end

        %% right boundary
    switch lower(BC.dx.type)
        case {"adiabatic","ad"}
            xf_sx = (x(end-2)+x(end-1))/2;
            xf_dx = x(end-1);
            dx = xf_dx - xf_sx;
            Ri_sx = abs(x(end-1)-xf_sx)/getLambda(mat(end-1),Y(end-1));
            Rsx_sx = abs(x(end-2)-xf_sx)/getLambda(mat(end-2),Y(end-2));
            Q_sx = (Y(end-2) - Y(end-1))/(Ri_sx + Rsx_sx)/dx;
            rhoCp = getRho(mat(end-1),Y(end-1))*getCp(mat(end-1),Y(end-1));
            dY(end-1) = Q_sx/rhoCp;
        case {"temperature","temp"}
            Y(end-1) = BC.dx.value;
            dY(end-1) = 0;

        case {"heat flux","hf"}
            xf_sx = (x(end-2)+x(end-1))/2;
            xf_dx = x(end-1);
            dx = xf_dx - xf_sx;
            Ri_sx = abs(x(end-1)-xf_sx)/getLambda(mat(end-1),Y(end-1));
            Rsx_sx = abs(x(end-2)-xf_sx)/getLambda(mat(end-2),Y(end-2));
            Q_sx = (Y(end-2) - Y(end-1))/(Ri_sx + Rsx_sx)/dx;
            rhoCp = getRho(mat(end-1),Y(end-1))*getCp(mat(end-1),Y(end-1));
            dY(end-1) = (Q_sx + BC.dx.value(Y(end-1))/dx)/rhoCp;
        otherwise
            error("boundary condition not recognized")
    end
    
    %% middle volumes
    for i = regIdx+2:length(x)-1    
        xf_sx = (x(i-1)+x(i))/2;
        xf_dx = (x(i+1)+x(i))/2;
        dx = xf_dx - xf_sx;
        Ri_dx = abs(x(i)-xf_dx)/getLambda(mat(i),Y(i));
        Rdx_dx = abs(x(i+1)-xf_dx)/getLambda(mat(i+1),Y(i+1));
        Ri_sx = abs(x(i)-xf_sx)/getLambda(mat(i),Y(i));
        Rsx_sx = abs(x(i-1)-xf_sx)/getLambda(mat(i-1),Y(i-1));
        Q_dx = (Y(i+1) - Y(i))/(Ri_dx + Rdx_dx)/dx;
        Q_sx = (Y(i-1) - Y(i))/(Ri_sx + Rsx_sx)/dx;
        % To compensate numerical errors when the cell becomes too
        % small due to ablation
        if enableAbl
            if i == regIdx + 1 && Q_sx + Q_dx < 0
                Q_sx = -Q_dx;
            end
        end
        rhoCp = getRho(mat(i),Y(i))*getCp(mat(i),Y(i));
        dY(i) = 1/rhoCp * (Q_dx + Q_sx);
    end
end
