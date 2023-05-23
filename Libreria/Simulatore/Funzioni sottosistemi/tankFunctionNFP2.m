%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  _____  __  __  ____  _   _  ____  _____                %
%                 |_   _||  ||  || __ \| | | |/ ___||_   _|               %
%                   | |  |  __  ||    /| |_| |\___ \  | |                 %
%                   |_|  |__||__||_|\_\ \___/ |____/  |_|                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Project: Propulsion Simulator

Property of THRUST, unauthorized distribution is not allowed

Description:
  This code implements the mathematical model for the emptying of a tank
  full of N2O from underneath. This indeed implements both the liquid
  discharge and the gas one (tailoff). It uses NFP to model the N2O

Changelog:
  > version: 1.4 - 05/02/2022 - Alessandro Rampazzo
    - decreased abs tolerance and relative tolerance to 1e-3 and 5e-6
    respectively

  > version: 1.3 - 29/01/2022 - Alessandro Rampazzo
    - removed patching of rho_out since another method to solve the problem 
    of discontinuities has been found (see propulsionSimulator changelog)

  > version: 1.2 - 27/01/2022 - Alessandro Rampazzo
    - added rho_out patching to remove discontinuity that caused problems
    with integration

  > version: 1.1 - 22/01/2022 - Alessandro Rampazzo
    - tranformed problem in ds and dm variational problem

  > version: 1.0 - 04/12/2022 - Alessandro Rampazzo
    - added heading, description, and changelog
    - refactoring
    - polishing some calculations
    - joining this function with the old Modello_Svuotamento_Bifase

  > version: 0.1 - 15/03/2022 - Luca Brovedani, Alberto Scomparin
    - created
%}

function [m_out,p_out,rho_out,T_out,x_out] = tankFunctionNFP(T0,m0_ox,V_tank,k)
%% Calcolo condizioni iniziali



if NFP("N2O",'psat_t',T0)*1e+5 > 70e+5
    warning('Pressione massima superiore ai 70 bar!')
end

rho_v = NFP("N2O",'rv_t',T0); %[kg/m^3]
rho_l = NFP("N2O",'rl_t',T0); %[kg/m^3]

x0 = (V_tank*rho_l*rho_v/m0_ox - rho_v)/(rho_l-rho_v); % Titolo di vapore iniziale

if x0 > 1 
    error('L''ossidante NON è in condizioni di saturazione!')
end

if x0 < 0
    error("Non c'è abbastanza volume per alloggiare la massa richiesta alla temperatura data")
end

%% Simulazione svuotamento ossidante liquido

mSpan = linspace(m0_ox,0,2);

% U0 = NFP("N2O","u_tx",T0,x0)*1000*m0_ox;
% 
% 
% n = 2;
% 
% S1 = solveOde(@(mm,Y) liquidOut(mm,Y,V_tank,T0),"RK45_A",mSpan,U0,...
%     "AbsTol",Inf,"RelTol",1e-5,"breakEvent",@(t,Y) noLiquid(t,Y,V_tank,T0));

% m_nominal = S1.t;
% U_nominal = S1.y;
% 
% T_nominal = zeros(size(U_nominal));
% p_nominal = zeros(size(U_nominal));
% x_nominal = zeros(size(U_nominal));
% 
% for i = 1:length(m_nominal)
%     T_nominal(i) = fzero(@(t) mUConditions_tx(U_nominal(i),m_nominal(i),t,V_tank),[260,T0+1e-5]);
%     [~,x_nominal(i)] = mUConditions_tx(U_nominal(i),m_nominal(i),T_nominal(i),V_tank);
%     p_nominal(i) = NFP("N2O",'psat_t',T_nominal(i));
% end

[m_nominal,y] = solveOde(@(mm,Y) liquidOut(mm,Y,V_tank),"RK45_A",mSpan,[x0,T0],...
    "AbsTol",1e-3,"RelTol",5e-6,"breakEvent",@(t,Y) noLiquid(t,Y));

x_nominal = y(:,1);
T_nominal = y(:,2);
p_nominal = NFP("N2O","psat_t",T_nominal)*1e5;

%% Patching results
pp = polyfit(x_nominal(end-1:end),m_nominal(end-1:end),1);
m_x1 = polyval(pp,1);
pp = polyfit(x_nominal(end-1:end),T_nominal(end-1:end),1);
T_x1 = polyval(pp,1);
pp = polyfit(x_nominal(end-1:end),p_nominal(end-1:end),1);
p_x1 = polyval(pp,1);

x_nominal = [x_nominal;1];
m_nominal = [m_nominal;m_x1];
p_nominal = [p_nominal;p_x1];
T_nominal = [T_nominal;T_x1];

%% Tailoff

% svuotamento adiabatico del serbatoio
% m_tailoff = linspace(m_nominal(end),m_nominal(end)*0.95);
% T_tailoff = T_nominal(end) * (m_tailoff/m_tailoff(1)).^(k-1);
% p_tailoff = p_nominal(end) * (m_tailoff/m_tailoff(1)).^(k);
% x_tailoff = ones(size(m_tailoff));

mSpan = linspace(m_nominal(end),0,2);

[m_tailoff,y2] = solveOde(@(mm,Y) gasOut(mm,Y,V_tank),"RK45_A",mSpan,[p_nominal(end)/1e5,T_nominal(end)],...
    "AbsTol",inf,"RelTol",1e-4,"breakEvent",@(t,Y) endTailoff(t,Y));

p_tailoff = y2(:,1);
T_tailoff = y2(:,2);
x_tailoff = ones(size(T_tailoff));

% figure
% hold on
% plot(S2.t,S2.y(:,1))
% plot(S2.t,S2.y(:,2))
%% Unione

T_cutoff = 183; % [K] fusion point of N2O (N2O solidifies at temperatures that are under this one)
I = find(T_tailoff > T_cutoff);

T_out = [T_nominal(:)',T_tailoff(I(2:end))'];
m_out = [m_nominal(:)',m_tailoff(I(2:end))'];
p_out = [p_nominal(:)',p_tailoff(I(2:end))']*1e5;
x_out = [x_nominal(:)',x_tailoff(I(2:end))'];

rho_out = zeros(size(x_out));

for i = 1:length(x_out)
    if x_out(i) == 1
        rho_out(i) = NFP("N2O",'rv_t',T_out(i));
    else
        rho_out(i) = NFP("N2O",'rl_t',T_out(i));
    end
end


%% Nested functions
function dY = liquidOut(m,Y,V_tank)

    x = Y(1);
    T = Y(2);

    if x >= 1
        dY = zeros(size(Y));
        return
    end

    dT = 1e-6;

    dv_dx = NFP("N2O","vv_t",T) - NFP("N2O","vl_t",T);
    dv_dT = (NFP("N2O","v_tx",T+dT,x) - NFP("N2O","v_tx",T-dT,x))/(2*dT);

    dm_dv = -m^2/V_tank;

    dm_dx = dm_dv * dv_dx;
    dm_dT = dm_dv * dv_dT;

    ds_dx = NFP("N2O","sv_t",T) - NFP("N2O","sl_t",T);
    ds_dT = (NFP("N2O","s_tx",T+dT,x) - NFP("N2O","s_tx",T-dT,x))/(2*dT);

    ds_dT = ds_dT*1000; % [kJ/kgK] -> [J/kgK]
    ds_dx = ds_dx*1000; % [kJ/kgK] -> [J/kgK]

    Jinv = 1/(dm_dx*ds_dT - ds_dx*dm_dT)* ...
        [ds_dT, -dm_dT; ...
        -ds_dx, dm_dx];

    ds = (NFP("N2O","sl_t",T) - NFP("N2O","s_tx",T,x))/m*1000;
    
    dY = Jinv * [1;ds];
    
    dY = dY(:).';
end

function [position,isterminal,direction] = noLiquid(m,Y)
    %     fprintf("%1.5f",mm)
    position = Y(1) < 1-1e-6; % continue? 1 (>0) = true, 0 = false
    isterminal = 1;  % Halt integration
    direction = 0;   % The zero can be approached from either direction
end

function dY = gasOut(m,Y,V_tank)

    p = Y(1);
    T = Y(2);

    if T <= 185
        dY = zeros(size(Y));
        return
    end

    dT = 1e-6;
    dp = 1e-6;

    dv_dp = (NFP("N2O","v_tp",T,p-dp) - NFP("N2O","v_tp",T,p-dp*2))/dp;
    dv_dT = (NFP("N2O","v_tp",T+dT*2,p) - NFP("N2O","v_tp",T+dT,p))/dT;

    dm_dv = -m^2/V_tank;

    dm_dp = dm_dv * dv_dp;
    dm_dT = dm_dv * dv_dT;

    ds_dp = (NFP("N2O","h_tp",T,p-dp) - NFP("N2O","h_tp",T,p-dp*2))/dp;
    ds_dT = (NFP("N2O","h_tp",T+dT*2,p) - NFP("N2O","h_tp",T+dT,p))/dT;

    ds_dT = ds_dT*1000; % [kJ/kgK] -> [J/kgK]
    ds_dp = ds_dp*1000; % [kJ/kgK] -> [J/kgK]

    s = NFP("N2O","h_tp",T+dT,p);

    dS_dT = ds_dT*m + dm_dT*s;
    dS_dp = ds_dp*m + dm_dp*s;

    Jinv = 1/(dm_dp*dS_dT - dS_dp*dm_dT)* ...
        [dS_dT, -dm_dT; ...
        -dS_dp, dm_dp];

    dS = s;
    
    dY = Jinv * [1;dS];
    
    dY = dY(:).';
end

function [position,isterminal,direction] = endTailoff(m,Y)
    fprintf("m=%1.5f | p=%1.5f | T=%1.12f | s=%1.12f\n",m,Y(1),Y(2),NFP("N2O","s_pt",Y(1),Y(2)))
%     fprintf("%1.10f %1.10f\n",Y(2),NFP("N2O","s_tp",Y(2)+1e-2,Y(1)))
    position = Y(2) > 185+1e-6; % continue? 1 (>0) = true, 0 = false
    isterminal = 1;  % Halt integration
    direction = 0;   % The zero can be approached from either direction
end

end