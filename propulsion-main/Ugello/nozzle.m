%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      %
%              THRUST                  %
%              Nozzle                  %
%                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc 
clear
close all

%% DATA
pcc = 1e5*linspace(41,30,10); % bar -> Pa
% mdot = linspace(0.65,0.35,10); % kg/s
mdot = 0.577;
cstar = 1440;
k = 1.2; % from the cc members
tb = 7; % s, burning time

%% THROAT
% FROM PAPER: A SEMIQUANTITATIVE PREDICTION OF THE EROSION OF THE GRAPHITE
% NOZZLE -> no buono
% % throat area regression (graphite) 
% A_init_ex = 0.12; % inch^2, initial throat area  % ex = example from the paper
% A_fin_ex = 0.065; % inch^2, final throat area
% r_init_ex = sqrt(A_init_ex/pi)*2.54/100; % m, initial area radius
% r_fin_ex = sqrt(A_fin_ex/pi)*2.54/100; % m, initial area radius
% r_dot_graphite = abs(r_fin_ex-r_init_ex)/tb; % m/s
% % const_init = pcc(1)*At(1)^(1/(1-n));
% % const_fin = pcc(end)*At(end)^(1/(1-n));
% 
% % considering the previous example and the found regression rate for the
% % graphite, our final throat area is computed as it follows:
% r_graphite = r_dot_graphite*tb; % m, regression of the radious


% FROM PAPER: ROCKET NOZZLE EROSION GNEM -> no buono

% hypothesis 
r_dot = 0.5e-3; % mm/s -> m/s. assumed constant - r_dot increases with p, our p decreases = no prob
r_erosion_vec = linspace(0,r_dot*tb,10);

At_vec = zeros(length(pcc),1); % m^2, throat area
rt_vec = zeros(length(pcc),1); % m, throat radius

for i = 1:length(pcc)
    
    At_vec(i) = mdot*cstar/pcc(i);
    rt_vec(i) = sqrt(At_vec(i)/pi) + r_erosion_vec(i);
end

%% out-conditions (s condition)
hmax = 900; % m, max altitude reachable
h_adapt = hmax*0.7; % m, altitude to adapt the nozzle
[Ts,as,ps,rhos] = atmosisa(hmax); % K, m/s, Pa, kg/m^3
[Ts_adapt,as_adapt,ps_adapt,rhos_adapt] = atmosisa(h_adapt);

% mach and pressure
Mcc = 0.1; % mach in CC -> ASSUMED -> ASK TO CC MEMBERS
Mout_vec = zeros(1,length(pcc));

for i = 1:length(pcc)
    f1 = @(Mout) pcc(i)./ps_adapt - ((1+(k-1)./2.*Mout.^2)./(1+(k-1)./2.*Mcc.^2)).^(k./(k-1));
    Mout_vec(i) = fzero(f1,[1.2,10]);
end

% areas law
eps_vec = 1./Mout_vec.*(2./(k+1).*(1+(k-1)./2.*Mout_vec.^2)).^((k+1)/(2.*(k-1)))
% d_out_vec = rt_vec*2*eps_vec; % outlet diameter
A_out_vec = eps_vec.*At_vec; % outlet area

% vectors
h_vec = linspace(0,hmax);
Ts_vec = zeros(1,length(h_vec)); 
as_vec = zeros(1,length(h_vec));
ps_vec = zeros(1,length(h_vec));
rhos_vec = zeros(1,length(h_vec));

for i = 1:length(h_vec)
    [Ts_vec(i),as_vec(i),ps_vec(i),rhos_vec(i)] = atmosisa(h_vec(i));
end

%% PLOTS
% plot
figure
plot(pcc./ps_adapt,eps_vec, 'Linewidth',1.5)
hold on 
title('p ratio - eps')
xlabel('p ratio')
ylabel('eps')
grid on

% Ts plot
figure
plot(h_vec,Ts_vec, 'Linewidth',1.5)
hold on 
title('Ts - altitude')
ylabel('Ts [K]')
xlabel('altitude [m]')
grid on

% ps plot
figure
plot(h_vec,ps_vec, 'Linewidth',1.5)
hold on 
title('ps - altitude')
ylabel('Pressure [Pa]')
xlabel('altitude [m]')
grid on

% rhos plot
figure
plot(h_vec,rhos_vec, 'Linewidth',1.5)
hold on 
title('\rho - altitude')
ylabel('\rho [kg/m^3]')
xlabel('altitude [m]')
grid on




