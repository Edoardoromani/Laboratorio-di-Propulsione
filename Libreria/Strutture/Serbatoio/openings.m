%   RIFERIMENTI IN ASME UG-37
clear
clc
close all

%% start

d = 20.34;   % da cambiare
t = 3.5;    % spessore testa
tr = 2.29;  % spessore foro
tn = (17.32-d/2) % spessore ugello

%% Tank Characteristics

metal.Ys = 288e6;   % [Pa] Yield strength
metal.FoS = 3;      % Factor of Safety

p = 6e6;       % Pressione interna serbatoio prevista (Pa)

%% Thickness of the tank calculated with the small thickness approximation
D = d/1000;
E = 1;

metal.t1 = D/2*p*metal.FoS/(metal.Ys+p*metal.FoS);                  % [m] Classic

metal.t2 = D/2*p*metal.FoS/(metal.Ys*E+(metal.FoS-0.6)*p);          % [m] Circumferential Stress (Longitudinal Joints)

metal.t3 = D/2*p*metal.FoS/(2*metal.Ys*E+(0.4+metal.FoS)*p)        % [m] Longitudinal Stress (Circumferential Joints) 

trn = metal.t3; % spessore ugello minimo

A1_a = d*(t-tr)
A1_b = 2*(t+tn)*(t-tr)

A1 = max(A1_a,A1_b);

% A2 = 5*(tn-trn)*t
A2 = 2*(tn-trn)*9.33

A = d*tr

if A < A1 + A2
    disp('eureka')
elseif A == A1 + A2
    disp("così così")
else
    disp("debacle")
end