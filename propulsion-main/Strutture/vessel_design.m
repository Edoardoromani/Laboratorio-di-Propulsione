% CALCOLO SPESSORE SERBATOIO N2O IN COMPOSITO
% Thrust2 project - Laboratorio di propulsione
% federicobortolato

clear all
close all
clc

%% Carbon Fiber Characteristics

sigma_lim = 3000; 
% Tensione massima ply UD direzione fibre (MPa)
% tens.max di una fibra media scarsa
p = 5.6; % Pressione interna serbatoio prevista (MPa)
f = 3; % Safety Factor
P = p*f; % Pressione di design serbatoio (MPa)
D = 145; % Diametro (interno) serbatoio (mm)

fibers_volume_fraction = 0.5; 
% Coefficiente rapporto volumetrico fibre/totale conservativo (si trova 0.54 per il filament winding)

rho_fibers = 1800; % Densità fibre carbonio [Kg/m^3]
rho_matrix = 1200; % Densità matrice epossidica [Kg/m^3]
len_cyl = 1.9; % Lunghezza serbatoio [m]

%% Aluminum 6082 Characteristics 

metal.Ys = 250e6; % [Pa] Yield strength
metal.rho = 2700; % [kg/m^(3)] density
metal.FoS = 2; % Factor of Safety
metal2.rho = 2700; % LINER

%%  Helical Plies

alpha = 30; % Angolo delle plies oblique (elicoidali) (°)
talpha = P*D/4/sigma_lim/(cosd(alpha))^2;
talpha = talpha/fibers_volume_fraction;
% Spessore plies oblique con angolo alpha (mm)

%%  Hoop Plies

thoop = (P*D*5/12/sigma_lim); 
thoop = thoop/fibers_volume_fraction;
% Spessore plies ortogonali (mm)

%%  Spessore totale

t = thoop + talpha;
% Spessore complessivo plies (mm)

%% VESSEL VOLUME COMPUTATIONS 

area_cyl_mm = ((D+2*t)^2-D^2)/4*pi; % Area di base cilindro [mm^2]
area_cyl = area_cyl_mm/10^6; % Area di base cilindro [m^2]
vol_cylinder = len_cyl*area_cyl; % Volume cilindro [m^3]
vol_sphere = pi/6*((D+2*t)^3-D^3); % Volume sfera [mm^3]
vol_sphere = vol_sphere/10^9; % Volume sfera [m^3]
vol_vessel = vol_sphere + vol_cylinder; % Volume tot [m^3]

%% VESSEL WEIGHT

vessel_weight = vol_vessel*fibers_volume_fraction*rho_fibers + vol_vessel*(1-fibers_volume_fraction)*rho_matrix;

%%  Printing Results

fprintf("FIBER VESSEL COMPUTATIONS\n")
fprintf("-------------------------\n")
% fprintf("Spessore Hoop: %2.4f mm\n",thoop)
% fprintf("Spessore Heli: %2.4f mm\n",talpha)
fprintf("Spessore  Tot: %2.4f mm\n",t)
fprintf("Peso      Tot: %2.4f Kg\n",vessel_weight)
fprintf("-------------------------\n\n")

%% COMPARISON WITH FULL ALUMINUM VESSEL

% export metal vessel design
fprintf("METAL VESSEL COMPUTATIONS\n")
fprintf("-------------------------\n")

%% Tank data

d = D/1000; % [m] Tank diameter
l = len_cyl; % [m] Tank length
p = p*1e6; % [Pa] Tank pressure
E = 0.65; % Efficienza ASME: TABLE UW-12 (2)(c) -> single-welded butt joint with backing strip, no radiography

%% Thickness of the tank calculated with the small thickness approximation

metal.t1 = d/2*p*metal.FoS/metal.Ys; % [m]
metal.t2 = d/2*p*metal.FoS/(metal.Ys*E-0.6*p); % [m] Circumferential Stress (Longitudinal Joints)
metal.t3 = d/2*p*metal.FoS/(2*metal.Ys*E+0.4*p); % [m] Longitudinal Stress (Circumferential Joints)
metal.t4 = d/2*p*metal.FoS/(2*metal.Ys*E-0.2*p); % [m] Spherical Shells - Hemispherical Heads

%% Verification that small thickness is allowed
thicknesses = [metal.t1 metal.t2 metal.t3];
metal.t = max(thicknesses); % smallest thickness needed
%metal.t = metal.t2;

if metal.t == metal.t1 
    disp("classic humble theory")
elseif metal.t == metal.t2
    disp("Circumferential Stress Governing")
elseif metal.t == metal.t3
    disp("Longitudinal Stress Governing")
elseif metal.t == metal.t4
    disp("Spherical Shells theroy")
end

metal.rate = d/metal.t; % rate between diameter and thickness of the pressure vessel

% if the rate is over 10 the approximation in allowed
VM = 1000e6; % [Pa] allocated for the while
if metal.rate > 10
    fprintf ("small thickness verified\n")
else
    while VM > metal.Ys % increase the thickness until ideal stress > yield stress
    metal.t = metal.t+1e-3; 
    fprintf("small thicknes is an invalid approximation")
    ro = d/2 + metal.t; % [mm] external radius
    ri = d/2; % [mm] internal radius
    c = p*(ri^(2))/(ro^(2)-ri^(2))*(1+(ro/ri)^(2)); % [Mpa] circonferential stress
    r = p*(ri^(2))/(ro^(2)-ri^(2))*(1-(ro/ri)^(2)); % [Mpa] circonferential stress
    l = p*(ri^(2))/(ro^(2)-ri^(2)); % [Mpa] longitudinal stress
    VM = metal.FoS*sqrt(l^(2)+r^(2)+c^(2)+l*c+l*r+r*c); % [Mpa] ideal equivalent stress
    end
end

%% Calculation of the mass of the tank

metal.do = d+2*metal.t;
metal.VolBottoms = pi/6*((metal.do)^(3)-(d)^(3)); % [m^(3)] internal volume of the 2 bottoms (a sphere with internal diameter = d)
metal.VolCilinder = ((metal.do)^(2)-(d)^(2))/4*pi*l; % [m^(3)] internal volume of cylinder
metal.Vol = metal.VolBottoms+metal.VolCilinder;
metal.mass = metal.rho* metal.Vol; % [kg] pressure vessel mass

%%  Printing Results 

fprintf("Spessore  Tot: %2.4f mm\n",metal.t*1000)
fprintf("Peso      Tot: %2.4f Kg\n",metal.mass)
fprintf("-------------------------\n\n")

%% LINER WEIGHT COMPUTATIONS

fprintf("VESSEL LINER COMPUTATIONS\n")
fprintf("-------------------------\n")

%% Liner data

d = (D-2)/1000; % [m] Tank diameter
l = len_cyl; % [m] Tank length
metal2.t = 1e-3; % Thickness [mm]

%% Calculation of the mass of the tank

metal2.do = d+2*metal2.t;
metal2.VolBottoms = pi/6*((metal2.do)^(3)-(d)^(3)); % [m^(3)] internal volume of the 2 bottoms (a sphere with internal diameter = d)
metal2.VolCilinder = ((metal2.do)^(2)-(d)^(2))/4*pi*l; % [m^(3)] internal volume of cylinder
metal2.Vol = metal2.VolBottoms+metal2.VolCilinder;
metal2.mass = metal2.rho* metal2.Vol; % [kg] pressure vessel mass

%%  Printing Results 

fprintf("Spessore  Tot: %2.4f mm\n",metal2.t*1000)
fprintf("Peso      Tot: %2.4f Kg\n",metal2.mass)
fprintf("-------------------------\n\n")

%%  FINAL VESSEL VALUES
fprintf("VESSEL FINAL DATA (fiber+liner)\n")
fprintf("-------------------------\n")
fprintf("Diam. interno: %2.2f mm\n",d*1000)
fprintf("Lunghezza:     %2.4f mm\n",l)
fprintf("Spessore  Tot: %2.4f mm\n",t+metal2.t*1000)
fprintf("Peso      Tot: %2.4f Kg\n",vessel_weight + metal2.mass)
fprintf("-------------------------\n\n")