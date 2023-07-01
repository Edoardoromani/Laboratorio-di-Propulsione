% CALCOLO SPESSORE SERBATOIO N2O IN COMPOSITO + LINER ALLUMINIO
% Thrust2 project - Laboratorio di propulsione
% federicobortolato

clear all
close all
% clc

%% Tank Characteristics

p = 6e6;        % Pressione interna serbatoio prevista (Pa)
d = 150;        % Diametro (esterno) serbatoio (mm)
d = d/1000;     % Diametro (esterno) serbatoio (m)
carbon.vol = 0.03;   % Volume serbatoio [m^3]

%% Carbon Fiber Characteristics

carbon.sigma = 3000; 
% Tensione massima ply UD direzione fibre (MPa)
% tens.max di una fibra media scarsa
carbon.FoS = 3;         % Factor of Safety 4 fibers
liner.thick = 1e-3;     % [m] Liner thickness 
liner.rho = 2700;       % [kg/m^3] Liner density
carbon.volume_fraction = 0.5;   % Coefficiente rapporto volumetrico ...
                                % fibre/totale conservativo (si trova 
                                % 0.54 per il filament winding)

carbon.rho_fibers = 1800;   % Densità fibre carbonio [Kg/m^3]
carbon.rho_matrix = 1200;   % Densità matrice epossidica [Kg/m^3]
carbon.rho = carbon.rho_fibers*carbon.volume_fraction + ...
    carbon.rho_matrix*(1-carbon.volume_fraction); 
                            % Densità totale composito [Kg/m^3]

%%  Helical Plies

alpha = 30; % Angolo delle plies oblique (elicoidali) (°)
talpha = p*(d/2-liner.thick)/(2*carbon.sigma*(cosd(alpha))^2+p*carbon.FoS);
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

%% LINER WEIGHT COMPUTATIONS

fprintf("VESSEL LINER COMPUTATIONS\n")
fprintf("-------------------------\n")

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

%% COMPARISON WITH FULL ALUMINUM VESSEL

% export metal vessel design
fprintf("METAL VESSEL COMPUTATIONS\n")
fprintf("-------------------------\n")

%% Tank data

d = d/1000; % [m] Tank diameter
E = 0.5; % Efficienza ASME: TABLE UW-12 (2)(c) -> single-welded butt joint with backing strip, no radiography

%% Thickness of the tank calculated with the small thickness approximation

metal.t1 = d/2*p*metal.FoS/(metal.Ys+p*metal.FoS);                  % [m] Classic

metal.t2 = d/2*p*metal.FoS/(metal.Ys*E+(metal.FoS-0.6)*p);          % [m] Circumferential Stress (Longitudinal Joints)

metal.t3 = d/2*p*metal.FoS/(2*metal.Ys*E+(0.4+metal.FoS)*p);        % [m] Longitudinal Stress (Circumferential Joints) 

metal.ts = d/2*p*metal.FoS/(2*metal.Ys*E+(metal.FoS-0.2)*p);        % [m] Spherical Shells - Hemispherical Heads

metal.tf = d/2*sqrt(metal.C*p)/(sqrt(metal.Ys*E)+sqrt(metal.C*p));  % [m] Flat Heads

%% Composition

%   user input ----------
metal.tc = metal.t1;
metal.th = metal.ts;
% -----------------------

switch metal.tc
    case metal.t1 
        disp("Classic theory")
    case metal.t2
        disp("Circ Stress")
    case metal.t3
        disp("Long Stress")
end
switch metal.th 
    case metal.ts
        disp("Hemi Heads")
    case metal.tf
        disp("Flat Heads")
end

%%  Volume + mass tank calculations

switch metal.th 
    case metal.ts 
        metal.head_vol = 2/3*pi*(d/2-metal.th)^3;                       % [m^3]     Heads Void Volume
        metal.cyl_vol = metal.vol-2*metal.head_vol;                     % [m^3]     Cylinder Void Volume
        metal.cyl_len = metal.cyl_vol/((d/2-metal.tc)^2*pi);            % [m]       Cylinder length
        metal.cyl_intv = ((d/2)^(2)-(d/2-metal.tc)^2)*pi*metal.cyl_len; % [m^(3)]   internal volume of cylinder
        metal.head_intv = 4/3*pi*((d/2)^(3)-(d/2-metal.th)^3);          % [m^(3)]   internal volume of heads
        metal.mass = (metal.cyl_intv+metal.head_intv)*metal.rho;        % [kg]      tank mass
        metal.length = metal.cyl_len + d;                               % [m]       tank tot length
    case metal.tf
        metal.cyl_len = metal.vol/((d/2-metal.tc)^2*pi);                % [m]       Cylinder length
        metal.cyl_intv = ((d/2)^(2)-(d/2-metal.tc)^2)*pi*metal.cyl_len; % [m^(3)]   internal volume of cylinder
        metal.head_intv = d^2*pi*metal.tf;                              % [m^(3)]   internal volume of heads
        metal.mass = (metal.cyl_intv+metal.head_intv)*metal.rho;        % [kg]      tank mass
        metal.length = metal.cyl_len + 2*metal.th;                      % [m]       tank tot length
end

%%  Printing Results 
fprintf("Spessore  cyl: %2.4f mm\n",metal.tc*1000)
fprintf("Spessore head: %2.4f mm\n",metal.th*1000)
fprintf("Peso      Tot: %2.4f Kg\n",metal.mass)
fprintf("-------------------------\n\n")
metal
