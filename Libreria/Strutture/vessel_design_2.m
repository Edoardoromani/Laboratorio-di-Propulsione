% CALCOLO SPESSORE SERBATOIO N2O IN COMPOSITO
% Thrust2 project - Laboratorio di propulsione
% federicobortolato

clear all
close all
clc

%% Tank Characteristics

metal.Ys = 288e6;   % [Pa] Yield strength
metal.rho = 2700;   % [kg/m^(3)] density
metal.FoS = 2;      % Factor of Safety
metal.C = 0.33;     % Coefficient Flat Heads
metal.vol = 0.03;   % Volume serbatoio [m^3]

p = 6e6;       % Pressione interna serbatoio prevista (Pa)
d = 150;        % Diametro (esterno) serbatoio (mm)

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
metal.tc = metal.t3;
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
    case metal.tf
        metal.cyl_len = metal.vol/((d/2-metal.tc)^2*pi);                % [m]       Cylinder length
        metal.cyl_intv = ((d/2)^(2)-(d/2-metal.tc)^2)*pi*metal.cyl_len; % [m^(3)]   internal volume of cylinder
        metal.head_intv = d^2*pi*metal.tf;                              % [m^(3)]   internal volume of heads
        metal.mass = (metal.cyl_intv+metal.head_intv)*metal.rho;        % [kg]      tank mass
end

%%  Printing Results 

fprintf("Peso      Tot: %2.4f Kg\n",metal.mass)
fprintf("-------------------------\n\n")

