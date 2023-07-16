% CALCOLO SPESSORE SERBATOIO N2O IN COMPOSITO + LINER ALLUMINIO
% Thrust2 project - Laboratorio di propulsione
% federicobortolato

clear all
close all


%% Tank Characteristics


tank.vol = 0.03;                % Volume serbatoio [m^3]
tank.FoS = 3;                   % Factor of Safety for fibers vessels

d = 150;                        % Diametro (esterno) serbatoio (mm)
d = d/1000;                     % Diametro (esterno) serbatoio (m)

p = 6e6;                        % Pressione interna serbatoio prevista (Pa)
p = p*tank.FoS;                 % Design burst pressure

%% Carbon Fiber Characteristics

carbon.sigma = 3000e6;             % Tensione massima ply UD direzione fibre (Pa)
                                % tens max di una fibra media scarsa
carbon.volume_fraction = 0.5;   % Coefficiente rapporto volumetrico fibre/totale conservativo (si trova 0.54 per il filament winding)
carbon.rho_fibers = 1800;       % Densità fibre carbonio [Kg/m^3]
carbon.rho_matrix = 1200;       % Densità matrice epossidica [Kg/m^3]
carbon.rho = carbon.rho_fibers*carbon.volume_fraction + ...
    carbon.rho_matrix*(1-carbon.volume_fraction); 
                                % Densità totale composito [Kg/m^3]

%%  Liner + Heads Characteristics

liner.thick = 1e-3;             % [m] Liner thickness
liner.rho = 2700;               % [kg/m^3] Liner density

metal.Ys = 288e6;   % [Pa] Yield strength
metal.rho = 2700;   % [kg/m^(3)] density
metal.C = 0.10;     % Coefficient Flat Heads
metal.E = 0.5;      % Coefficient Efficiency Welded Joint

%%  Calculations start

only_heli = "no";   
                    % "yes" or "no" 
                    % "yes" if intend to use only helical plies
                    %
                    % !!!
                    % non cambia lo spessore totale del carbonio, cambia
                    % solo la ripartizione tra spessore heli e spessore
                    % hoop
                    % !!!

%%  Cylinder Plies calcs

switch only_heli
    case "yes"
        alpha = 54.7; % Angolo delle plies oblique (elicoidali) (°) - necessario se only heli
        talpha = p*(d/2-liner.thick)/(2*carbon.sigma*(cosd(alpha))^2+p);
        carbon.talpha = talpha/carbon.volume_fraction;     % Spessore plies oblique con angolo alpha (mm)
        carbon.thick = carbon.talpha;  % Spessore totale carbon
        carbon.thoop = 0;
    case "no"
        alpha = 0; % Angolo delle plies oblique (elicoidali) (°)
        syms talpha thoop
        eqnalpha = talpha == p*(d/2-liner.thick-thoop)/(2*carbon.sigma*(cosd(alpha))^2+p);
        eqnhoop = thoop == p*(1-0.5*(tand(alpha))^2)*(d/2-liner.thick-talpha)/(carbon.sigma+p*(1-0.5*(tand(alpha))^2));
        [carbon.talpha, carbon.thoop] = solve([eqnalpha,eqnhoop], [talpha, thoop]);  
        carbon.talpha = double(carbon.talpha/carbon.volume_fraction);       % Spessore plies oblique con angolo alpha (mm)
        carbon.thoop = double(carbon.thoop/carbon.volume_fraction);         % Spessore plies ortogonali (mm)
        carbon.thick = carbon.talpha + carbon.thoop;                        % Spessore totale carbon
    otherwise 
        error('variable only_heli must be "yes" or "no"')
end


%%  Heads calcs

metal.ts = d/2*p/(2*metal.Ys*metal.E+(1-0.2)*p);                                    % [m] Hemispherical Heads

metal.tf = (d-2*(carbon.thick+liner.thick))*sqrt(metal.C*p/metal.E/metal.Ys);       % [m] Flat Heads

%% Composition

%   choose head: tf or ts
%   user input ----------
metal.th = metal.ts;
% -----------------------

%%  Volume + mass tank calculations

switch metal.th 
    case metal.ts 
        metal.head_voiv = 4/3*pi*(d/2-metal.th)^3;                          % [m^3]     Heads Void Volume
        metal.head_intv = 4/3*pi*((d/2)^(3)-(d/2-metal.th)^3);              % [m^(3)]   Heads Internal volume (material)

        carbon.cyl_voiv = tank.vol-metal.head_voiv;                         % [m^3]     Cylinder Void Volume
        carbon.cyl_len = carbon.cyl_voiv/...
            ((d/2-carbon.thick-liner.thick)^2*pi);                          % [m]       Cylinder length
        carbon.cyl_intv = ((d/2)^(2)-(d/2-carbon.thick)^2)*...
            pi*carbon.cyl_len;                                              % [m^(3)]   Carbon cylinder Internal volume
        
        liner.cyl_len = carbon.cyl_len;                                     % [m]       Cylinder length
        liner.cyl_intv = ((d/2-carbon.thick)^2-...
            (d/2-carbon.thick-liner.thick)^2)*pi*carbon.cyl_len;            % [m^(3)]   Liner cylinder Internal volume

        tank.mass = liner.cyl_intv*liner.rho + metal.head_intv*metal.rho + ...
            carbon.cyl_intv*carbon.rho;                                     % [kg]      tank mass
        tank.length = carbon.cyl_len + d;                                   % [m]       tank tot length

    case metal.tf

        carbon.cyl_voiv = tank.vol;                                         % [m^3]     Cylinder Void Volume
        carbon.cyl_len = carbon.cyl_voiv/...
            ((d/2-carbon.thick-liner.thick)^2*pi);                          % [m]       Cylinder length
        carbon.cyl_intv = ((d/2)^(2)-(d/2-carbon.thick)^2)*...
            pi*carbon.cyl_len;                                              % [m^(3)]   Carbon cylinder Internal volume

        liner.cyl_len = carbon.cyl_len;                                     % [m]       Cylinder length
        liner.cyl_intv = ((d/2-carbon.thick)^2-...
            (d/2-carbon.thick-liner.thick)^2)*pi*carbon.cyl_len;            % [m^(3)]   Liner cylinder Internal volume

        metal.head_intv = d^2*pi*metal.th;                                  % [m^(3)]   internal volume of heads
        
        tank.mass = liner.cyl_intv*liner.rho + metal.head_intv*metal.rho + ...
            carbon.cyl_intv*carbon.rho;                                     % [kg]      tank mass
        tank.length = carbon.cyl_len + 2*metal.th;                          % [m]       tank tot length

end

%%  Printing Results 
fprintf("--------------------------\n")
fprintf("---CARBON VESSEL DESIGN---\n")
fprintf("--------------------------\n")

switch metal.th 
    case metal.ts
        disp("Hemi Heads")
    case metal.tf
        disp("Flat Heads")
end
fprintf("Lunghezza tot: %2.4f mm\n",tank.length*1000)
fprintf("Spessore carb: %2.4f mm\n",carbon.thick*1000)
fprintf("Spessore alph: %2.4f mm\n",carbon.talpha*1000)
fprintf("Spessore hoop: %2.4f mm\n",carbon.thoop*1000)
fprintf("Spessore head: %2.4f mm\n",metal.th*1000)
fprintf("Peso      Tot: %2.4f Kg\n",tank.mass)
fprintf("--------------------------\n\n")

format compact
format shortg
