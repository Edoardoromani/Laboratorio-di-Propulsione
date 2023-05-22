% estimate garin behaviour from a and n factors

clc
clear
close all

%% Data
% oxidizer fuel mass ratio
O_F_vec = [7.5,5.5];
% total mass flow rate
mDot_prop_vec = [0.63,0.34];    % [kg/s]
% denisty of the fuel
rho_fuel = 1107;    % [kg/m^3]
r_ext = 0.04;       % [m]
r_int = 0.02;     % [m]

% case 1
a(1) = 0.6e-3;      % [mm/s]
n(1) = 0.5;
m(1) = 0;

% case 2
a(2) = 0.029e-3;    % [mm/s]  https://www.mdpi.com/2226-4310/6/8/89/pdf
n(2) = 0.697;
m(2) = 0.398;

% case 3
a(3) = 0.071e-3;    % [mm/s]  https://www.mdpi.com/2226-4310/6/8/89/pdf
n(3) = 0.795;
m(3) = 0;

% case 4
a(4) = 0.227e-3;    % [mm/s]  https://www.mdpi.com/2226-4310/6/8/89/pdf
n(4) = 0.593;
m(4) = 0;

% case 5
a(5) = 0.17e-3;     % [mm/s] 
n(5) = 0.5;
m(5) = 0;

% case 6 Fr5560 C32H66
a(6) = 0.155e-3;    % [mm/s]
n(6) = 0.5;
m(6) = 0;

% case 7 SasolWax C50H102
a(7) = 0.132e-3;    % [mm/s]
n(7) = 0.555;
m(7) = 0;


%% analysis
for i = 1:length(a)
    % print case
    fprintf("case %d\n",i)

    % at the beginning of the burning start
    fprintf("start:\n")
    
    % calculate mass flow rates
    O_F = O_F_vec(1);
    mDot_prop = mDot_prop_vec(1);
    mDot_fuel = mDot_prop / (1+O_F);
    mDot_ox = mDot_prop *O_F/ (1+O_F);
    
    % calculate rDot
    r = r_int;
    G = mDot_ox./(pi.*r.^2);
    rDot_start = a(i) .* G.^n(i) .* (2.*(r_int + r_ext)).^m(i);
    
    % calculate the length necessary to match the required fuel mass flow
    % rate
    L = fzero(@(L) rho_fuel .* L .* 2 .* pi .* r .* rDot_start - mDot_fuel,[0,2]);
    
    % print things
    fprintf("  L: %1.3f cm\n",100*L)
    fprintf("  rDot: %1.3f mm/s\n",rDot_start*1e+3)

    % mean
    fprintf("mean:\n")
    
    % everything as before for mean values and end of burning
    O_F = mean(O_F_vec);
    mDot_prop = mean(mDot_prop_vec);
    mDot_fuel = mDot_prop / (1+O_F);
    mDot_ox = mDot_prop *O_F/ (1+O_F);

    r = (r_int + r_ext)/2;
    G = mDot_ox./(pi.*r.^2);
    rDot_mean = a(i) .* G.^n(i) .* (2.*(r_int + r_ext)).^m(i);

    L = fzero(@(L) rho_fuel .* L .* 2 .* pi .* r .* rDot_mean - mDot_fuel,[0,2]);

    fprintf("  L: %1.3f cm\n",100*L)
    fprintf("  rDot: %1.3f mm/s\n",rDot_mean*1e+3)

    % end
    fprintf("end:\n")

    O_F = O_F_vec(2);
    mDot_prop = mDot_prop_vec(2);
    mDot_fuel = mDot_prop / (1+O_F);
    mDot_ox = mDot_prop *O_F/ (1+O_F);

    r = r_ext;
    G = mDot_ox./(pi.*r.^2);
    rDot_end = a(i) .* G.^n(i) .* (2.*(r_int + r_ext)).^m(i);

    L = fzero(@(L) rho_fuel .* L .* 2 .* pi .* r .* rDot_end - mDot_fuel,[0,2]);

    fprintf("  L: %1.3f cm\n",100*L)
    fprintf("  rDot: %1.3f mm/s\n",rDot_end*1e+3)
    
    % calculate burning time as the time took by the mean regression rate
    % to consume all the fuel
    fprintf("\n>> tb: %1.3f s\n",(r_ext - r_int)/rDot_mean)

    fprintf("\n")
end


