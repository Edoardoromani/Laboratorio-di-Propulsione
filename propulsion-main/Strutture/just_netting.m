% CALCOLO SPESSORE SERBATOIO N2O IN COMPOSITO
% Thrust2 project - Laboratorio di propulsione

%% Init Variables

rho_lim = 3000; 
% Tensione massima ply UD direzione fibre (MPa)
% tens.max di una fibra media scarsa
p = 5.6; % Pressione interna serbatoio prevista (MPa)
f = 3; % Safety Factor
P = p*f; % Pressione di design serbatoio (MPa)
D = 140; % Diametro (interno) serbatoio (mm)

%%  Helical Plies

alphas = 10:10:80; % Angolo delle plies oblique (elicoidali) (°)
tvec = zeros(1,length(alphas));

for i = 1:1:length(alphas)

alpha = alphas(i);

talpha = P*D/4/rho_lim/(cosd(alpha))^2;  
% Spessore plies oblique con angolo alpha (mm)

%%  Hoop Plies

thoop = (P*D*5/12/rho_lim);  
% Spessore plies ortogonali (mm)

%%  Spessore totale

t = thoop + talpha;
% Spessore complessivo plies (mm)

tvec(i) = t;
end

%%  Printing Results

% fprintf("\nSpessore Hoop: %2.4f mm\n",thoop)
% fprintf("Spessore Heli: %2.4f mm\n",talpha)
% fprintf("Spessore  Tot: %2.4f mm\n\n",t)

%% Graphs

figure()
plot(alphas,tvec)
title("Winding Angles vs Spessore")
xlabel("Angolo [°]")
ylabel("Spessore [mm]")
