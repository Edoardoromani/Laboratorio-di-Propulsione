clear all
close all
clc
%% VERIFICHE TENSIONI AMMISSIBILI

% VERIFICA TAGLIO VITI
n_fori = 10;                          % numero fori sulla flangia
P = 40*10^5;                          % pressione di progetto
FOS = 1.25;                           % safety factor
rag = 0.075;                          % raggio piastra
F = P*FOS*pi*rag^2/n_fori;            % carico di taglio su ogni vite
ARES = 8.78;                          % area resistente
tau = F/ARES;                         % tensione a taglio su ogni vite
tau_adm = 330;                        % tensione ammisibile di taglio materiale vite
if tau<tau_adm
    disp('verifica a taglio viti OK')
else
    disp('verifica a taglio viti NO BUENO')
end

% VERIFICA A TRAZIONE GUSCIO FLANGIA
r = 75;                  % raggio esterno [mm]
t = 4;                   % spessore guscio 
d = 6;                   % diametro fori
A = pi*r^2 - pi*(r-t)^2 - n_fori*d*t ; % area resistente alla trazione
sigma_fl = F*n_fori/A;           % tensione di trazione/compressione sulla flangia
sigma_all_adm = 290;             % tensione a rottura a trazione N/mm2 della flangia (NON DELLE VITI !!!)
if sigma_fl<sigma_all_adm
    disp('verifica a trazione guscio flangia OK')
else
    disp('verifica a trazione guscio flangia NO BUENO')
end

% VERIFICA RIFOLLAMENTO
t = 8;                       % spessore laterale flangia
a = 8.35;                    % distanza centro foro - lato superiore 
alpha = a/d;      
sigma_rif = F/t/d;             % tensione di trazione sul bordo del foro
if sigma_rif<sigma_all_adm
    disp('verifica rifollameto OK')
else
    disp('verifica rifollamento NO BUENO')
end

% VERIFICA A TAGLIO GUSCIO 
dist = 15;                   % distanza dal bordo superiore e linea di contatto vite
AreaRes = dist*t;            % area resistente
tau_g = F/AreaRes;
tau_admg = sigma_all_adm/sqrt(3);
if tau_g<tau_admg
    disp ('verifica taglio guscio OK')
else 
    disp ('verifica taglio guscio NO BUENO')
end
%% VERIFICHE STATI LIMITE

% VERIFICA VITI A TAGLIO
sigma_rott = 1000;
Fv = 0.5*sigma_rott*ARES/FOS;        % carico a taglio ammisibile vite 0.5 PER CLASSE 6.8 E 10.9, ALTRIMENTI 0.6
if F<Fv
    disp('verifica a taglio viti LIMITE OK')
else
    disp('verifica a taglio viti LIMITE NO BUENO')
end

% VERIFICA A TRAZIONE FLANGIA
Nt = 0.9*sigma_all_adm*A/FOS;       % carico a tr/com ammissibile flangia
if F*8<Nt
    disp('verifica a trazione flangia LIMITE OK')
else
    disp('verifica a trazione flangia LIMITE NO BUENO')
end

% VERIFICA RIFOLLAMENTO
e1 = a;
p2 = 57.3341;   % distanza tra centri dei fori
if d < 15
    d0 = d + 1;
elseif  d > 27
    d0 = d + 3;
    if d < 25 && d > 15
    d0 = d + 2;
    end

end

ALPHA = [e1/3/d0, tau_adm/sigma_all_adm, 1];
K = [1.4*p2/d0 - 1.7, 2.5];
alphaL = min(ALPHA);
k = min(K);

Fbrd = k*alphaL*sigma_all_adm*d*t/1.25;         % carico ammissibile sul bordo fori
if F<Fbrd
    disp('verifica  rifollamento LIMITE OK')
else
    disp('verifica rifollamento LIMITE NO BUENO')
end
%% RISULTATI ADM
TV = tau/tau_adm
TF = sigma_fl/sigma_all_adm
TR = sigma_rif/sigma_all_adm
TVL = F/Fv
TFL = F*8/Nt
TRL = F/Fbrd
PP = tau_g/tau_admg