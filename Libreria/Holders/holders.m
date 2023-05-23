%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      %
%              THRUST                  %
%    Holder superiore ed inferiore     %
%                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc 
clear
close all

% book: Pressure Vessel Design Manual 3rd Ed.
% pdf: cap16

% Data
p_tank = 60e5; % [bar] -> [Pa]
fs = 2; % safety factor
p1 = (fs * p_tank) * 1e-6; % [Pa] -> [N/mm^2], pressione operativa nel caso peggiore (all'accensione)

% Avendo un'area A1 dove agisce la pressione p1 è possibile trovare la forza F1 agente su di essa
% Conoscendo la forza si trova la pressione effettiva in cc dato il rapporto tra le aree

d1 = 20; % [mm], diametro
r1 = d1 / 2;
A1 = r1^2 * pi; % [mm^2], area raccordo valvola
F1 = p1 * A1; % [N], (= Ht nel book)

d2 = 25; % [mm], raggio interno dell'area effettiva dove agisce la pressione (= Bn nel book)
r2 = d2 / 2;
d3 = 50; % [mm], raggio esterno dell'area effettiva dove agisce la pressione 
r3 = d3 / 2;
A2 = pi*(r3^2 - r2^2); % [mm^2], area effettiva dove agisce la pressione

p2 = (F1 / A2) ; % [N/mm^2]

%% HOLDER SUPERIORE in alluminio AL6082 - apertura senza nozzle (= loose), (VEDI BOOK P.78)
% dimensioni -> vedi immagini pag.78 (book)
t = 9; % [mm], spessore verticale holder
g0 = 4; % [mm], spessore dell'holder
g1 = 4; % [mm], parte orizzontale raccordino
h = 3; % [mm], parte verticale raccordino
d_ext = 112; % [mm], diametro esterno holder (= A nel book)
d_int = d_ext - 2*g0; % [mm], diametro interno holder (= Bs nel book)

K = d_ext / d2; % rapporto diametri
Z = (K^2 + 1) / (K^2 - 1); % fattore K
hd = (d_ext - d2) / 2;
ht = (d_int - d2) / 4 * g0/2;

Hd = (pi * d2^2 * p2) / 4; % [N], forza idrostatica sulla parte centrale
H = (pi * d3^2 * p2) / 4; % [N], forza idrostatica
Ht = H - Hd; % [N], uguale a F1 -> ok

% Momenti
Md = hd * Hd; % [Nmm]
Mt = ht * Ht; % [Nmm]
M0 = Mt + Md; % [Nmm]

% scelta materiale -> alluminio o acciaio
ni = 0.33; % fattore di Poisson alluminio 

% alluminio Al6082
sigma_adm = 250; % [Mpa]
tau_adm = 210; % [MPa]

% step 1 p.78 -> calcolo fattori adimensionali
U = (K^2 * (1 + 4.6052 * (1 + ni) / (1 - ni) * log10(K))-1) / (1.0472 * (K^2 - 1) * (K - 1) * (1 + ni)); 
Y = (1 - ni^2) * U;
h0 = sqrt(d2 * g0);
rapp1 = g1 / g0; % s2/s1 = 1
rapp2 = h / h0; % 0.3

% step 2 p.78 
V = 0.56; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
F = 0.91; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
f = 0.5; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
e =  F / h0; % fattore pag 38 (book)
d = (U * h0 * g0^2) / V;
T = ((1 - ni^2) * (K^2 - 1 )) * U / ((1 - ni) + (1 + ni) * K^2); % fattore pag 38 (book)
L = (t * e + 1) / T + t^3 / d;

St = Y * M0 / (t^2 * d2); % [N/mm^2]
Sr = ((1.33 * e * t + 1) * M0) / (L * t^2 * d2);
Sh = f * M0 / (L * g1^2 * d2);
theta = d2 * St / t; % fattore
B1 = d2 + g0; % [mm]
Mh = theta / (1.74 * h0 * V / (g0^3 * B1) + theta / M0 * (1 + F * t / h0)); % [Nmm], momento nella giunzione tra holder e cc
X1 = (M0 - Mh * (1 + F * t / h0)) / M0; % fattore

% calcolo stress finali
S_hs = (1.1 * X1 * theta * f * h0) / ((g1/g0)^2 * d_int * V); % [N/mm^2], stress longitudinale
S_rs = (1.91 * Mh * (1 + (F * t) / h0)) / (d_int * t^2) + (0.64 * F * Mh) / (d_int * h0 * t); % [N/mm^2], stress radiale
S_ts = (X1 * theta * t) / d_int - (0.57 * Mh * (1 + (F * t) / h0)) / (d_int * t^2) + (0.64 * F * Z * Mh) / (d_int * h0 * t); % [N/mm^2], stress tangenziale

% MoS
FoS = 1.2;
MoS_hs_sigma = sigma_adm / (S_hs * FoS) - 1; % margine di sicurezza head-shell juncture
MoS_s_sigma = sigma_adm / (Sh * FoS) - 1; % margine di sicurezza shell
MoS_h_tau = tau_adm / ((St+Sr) * FoS) - 1; % margine di sicurezza
MoS_hs_tau = tau_adm / ((S_ts+S_rs) * FoS) - 1; % margine di sicurezza

% output
disp('<< Holder Superiore - Alluminio Al6082 >>')
fprintf("t: %1.2f [mm]\n", t)
fprintf("g0: %1.2f [mm]\n", g0)
fprintf("g1: %1.2f [mm]\n\n", g1)
fprintf("Sigma adm: %1.2f [MPa]\n", sigma_adm) % per verifica a trazione/compressione
fprintf("Tau adm: %1.2f [MPa]\n\n", tau_adm) % per verifica a taglio

disp('Head')
fprintf("Tau head: %1.2f [MPa]\n\n", St+Sr)

disp('Head-Shell juncture')
fprintf("Sigma h-s juncture: %1.2f [MPa]\n", S_hs)
fprintf("Tau h-s juncture: %1.2f [MPa]\n\n", S_ts+S_rs)

disp('Shell')
fprintf("Sigma shell: %1.2f [MPa]\n\n\n", Sh)


%% HOLDER SUPERIORE in acciaio 304L - apertura senza nozzle (= loose), (VEDI BOOK P.78)
% dimensioni -> vedi immagini pag.78 (book)
t = 13; % [mm], spessore verticale holder
g0 = 8; % [mm], spessore dell'holder
g1 = 10; % [mm], parte orizzontale raccordino
h = 3; % [mm], parte verticale raccordino
d_ext = 112; % [mm], diametro esterno holder (= A nel book)
d_int = d_ext - 2*g0; % [mm], diametro interno holder (= Bs nel book)

K = d_ext / d2; % rapporto diametri
Z = (K^2 + 1) / (K^2 - 1); % fattore K
hd = (d_ext - d2) / 2;
ht = (d_int - d2) / 4 * g0/2;

Hd = (pi * d2^2 * p2) / 4; % [N], forza idrostatica sulla parte centrale
H = (pi * d3^2 * p2) / 4; % [N], forza idrostatica
Ht = H - Hd; % [N], uguale a F1 -> ok

% Momenti
Md = hd * Hd; % [Nmm]
Mt = ht * Ht; % [Nmm]
M0 = Mt + Md; % [Nmm]

% scelta materiale -> alluminio o acciaio
ni = 0.3; % fattore di Poisson acciaio

% acciaio 304L
fs_sigma = 1.5; % fattore di sicurezza sulla sigma
sigma_adm = 882; % [MPa]
tau_adm = sigma_adm * 0.577 * fs_sigma; % [MPa]

% step 1 p.78 -> calcolo fattori adimensionali
U = (K^2 * (1 + 4.6052 * (1 + ni) / (1 - ni) * log10(K))-1) / (1.0472 * (K^2 - 1) * (K - 1) * (1 + ni)); 
Y = (1 - ni^2) * U;
h0 = sqrt(d2 * g0);
rapp1 = g1 / g0;
rapp2 = h / h0;

% step 2 p.78 
V = 0.48; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
F = 0.9; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
f = 1.2; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
e =  F / h0; % fattore pag 38 (book)
d = (U * h0 * g0^2) / V;
T = ((1 - ni^2) * (K^2 - 1 )) * U / ((1 - ni) + (1 + ni) * K^2); % fattore pag 38 (book)
L = (t * e + 1) / T + t^3 / d;

St = Y * M0 / (t^2 * d2); % [N/mm^2]
Sr = ((1.33 * e * t + 1) * M0) / (L * t^2 * d2);
Sh = f * M0 / (L * g1^2 * d2);
theta = d2 * St / t; % fattore
B1 = d2 + g0; % [mm]
Mh = theta / (1.74 * h0 * V / (g0^3 * B1) + theta / M0 * (1 + F * t / h0)); % [Nmm], momento nella giunzione tra holder e cc
X1 = (M0 - Mh * (1 + F * t / h0)) / M0; % fattore

% calcolo stress finali
S_hs = (1.1 * X1 * theta * f * h0) / ((g1/g0)^2 * d_int * V); % [N/mm^2], stress longitudinale
S_rs = (1.91 * Mh * (1 + (F * t) / h0)) / (d_int * t^2) + (0.64 * F * Mh) / (d_int * h0 * t); % [N/mm^2], stress radiale
S_ts = (X1 * theta * t) / d_int - (0.57 * Mh * (1 + (F * t) / h0)) / (d_int * t^2) + (0.64 * F * Z * Mh) / (d_int * h0 * t); % [N/mm^2], stress tangenziale

% MoS
FoS = 1.2;
MoS_hs_sigma = sigma_adm / (S_hs * FoS) - 1; % margine di sicurezza head-shell juncture
MoS_s_sigma = sigma_adm / (Sh * FoS) - 1; % margine di sicurezza shell
MoS_h_tau = tau_adm / ((St+Sr) * FoS) - 1; % margine di sicurezza
MoS_hs_tau = tau_adm / ((S_ts+S_rs) * FoS) - 1; % margine di sicurezza

% output
disp('<< Holder Superiore - Acciaio >>')
fprintf("t: %1.2f [mm]\n", t)
fprintf("g0: %1.2f [mm]\n", g0)
fprintf("g1: %1.2f [mm]\n\n", g1)
fprintf("Sigma adm: %1.2f [MPa]\n", sigma_adm) % per verifica a trazione/compressione
fprintf("Tau adm: %1.2f [MPa]\n\n", tau_adm) % per verifica a taglio

disp('Head')
fprintf("Tau head: %1.2f [MPa]\n\n", St+Sr)

disp('Head-Shell juncture')
fprintf("Sigma h-s juncture: %1.2f [MPa]\n", S_hs)
fprintf("Tau h-s juncture: %1.2f [MPa]\n\n", S_ts+S_rs)

disp('Shell')
fprintf("Sigma shell: %1.2f [MPa]\n\n\n", Sh)


%% HOLDER INFERIORE in alluminio AL6082 - apertura con nozzle (= si suppone loose pure qui), (VEDI BOOK P.78)
% dimensioni -> vedi immagini pag.78 (book)
t = 14; % [mm], spessore verticale holder
g0 = 14; % [mm], spessore dell'holder
g1 = 23; % [mm], parte orizzontale raccordino
h = 11; % [mm], parte verticale raccordino
d_ext = 112; % [mm], diametro esterno holder (= A nel book)
d_int = d_ext - 2*g0; % [mm], diametro interno holder (= Bs nel book)

K = d_ext / d2; % rapporto diametri
Z = (K^2 + 1) / (K^2 - 1); % fattore K
hd = (d_ext - d2) / 2;
ht = (d_int - d2) / 4 * g0/2;

Hd = (pi * d2^2 * p2) / 4; % [N], forza idrostatica sulla parte centrale
H = (pi * d3^2 * p2) / 4; % [N], forza idrostatica
Ht = H - Hd; % [N], uguale a F1 -> ok

% Momenti
Md = hd * Hd; % [Nmm]
Mt = ht * Ht; % [Nmm]
M0 = Mt + Md; % [Nmm]

% scelta materiale -> alluminio o acciaio
ni = 0.33; % fattore di Poisson alluminio 

% alluminio Al6082
sigma_adm = 250; % [Mpa]
tau_adm = 210; % [MPa]

% step 1 p.78 -> calcolo fattori adimensionali
U = (K^2 * (1 + 4.6052 * (1 + ni) / (1 - ni) * log10(K))-1) / (1.0472 * (K^2 - 1) * (K - 1) * (1 + ni)); 
Y = (1 - ni^2) * U;
h0 = sqrt(d2 * g0);
rapp1 = g1 / g0;
rapp2 = h / h0;

% step 2 p.78 
V = 0.34; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
F = 0.86; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
f = 1; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
e =  F / h0; % fattore pag 38 (book)
d = (U * h0 * g0^2) / V;
T = ((1 - ni^2) * (K^2 - 1 )) * U / ((1 - ni) + (1 + ni) * K^2); % fattore pag 38 (book)
L = (t * e + 1) / T + t^3 / d;

St = Y * M0 / (t^2 * d2); % [N/mm^2]
Sr = ((1.33 * e * t + 1) * M0) / (L * t^2 * d2);
Sh = f * M0 / (L * g1^2 * d2);
theta = d2 * St / t; % fattore
B1 = d2 + g0; % [mm]
Mh = theta / (1.74 * h0 * V / (g0^3 * B1) + theta / M0 * (1 + F * t / h0)); % [Nmm], momento nella giunzione tra holder e cc
X1 = (M0 - Mh * (1 + F * t / h0)) / M0; % fattore

% calcolo stress finali
S_hs = (1.1 * X1 * theta * f * h0) / ((g1/g0)^2 * d_int * V); % [N/mm^2], stress longitudinale
S_rs = (1.91 * Mh * (1 + (F * t) / h0)) / (d_int * t^2) + (0.64 * F * Mh) / (d_int * h0 * t); % [N/mm^2], stress radiale
S_ts = (X1 * theta * t) / d_int - (0.57 * Mh * (1 + (F * t) / h0)) / (d_int * t^2) + (0.64 * F * Z * Mh) / (d_int * h0 * t); % [N/mm^2], stress tangenziale

% MoS
FoS = 1.2;
MoS_hs_sigma = sigma_adm / (S_hs * FoS) - 1; % margine di sicurezza head-shell juncture
MoS_s_sigma = sigma_adm / (Sh * FoS) - 1; % margine di sicurezza shell
MoS_h_tau = tau_adm / ((St+Sr) * FoS) - 1; % margine di sicurezza
MoS_hs_tau = tau_adm / ((S_ts+S_rs) * FoS) - 1; % margine di sicurezza

% output
disp('<< Holder Inferiore - Alluminio Al6082 >>')
fprintf("t: %1.2f [mm]\n", t)
fprintf("g0: %1.2f [mm]\n", g0)
fprintf("g1: %1.2f [mm]\n\n", g1)
fprintf("Sigma adm: %1.2f [MPa]\n", sigma_adm) % per verifica a trazione/compressione
fprintf("Tau adm: %1.2f [MPa]\n\n", tau_adm) % per verifica a taglio

disp('Head')
fprintf("Tau head: %1.2f [MPa]\n\n", St+Sr)

disp('Head-Shell juncture')
fprintf("Sigma h-s juncture: %1.2f [MPa]\n", S_hs)
fprintf("Tau h-s juncture: %1.2f [MPa]\n\n", S_ts+S_rs)

disp('Shell')
fprintf("Sigma shell: %1.2f [MPa]\n\n\n", Sh)


%% HOLDER INFERIORE in acciaio 304L - apertura con nozzle (= si suppone loose pure qui), (VEDI BOOK P.78)
% dimensioni -> vedi immagini pag.78 (book)
t = 28; % [mm], spessore verticale holder
g0 = 7; % [mm], spessore dell'holder
g1 = 23; % [mm], parte orizzontale raccordino
h = 11; % [mm], parte verticale raccordino
d_ext = 112; % [mm], diametro esterno holder (= A nel book)
d_int = d_ext - 2*g0; % [mm], diametro interno holder (= Bs nel book)

K = d_ext / d2; % rapporto diametri
Z = (K^2 + 1) / (K^2 - 1); % fattore K
hd = (d_ext - d2) / 2;
ht = (d_int - d2) / 4 * g0/2;

Hd = (pi * d2^2 * p2) / 4; % [N], forza idrostatica sulla parte centrale
H = (pi * d3^2 * p2) / 4; % [N], forza idrostatica
Ht = H - Hd; % [N], uguale a F1 -> ok

% Momenti
Md = hd * Hd; % [Nmm]
Mt = ht * Ht; % [Nmm]
M0 = Mt + Md; % [Nmm]

% scelta materiale -> alluminio o acciaio
ni = 0.3; % fattore di Poisson acciaio 

% acciaio 304L
fs_sigma = 1.5; % fattore di sicurezza sulla sigma
sigma_adm = 176; % [Mpa]
tau_adm = sigma_adm * 0.577 * fs_sigma; % [MPa]

% step 1 p.78 -> calcolo fattori adimensionali
U = (K^2 * (1 + 4.6052 * (1 + ni) / (1 - ni) * log10(K))-1) / (1.0472 * (K^2 - 1) * (K - 1) * (1 + ni)); 
Y = (1 - ni^2) * U;
h0 = sqrt(d2 * g0);
rapp1 = g1 / g0;
rapp2 = h / h0;

% step 2 p.78 
V = 0.11; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
F = 0.76; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
f = 2.4; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
e =  F / h0; % fattore pag 38 (book)
d = (U * h0 * g0^2) / V;
T = ((1 - ni^2) * (K^2 - 1 )) * U / ((1 - ni) + (1 + ni) * K^2); % fattore pag 38 (book)
L = (t * e + 1) / T + t^3 / d;

St = Y * M0 / (t^2 * d2); % [N/mm^2]
Sr = ((1.33 * e * t + 1) * M0) / (L * t^2 * d2);
Sh = f * M0 / (L * g1^2 * d2);
theta = d2 * St / t; % fattore
B1 = d2 + g0; % [mm]
Mh = theta / (1.74 * h0 * V / (g0^3 * B1) + theta / M0 * (1 + F * t / h0)); % [Nmm], momento nella giunzione tra holder e cc
X1 = (M0 - Mh * (1 + F * t / h0)) / M0; % fattore

% calcolo stress finali
S_hs = (1.1 * X1 * theta * f * h0) / ((g1/g0)^2 * d_int * V); % [N/mm^2], stress longitudinale
S_rs = (1.91 * Mh * (1 + (F * t) / h0)) / (d_int * t^2) + (0.64 * F * Mh) / (d_int * h0 * t); % [N/mm^2], stress radiale
S_ts = (X1 * theta * t) / d_int - (0.57 * Mh * (1 + (F * t) / h0)) / (d_int * t^2) + (0.64 * F * Z * Mh) / (d_int * h0 * t); % [N/mm^2], stress tangenziale

% MoS
FoS = 1.2;
MoS_hs_sigma = sigma_adm / (S_hs * FoS) - 1; % margine di sicurezza head-shell juncture
MoS_s_sigma = sigma_adm / (Sh * FoS) - 1; % margine di sicurezza shell
MoS_h_tau = tau_adm / ((St+Sr) * FoS) - 1; % margine di sicurezza
MoS_hs_tau = tau_adm / ((S_ts+S_rs) * FoS) - 1; % margine di sicurezza

% output
disp('<< Holder Inferiore - Acciaio 304L >>')
fprintf("t: %1.2f [mm]\n", t)
fprintf("g0: %1.2f [mm]\n", g0)
fprintf("g1: %1.2f [mm]\n\n", g1)
fprintf("Sigma adm: %1.2f [MPa]\n", sigma_adm) % per verifica a trazione/compressione
fprintf("Tau adm: %1.2f [MPa]\n\n", tau_adm) % per verifica a taglio

disp('Head')
fprintf("Tau head: %1.2f [MPa]\n\n", St+Sr)

disp('Head-Shell juncture')
fprintf("Sigma h-s juncture: %1.2f [MPa]\n", S_hs)
fprintf("Tau h-s juncture: %1.2f [MPa]\n\n", S_ts+S_rs)

disp('Shell')
fprintf("Sigma shell: %1.2f [MPa]\n\n\n", Sh)

%% HOLDER INFERIORE in acciaio - apertura con nozzle (= si suppone loose pure qui), (VEDI BOOK P.78)
% dimensioni -> vedi immagini pag.78 (book)
t = 20; % [mm], spessore verticale holder
g0 = 7; % [mm], spessore dell'holder
g1 = 23; % [mm], parte orizzontale raccordino
h = 11; % [mm], parte verticale raccordino
d_ext = 112; % [mm], diametro esterno holder (= A nel book)
d_int = d_ext - 2*g0; % [mm], diametro interno holder (= Bs nel book)

K = d_ext / d2; % rapporto diametri
Z = (K^2 + 1) / (K^2 - 1); % fattore K
hd = (d_ext - d2) / 2;
ht = (d_int - d2) / 4 * g0/2;

Hd = (pi * d2^2 * p2) / 4; % [N], forza idrostatica sulla parte centrale
H = (pi * d3^2 * p2) / 4; % [N], forza idrostatica
Ht = H - Hd; % [N], uguale a F1 -> ok

% Momenti
Md = hd * Hd; % [Nmm]
Mt = ht * Ht; % [Nmm]
M0 = Mt + Md; % [Nmm]

% scelta materiale -> alluminio o acciaio
ni = 0.3; % fattore di Poisson acciaio 

% acciaio
fs_sigma = 1.5; % fattore di sicurezza sulla sigma
sigma_adm = 300; % [Mpa]
tau_adm = sigma_adm * 0.577 * fs_sigma; % [MPa]

% step 1 p.78 -> calcolo fattori adimensionali
U = (K^2 * (1 + 4.6052 * (1 + ni) / (1 - ni) * log10(K))-1) / (1.0472 * (K^2 - 1) * (K - 1) * (1 + ni)); 
Y = (1 - ni^2) * U;
h0 = sqrt(d2 * g0);
rapp1 = g1 / g0;
rapp2 = h / h0;

% step 2 p.78 
V = 0.11; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
F = 0.76; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
f = 2.4; % fattore ottenuto da grafico in pdf inserendo rapp1 e rapp2
e =  F / h0; % fattore pag 38 (book)
d = (U * h0 * g0^2) / V;
T = ((1 - ni^2) * (K^2 - 1 )) * U / ((1 - ni) + (1 + ni) * K^2); % fattore pag 38 (book)
L = (t * e + 1) / T + t^3 / d;

St = Y * M0 / (t^2 * d2); % [N/mm^2]
Sr = ((1.33 * e * t + 1) * M0) / (L * t^2 * d2);
Sh = f * M0 / (L * g1^2 * d2);
theta = d2 * St / t; % fattore
B1 = d2 + g0; % [mm]
Mh = theta / (1.74 * h0 * V / (g0^3 * B1) + theta / M0 * (1 + F * t / h0)); % [Nmm], momento nella giunzione tra holder e cc
X1 = (M0 - Mh * (1 + F * t / h0)) / M0; % fattore

% calcolo stress finali
S_hs = (1.1 * X1 * theta * f * h0) / ((g1/g0)^2 * d_int * V); % [N/mm^2], stress longitudinale
S_rs = (1.91 * Mh * (1 + (F * t) / h0)) / (d_int * t^2) + (0.64 * F * Mh) / (d_int * h0 * t); % [N/mm^2], stress radiale
S_ts = (X1 * theta * t) / d_int - (0.57 * Mh * (1 + (F * t) / h0)) / (d_int * t^2) + (0.64 * F * Z * Mh) / (d_int * h0 * t); % [N/mm^2], stress tangenziale

% MoS
FoS = 1.2;
MoS_hs_sigma = sigma_adm / (S_hs * FoS) - 1; % margine di sicurezza head-shell juncture
MoS_s_sigma = sigma_adm / (Sh * FoS) - 1; % margine di sicurezza shell
MoS_h_tau = tau_adm / ((St+Sr) * FoS) - 1; % margine di sicurezza
MoS_hs_tau = tau_adm / ((S_ts+S_rs) * FoS) - 1; % margine di sicurezza

% output
disp('<< Holder Inferiore - Acciaio >>')
fprintf("t: %1.2f [mm]\n", t)
fprintf("g0: %1.2f [mm]\n", g0)
fprintf("g1: %1.2f [mm]\n\n", g1)
fprintf("Sigma adm: %1.2f [MPa]\n", sigma_adm) % per verifica a trazione/compressione
fprintf("Tau adm: %1.2f [MPa]\n\n", tau_adm) % per verifica a taglio

disp('Head')
fprintf("Tau head: %1.2f [MPa]\n\n", St+Sr)

disp('Head-Shell juncture')
fprintf("Sigma h-s juncture: %1.2f [MPa]\n", S_hs)
fprintf("Tau h-s juncture: %1.2f [MPa]\n\n", S_ts+S_rs)

disp('Shell')
fprintf("Sigma shell: %1.2f [MPa]\n\n\n", Sh)




