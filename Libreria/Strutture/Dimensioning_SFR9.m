%% Codice per ottenere le dimensioni del razzo da lunghezza e
% e diametro. preso come riferimento Paper:Design Methodology and Performance Evaluation of New
%Generation Sounding Rockets

clc
clear all

l_SFR =4.5;     % [m] lunghezza SFR9
D_SFR = 150/1e3;        % [mm] diametro razz o SFR9
Cogive = 0.17;      %Coefficente da paper

l_nosec = Cogive*l_SFR;
l_body = l_SFR-l_nosec;

% Determina superfici aereodinamiche 
Ap = (l_SFR-l_nosec)*D_SFR;
Ac = Ap/45;
At = Ap/12;
h_rott = sqrt((4/3)*At);      % [m] Corda delle alette alla radice
h_tip = h_rott/2;    % [m] Corda alla tip

dimension_SFR = {l_SFR;D_SFR;l_nosec;l_body;h_rott;h_tip};
%Conversione in pollici per rasaero
dimension_SFR_in = [l_SFR;D_SFR;l_nosec;l_body;h_rott;h_tip]*39.37;

% l_SFR_in = l_SFR*39.37;
% l_body = l_body*39.37;
% l_nosec_in = l_nosec*39.37;






