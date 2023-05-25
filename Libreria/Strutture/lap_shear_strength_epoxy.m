clear all
close all
clc

P = 6e6; 		    % [Pa] pressione in serbatoio
r_int = 0.075; 		    % [m] raggio sferico interno testa 
A_int = pi*r_int^2;     % [m^2] area sezione testa sferica, ovvero l'area 
                        % su cui la pressione agente genera una tensione longitudinale
sigma = 6e6; 		    % [Pa] lap shear strength epoxy resin
F_long = A_int*P; 	    % [N] Forza "stappante"
FoS = 1.25;                % Fattore di sicurezza
A_req = FoS*F_long/sigma; 	% [m^2] Area necessaria per la tenuta dell'incollaggio
r_ext = 0.075; 		    % [m] raggio cilindrico esterno testa 
Circ = 2*pi*r_int; 	    % [m] circonferenza cilindrica esterna testa 
H = A_req/Circ; 	    % [m] altezza minima di incollaggio a partire da saldatura testa

fprintf("Altezza di incollaggio: H = %2.4f mm", H*1000)