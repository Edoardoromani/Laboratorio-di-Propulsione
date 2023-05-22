function [S] = HGSs(a,T,P,state)
% HGS - Thermochemistry of gas mixture
%
% Property of THRUST, unauthorized distribution is not allowed
% version: HGSs 2.1
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % Branched and modified from the original ManelSoria/HGS repository %
% % Original head:                                                    %
% % *HGS 2.1                                                          % 
% % *By Caleb Fuster, Manel Soria and Arnau Mir�                      %
% % *ESEIAAT UPC                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%**************************************************************************
%
% [S] = HGSs(a,T,P,state)
%
%**************************************************************************
%
% HGSh calculates the enthropy of a species using his Burcat coeficients
% and temperature
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% a -->  Burcat coefficients
% T --> [K] Temperature
% P --> [bar] Pressure
% state -->  State of the species
%
% Outputs:
%--------------------------------------------------------------------------
% S --> [kJ/(mol*K)] Entropy
%
%**************************************************************************

%{
Changelog:
  > version: 2.1 - 14/12/2022 - Alessandro Rampazzo
     - branched from the original ManelSoria/HGS repository
%}

global R; HGSr

Pref = 1; % [bar]

S = R * (a(7) + a(1)*log(T) + sum(a(2:5).*(T.^(1:4))./(1:4))); % [kJ/mol*K]

if strcmp('G',state) &&  P ~= 0
    S = S-R*log(P/Pref); 
end

end