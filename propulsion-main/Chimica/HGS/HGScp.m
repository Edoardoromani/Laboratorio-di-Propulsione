function [Cp] = HGScp(a,T)
% HGS - Thermochemistry of gas mixture
%
% Property of THRUST, unauthorized distribution is not allowed
% version: HGScp 2.1
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % Branched and modified from the original ManelSoria/HGS repository %
% % Original head:                                                    %
% % *HGS 2.1                                                          % 
% % *By Caleb Fuster, Manel Soria and Arnau Miró                      %
% % *ESEIAAT UPC                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%**************************************************************************
%
% [Cp] = HGScp(a,T)
%
%**************************************************************************
%
% HGScp calculates the species Constant pressure coeficient using Burcat
% coeficients  and temperature
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% a --> Burcat coefficients
% T --> [K] Temperature
%
% Outputs:
%--------------------------------------------------------------------------
% Cp --> [kJ/(mol*K)] Constant pressure coeficient
%
%**************************************************************************

%{
Changelog:
  > version: 2.1 - 14/12/2022 - Alessandro Rampazzo
     - branched from the original ManelSoria/HGS repository
%}

global R; HGSr;

Cp = R * (a(1) + sum(a(2:5).*(T.^(1:4)))); % [kJ/mol*K]


end