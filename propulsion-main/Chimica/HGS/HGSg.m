function [G] = HGSg(S,H,T)
% HGS - Thermochemistry of gas mixture
%
% Property of THRUST, unauthorized distribution is not allowed
% version: HGSg 2.1
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
% [G] = HGSg(S,H,T)
%
%**************************************************************************
%
% HGSg calculates the species free Gibbs energy using enthalpy, enthropy 
% and temperature
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% S --> [kJ/(mol*K)] Entropy
% H --> [kJ/mol] Enthalpy
% T --> [K] Temperature
%
% Outputs:
%--------------------------------------------------------------------------
% G --> [kJ/mol] Free Gibbs energy
%
%**************************************************************************

%{
Changelog:
  > version: 2.1 - 14/12/2022 - Alessandro Rampazzo
     - branched from the original ManelSoria/HGS repository
%}

G = H - T*S; % [kJ/mol]

end