function [H] = HGSh(a,T)
% HGS - Thermochemistry of gas mixture
%
% Property of THRUST, unauthorized distribution is not allowed
% version: HGSh 2.1
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
% [H] = HGSh(a,T)
%
%**************************************************************************
% 
% HGSh calculates the enthalpy of a species using his Burcat coeficients and
% temperature
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% a --> Burcat coefficients
% T --> [K] Temperature
%
% Outputs:
%--------------------------------------------------------------------------
% H --> [kJ/mol] Enthalpy
%
%**************************************************************************

%{
Changelog:
  > version: 2.1 - 14/12/2022 - Alessandro Rampazzo
     - branched from the original ManelSoria/HGS repository
%}  

global R; HGSr

H = R * (a(6) + sum(a(1:5).*(T.^(1:5))./(1:5))); % [kJ/mol]

end