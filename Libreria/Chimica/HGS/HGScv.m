function [Cv] = HGScv(Cp)
% HGS - Thermochemistry of gas mixture
%
% Property of THRUST, unauthorized distribution is not allowed
% version: HGScv 2.1
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
% [Cv] = HGScv(Cp)
%
%**************************************************************************
%
% HGScv calculates the species Constant volume coeficient using constant
% pressure coeficient
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% Cp --> [kJ/(mol*K)] Constant pressure coefficient
%
% Outputs:
%--------------------------------------------------------------------------
% Cv --> [kJ/(mol*K)] Constant volume coefficient
%
%**************************************************************************

%{
Changelog:
  > version: 2.1 - 14/12/2022 - Alessandro Rampazzo
     - branched from the original ManelSoria/HGS repository
%}

global R

Cv = Cp - R;  

end