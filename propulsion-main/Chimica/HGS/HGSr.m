function HGSr
% HGS - Thermochemistry of gas mixture
%
% Property of THRUST, unauthorized distribution is not allowed
% version: HGSr 2.1
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
% HGSr
%
%**************************************************************************
% 
% HGSr loads as global variable Universal R
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
%
% Outputs:
%--------------------------------------------------------------------------
% R is load in the global workspace
%
%**************************************************************************

%{
Changelog:
  > version: 2.1 - 14/12/2022 - Alessandro Rampazzo
     - branched from the original ManelSoria/HGS repository
%}   

global R;
if isempty(R) 
    R = 8.3144621/1000;  % [kJ/(mol*K)] 
end

end