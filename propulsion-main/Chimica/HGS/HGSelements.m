function [eln,elq] = HGSelements(name)
% HGS - Thermochemistry of gas mixture
%
% Property of THRUST, unauthorized distribution is not allowed
% version: HGSelements 2.1
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
% HGSelements(name)
%
%**************************************************************************
%
% HGSelements returns the elements present in the species name and their
% number
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% name --> Name or code of the species
%
% Outputs:
% eln --> List with the name of the elements present in the species
% elq --> Vector with their quantities
%--------------------------------------------------------------------------
% Examples:
% [a,b]=HGSelements('H2O');
% [a,b]=HGSelements(2247);
%**************************************************************************

%{
Changelog:
  > version: 2.1 - 14/12/2022 - Alessandro Rampazzo
     - branched from the original ManelSoria/HGS repository
%}

if isnumeric(name)
    id=name;
else
    id=HGSid(name);
end

global HGSdata;HGSload;

eln=HGSdata.ena(id);
eln=eln{1};
elq=HGSdata.nat(id);
elq=elq{1};
end