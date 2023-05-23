function  HGSsubtMixture(name)
% HGS - Thermochemistry of gas mixture
%
% Property of THRUST, unauthorized distribution is not allowed
% version: HGSsubtMixture 2.1
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
% HGSsubtMixture(name)
%
%**************************************************************************
%
% HGSsubtMixture eliminates a mixture from the HGSdata database
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% name --> Name of the mixture f.e. 'Air'
% 
% Outputs:
%--------------------------------------------------------------------------
% HGSdata update
%
%**************************************************************************

%{
Changelog:
  > version: 2.1 - 14/12/2022 - Alessandro Rampazzo
     - branched from the original ManelSoria/HGS repository
%}   

global HGSdata; HGSload


% Check if no species or mixture use the name
try
    [id] = HGSid(name);
catch
    error('Ups,.. this name is not used in HGSdata')
end

del = id - length(HGSdata.Mm);

HGSdata.comb(del) = [];
HGSdata.cspec(del) = [];
HGSdata.cper(del) = [];

end