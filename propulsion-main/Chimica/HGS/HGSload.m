function HGSload
% HGS - Thermochemistry of gas mixture
%
% Property of THRUST, unauthorized distribution is not allowed
% version: HGSload 2.1
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
% HGSload 
%
%**************************************************************************
% 
% HGSload loads as global variable HGSdata
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
%
% Outputs:
%--------------------------------------------------------------------------
% HGSdata is load in the global workspace
%
%**************************************************************************

%{
Changelog:
  > version: 2.1 - 14/12/2022 - Alessandro Rampazzo
     - branched from the original ManelSoria/HGS repository
%}  

global HGSdata;
if ~isstruct(HGSdata) 
    try
        load('HGSdata')
    catch
       error('HGSdata.mat can not be find in the current paths.\n Use HGSdataDownload to create the data') 
    end
end

end

