%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  _____  __  __  ____  _   _  ____  _____                %
%                 |_   _||  ||  || __ \| | | |/ ___||_   _|               %
%                   | |  |  __  ||    /| |_| |\___ \  | |                 %
%                   |_|  |__||__||_|\_\ \___/ |____/  |_|                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Project: 1D Thermal Transient Simulator

Property of THRUST, unauthorized distribution is not allowed

Description:
  This function returns the regression rate in [m/s] given the material ID 
  and the temperature

Changelog:
  > version: 1.0 - 13/11/2022 - Alessandro Rampazzo
    - created
%}

function rDot = getRegressionRate(materialID,t,T)

switch materialID
    case 0
        %% graphite

        rDot = 0;

    case 1
        %% steel

        rDot = 0;

    case 2
        %% phenolic resin Si based

        rDot = 0.3e-3;

    case 3
        %% aluminum alloy AA 7075 - T6 (Ergal)

        rDot = 0;

    case 4
        %% paraffin sasolwax 0907

        rDot = 1.34e-3;

    case 5
        %% high density polyethylene (HDPE)

        rDot = 0.5e-3;

    case 6
        %% aluminum alloy AA 6082 - T6

        rDot = 0;

    case 7
        %% aluminum alloy AA 6061 - T6

        rDot = 0;
        
end
end