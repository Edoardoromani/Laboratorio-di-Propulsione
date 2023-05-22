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
  This code returns the name of the material given the ID

Changelog:
  > version: 1.0 - 13/11/2022 - Alessandro Rampazzo
    - created
%}

function name = getName(materialID)

switch materialID
    case 0
        name = "graphite";
    case 1
        name = "steel";
    case 2
        name = "phenolic resin";
    case 3
        name = "AA 7075 - T6 ";
    case 4
        name = "sasolwax 0907";
    case 5
        name = "HDPE";
    case 6
        name = "AA 6082 - T6";
    case 7
        name = "AA 6061 - T6";
    otherwise
        name = "unknown";
end
end

