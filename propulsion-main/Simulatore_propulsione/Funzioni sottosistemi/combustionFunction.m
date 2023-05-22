%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  _____  __  __  ____  _   _  ____  _____                %
%                 |_   _||  ||  || __ \| | | |/ ___||_   _|               %
%                   | |  |  __  ||    /| |_| |\___ \  | |                 %
%                   |_|  |__||__||_|\_\ \___/ |____/  |_|                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Project: Propulsion Simulator

Property of THRUST, unauthorized distribution is not allowed

Description:
  This code implements the model for combustion in a combustion chamber

Changelog:
  > version: 1.1 - 05/12/2022 - Alessandro Rampazzo
    - added T_ox and T_f as imputs

  > version: 1.0 - 04/12/2022 - Alessandro Rampazzo
    - created
%}

function gasProp = combustionFunction(O_F, T_ox, T_f, rocket)
    gasProp.T = rocket.cc.T;
    gasProp.gamma = rocket.cc.gamma;
    gasProp.R = rocket.cc.R;
end