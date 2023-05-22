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
  This code implements the function to simulate the grain regression rate

Changelog:
  > version: 1.0 - 27/11/2022 - Alessandro Rampazzo
    - added heading, description, and changelog

  > version: 0.1 - 15/03/2022 - Alessandro Rampazzo
    - created
%}

function out = grainFunction(mDot_ox,r_port,rocket)
    % port area
    A_port = pi*r_port^2;
    % oxidizer mass flux
    G_ox = mDot_ox/A_port;
    % average fuel regression rate
    rDot_port = rocket.cc.a * G_ox.^rocket.cc.n;
    % burning area
    Ab = rocket.cc.L_grain * 2 * pi * r_port;
    % fuel mass flow rate
    mDot_fuel = rDot_port * Ab * rocket.cc.rho_fuel;
    
    % output variables
    out.mDot_fuel = mDot_fuel;
    out.rDot_port = rDot_port;
end

