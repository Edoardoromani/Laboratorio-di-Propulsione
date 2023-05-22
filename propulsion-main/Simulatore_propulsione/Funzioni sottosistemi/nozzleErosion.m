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
  This code simulates the erosion of a graphite nozzle using some empirical
  formulas obtained by online research papers on the matter

Changelog:
  > version: 1.0 - 27/11/2022 - Alessandro Rampazzo
    - added heading, description, and changelog

  > version: 0.1 - 15/03/2022 - Matteo Fiorio
    - created
%}

function rDot_t = nozzleErosion(t,p_cc,A_t,mDot_tot,O_F,rocket,env)

% stoichiometric equivalence ratio
O_F_st = 7.5;

% equivalence ratio
phi = O_F / O_F_st;

if t >= rocket.nozzle.onsetTime && rocket.nozzle.erosionModel >= 0
    if rocket.nozzle.erosionModel == 0
        % best interpolation of data from "Erosion Rate Investigation of
        % Various Nozzle Materials in Hybrid Rocket Motors"
        m = 0.6316;
        a = 0.0684e-3;
        n = 0;
        rDot_t = a * phi.^m * (mDot_tot./A_t).^n;

    elseif rocket.nozzle.erosionModel == 1
        % max value from interpolation of data from "Erosion Rate 
        % Investigation of Various Nozzle Materials in Hybrid Rocket Motors"
        m = 0.6316;
        a = 0.12e-3;
        n = 0;
        rDot_t = a * phi.^m * (mDot_tot./A_t).^n;
    elseif rocket.nozzle.erosionModel == 2
        % from "Pyrolytic Graphite and Boron Nitride as Low-Erosion Nozzle
        % Materials for Long-Duration Hybrid Rocket Testing"
        rDot_t = 0.05e-3;
    elseif rocket.nozzle.erosionModel == 3
        % correlation for solid rockets
        a = 5.234e-6;
        n = 0.8;
        rDot_t = a* O_F * (p_cc/1e+5)^n;
    else
        error("erosion model not implemented")
    end
else
    rDot_t = 0;
end
end

