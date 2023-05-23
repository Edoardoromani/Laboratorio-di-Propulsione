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
  This function returns the density in [kg/m^3] given the material ID and
  the temperature

Changelog:
  > version: 1.0 - 13/11/2022 - Alessandro Rampazzo
    - created
%}

function rho = getRho(materialID,T)

switch materialID
    case 0
        %% graphite
        % valid between 200 K ÷ 3000 K

        % from "Thermomechanical Response of Advanced Materials under
        % Quasi-Instantaneous Heating", Least square interpolation
        rho = 1.7971e+03 - 0.0296*T;
    case 1
        %% steel
        % valid between 200 K ÷ 3000 K

        % from https://it.wikipedia.org/wiki/Acciaio_inossidabile
        rho = 7750;
    case 2
        %% phenolic resin Si based
        % valid between ? K ÷ ? K

        % from https://www.researchgate.net/figure/Properties-of-Si-phenolic-resin_tbl1_276225652
        rho = 1730;
    case 3
        %% aluminum alloy AA 7075 - T6 (Ergal)
        % valid between ? K ÷ ? K

        % from https://www.comefimetalli.it/lega7075.asp
        rho = 2810;
    case 4
        %% paraffin sasolwax 0907
        % valid between ? K ÷ ? K

        % from "da An Innovative Strategy for Paraffin-based Fuels Reinforcement:
        % Part I, Mechanical and Pre-Burning Characterization"
        rho = 924;
    case 5
        %% high density polyethylene (HDPE)
        % valid between ? K ÷ ? K

        % from Alessandro Englaro's master thesis
        if T <= 407
            rho = 950;
        else
            rho = 1010 - 0.56*T;
        end
    case 6
        %% aluminum alloy AA 6082 - T6
        % valid between ? K ÷ ? K

        % from https://www.makeitfrom.com/
        rho = 2700;
    case 7
        %% aluminum alloy AA 6061 - T6
        % valid between ? K ÷ ? K

        % from https://www.makeitfrom.com/
        rho = 2700;
end
end

