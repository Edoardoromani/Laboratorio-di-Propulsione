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
  This function returns the thermal conductivity in [W/kgm^2] given the
  material ID and the temperature

Changelog:
  > version: 1.0 - 13/11/2022 - Alessandro Rampazzo
    - created
%}

function lambda = getLambda(materialID,T)

switch materialID
    case 0
        %% graphite
        % valid between 200 K ÷ 3000 K

        % from "A new numerical method and modified apparatus for the simultaneous
        % evaluation of thermo-physical properties above 1500 K: A case study on
        % isostatically pressed graphite", Least square interpolation
        % lambda = 24.8456 + 94.0797*(T/250).^0.4280 .* exp(-0.3983*T/250);

        % from "Computational Fluid-dynamic Simulations of Hybrid Rocket
        % Internal Flow Including Discharge Nozzle"
        lambda = 104;
    case 1
        %% steel
        % valid between ? K ÷ ? K

        % from https://it.wikipedia.org/wiki/Acciaio_inossidabile
        lambda = 17;
    case 2
        %% phenolic resin Si based
        % valid between ? K ÷ ? K

        % from https://www.researchgate.net/figure/Properties-of-Si-phenolic-resin_tbl1_276225652
        lambda = 0.485;
    case 3
        %% aluminum alloy AA 7075 - T6 (Ergal)
        % valid between ? K ÷ ? K

        % from https://www.comefimetalli.it/lega7075.asp
        lambda = 155;
    case 4
        %% paraffin sasolwax 0907
        % valid between ? K ÷ ? K

        % da https://material-properties.org
        lambda = 0.2;
    case 5
        %% high density polyethylene (HDPE)
        % valid between ? K ÷ ? K

        % from Alessandro Englaro's master thesis
        if T <= 407
            lambda = 0.15 + 1.1e-4*T;
        else
            lambda = 0.08 + 4.2e-4*T;
        end
    case 6
        %% aluminum alloy AA 6082 - T6
        % valid between ? K ÷ ? K

        % from https://www.makeitfrom.com/
        lambda = 160;
    case 7
        %% aluminum alloy AA 6061 - T6
        % valid between ? K ÷ ? K

        % from https://www.makeitfrom.com/
        lambda = 170;
end
end

