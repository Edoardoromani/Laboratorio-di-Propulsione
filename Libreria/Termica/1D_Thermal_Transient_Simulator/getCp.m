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
  This function returns the specific heat in [J/kgK] given the material ID 
  and the temperature

Changelog:
  > version: 1.0 - 13/11/2022 - Alessandro Rampazzo
    - created
%}

function Cp = getCp(materialID,T)

switch materialID
    case 0
        %% graphite
        % valid between 200 K ÷ 3000 K

        % from NIST - Graphite - Butland and Maddison, 1973
        % Cp = 4184*(0.538657 + 9.11129*1e-6*T - 90.2725./T - 43449.3./T.^2 + 1.59309*1e+7./T.^3 - 1.43688*1e+9./T.^4);

        % from "Computational Fluid-dynamic Simulations of Hybrid Rocket
        % Internal Flow Including Discharge Nozzle"
        Cp = 710;
    case 1
        %% steel
        % valid between ? K ÷ ? K

        % da https://it.wikipedia.org/wiki/Calore_specifico
        Cp = 502;
    case 2
        %% phenolic resin Si based
        % valid between ? K ÷ ? K

        % from https://www.researchgate.net/figure/Properties-of-Si-phenolic-resin_tbl1_276225652
        Cp = 1256;
    case 3
        %% aluminum alloy AA 7075 - T6 (Ergal)
        % valid between ? K ÷ ? K

        % from http://www.manenti.biz/it/prodotti/lastre/item/3-lastre-alluminio-7075.html
        Cp = 915;
    case 4
        %% paraffin sasolwax 0907
        % valid between ? K ÷ ? K

        % da https://material-properties.org
        Cp = 2200;
    case 5
        %% high density polyethylene (HDPE)
        % valid between ? K ÷ ? K

        % from Alessandro Englaro's master thesis
        if T <= 407
            Cp = -1040 + 9*T;
        else
            Cp = 370 + 5.1*T;
        end
    case 6
        %% aluminum alloy AA 6082 - T6
        % valid between ? K ÷ ? K

        % from https://www.makeitfrom.com/
        Cp = 900;
    case 7
        %% aluminum alloy AA 6061 - T6
        % valid between ? K ÷ ? K

        % from https://www.makeitfrom.com/
        Cp = 900;
end
end

