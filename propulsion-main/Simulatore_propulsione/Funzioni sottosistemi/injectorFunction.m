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
  This code implements the different models with whom we can simulate the
  oxidant mass flow and the discharged coefficient Cd.

Changelog:
  > Version: 1.6 - 25/03/2023 - Giulia Quinzi & Riccardo Dei Rossi
    - aggiornamento della pressione a monte dell'iniettore (p_tank ->
    p_tank - p_loss)
  > Version: 1.5 - 14/03/2023 - Alessandro Rampazzo & Giulia Quinzi & Riccardo Dei Rossi & Luca
    Dituri
    - abbiamo creato tre funzioni distinte: injectorFunction per il
    calcolo della potata massica di ossidante, P_LOSS per il calcolo delle perdite
    fluidiche e fluidcFunction come funzione di collegamento

  > Version: 1.4.2 - 04/03/2023 - Alessandro Rampazzo & Giulia Quinzi & Riccardo Dei Rossi & Luca
    Dituri
    - implemented equation for c_d

  > Version: 1.4.1 - 04/03/2023 - Alessandro Rampazzo & Giulia Quinzi & Riccardo Dei Rossi & Luca
    Dituri
    - modified Omega to make it work 

  > Version: 1.3 - 01/03/2023 - Alessandro Rampazzo & Giulia Quinzi & Riccardo Dei Rossi & Luca
    Dituri
    - modified Omega to make it work 

  > version: 1.2 - 25/02/2023 - Giulia Quinzi & Riccardo Dei Rossi & Luca
    Dituri
    - added p_loss: (not fully functioning) function for calculating the
    fluidic pressure losses

  > version: 1.1 - 23/02/2023 - Giulia Quinzi & Riccardo Dei Rossi
    - added Omega, Omega_mod (not fully functioning)

  > version: 1.0 - 27/11/2022 - Alessandro Rampazzo
    - added heading, description, and changelog

  > version: 0.1 - 15/03/2022 - Alessandro Rampazzo
    - created
%}

function [mDot_ox,c_d,v,m_vap,m_liq] = injectorFunction(p_cc,p_tank,m_tank,p_loss,rocket)

% DICHIARO LE VARIABILI LOCALI IN FUNZIONE DI m_tank
rho_out = rocket.tank.rho_out_func(m_tank);
x_tank = rocket.tank.x_func(m_tank);
t_tank = rocket.tank.T_func(m_tank);
p_tank = rocket.tank.p_func(m_tank);

% MODIFICHE ALL'ISOENTALPICA
% utilizzato per risolvere il problema del raccordo
xx = (x_tank - 0.95)/(1-0.95);                % non posso prendere valori maggiori di 0.95 (ho provato anche con 0.9501!!!)

% SWITCH CASE

% ELENCO DEI CASE
switch rocket.inj.model

    case "SPI"
        % Bernulli's law for injector plates
        mDot_ox = SPI(p_cc,p_tank,rho_out);

    case "SPI-ISOH0"
        % Bernulli's law for injector plates
        if x_tank < 0.95
            mDot_ox = SPI(p_cc,p_tank,rho_out);
        elseif x_tank >= 0.95 && x_tank < 1
            mDot_ox = SPI(p_cc,p_tank,rho_out)*(1-xx) + ISOH0(p_cc,p_tank)*(xx);
        else
            mDot_ox = ISOH0(p_cc,p_tank);
        end

    case "DYER-ISOH0"
        % Metodo di Dyer
        % single phase incompressible model + Two phase homogeneous equilibrium model
        if x_tank < 0.95
            psat = NFP("N2O","psat_t",t_tank)*1e5;
            if psat > p_cc
                k = sqrt(((p_tank - p_cc - p_loss))/(psat - p_cc)); % parametro di non equilibrio
                mDot_ox = (k/(k+1))*SPI(p_cc,p_tank,rho_out) + (1/(k+1))*HEM(p_cc,p_tank);
            else
                error("psat < p_cc !!!")      % se psat < p_cc si chiude la valvola e si blocca tutto ?????
            end
            % raccordo per l'integratore
        elseif x_tank >=0.95 && x_tank < 1
            psat = NFP("N2O","psat_t",t_tank)*1e5;
            if psat > p_cc
                k = sqrt(((p_tank - p_cc - p_loss))/(psat - p_cc)); % parametro di non equilibrio
                mDot_ox = (k/(k+1))*SPI(p_cc,p_tank,rho_out) + (1/(k+1))*HEM(p_cc,p_tank);
                mDot_ox = mDot_ox * (1-xx) + ISOH0(p_cc,p_tank)*(xx);
            else
                error("psat < p_cc nel raccordo !!!")
            end
        else
            mDot_ox = ISOH0(p_cc,p_tank);
        end

    case "OMEGA-ISOH0"
        if x_tank < 0.95
            mDot_ox = OMEGA(p_cc,p_tank);
        elseif x_tank >= 0.95 && x_tank < 1
            mDot_ox = OMEGA(p_cc,p_tank) * (1 - xx) + ISOH0(p_cc,p_tank) *(xx);
        else
            mDot_ox = ISOH0(p_cc,p_tank);
        end

end

% Calcolo del Cd 
c_d = C_D(rho_out,p_cc,p_tank);

% Calcolo della massa di liquido e di vapore nel tank
rhoL_tank = (NFP("N2O","vl_p",p_tank/1e5))^-1;
rhoG_tank = (NFP("N2O","vv_p",p_tank/1e5))^-1;
m_vap = x_tank*m_tank;
m_liq = (1-x_tank)*m_vap/x_tank;


%_______________________________________ FUNCTION ________________________________________________________________________

% ---------> SPI
function [mDot_ox] = SPI(p_cc,p_tank,rho_out)
    mDot_ox = C_D(rho_out,p_cc,p_tank) * rocket.inj.A_holes * sqrt(2*rho_out*(p_tank - p_cc - p_loss));
end

% ---------> HEM
function [mDot_hem_ox] = HEM(p_cc,p_tank)
    h1 = NFP("N2O","h_xp",x_tank,p_tank/1e5);
    h2 = NFP("N2O","hv_p",p_cc/1e5);
    Delta_h = h1 - h2;
    mDot_hem_ox = C_D(rho_out,p_cc,p_tank)*rho_out*rocket.inj.A_holes*sqrt(2*abs(Delta_h));
end

% ---------> ISOH0
function [mDot_isoh] = ISOH0(p_cc,p_tank)
    A_t = rocket.inj.A_holes;
    gamma = 1.33;
    R = 8314.5/44;
    M_e = sqrt(2/(gamma-1)*(((p_tank-p_loss)/p_cc)^((gamma-1)/gamma)-1));
    mDot_isoh = A_t * sqrt(gamma/R) * (p_tank-p_loss)/sqrt(t_tank) * M_e/(1+(gamma-1)/2*M_e^2)^((gamma+1)/(2*(gamma-1)));
end

% ----------> OMEGA 
function [mDot_omega] = OMEGA(p_cc,p_tank)
    psat = NFP("N2O","psat_t",t_tank)*1e5;
    tsat = NFP("N2O","tsat_p",p_tank/1e5);
    ni_tank = (NFP("N2O","v_px",p_tank/1e5,x_tank));
    niL_tank = (NFP("N2O","vl_p",p_tank/1e5));         % volume specifico del liquido nel serbatoio in condizioni di saturazione
    niG_tank = (NFP("N2O","vv_p",p_tank/1e5));         % volume specifico del vapore nel serbatoio in condizioni di saturazione
    niGL_tank = -niL_tank+niG_tank;
    hL_tank = NFP("N2O","hl_p",p_tank/1e5)*1000; % [kJ/kgK] -> [J/kgK]
    hG_tank = NFP("N2O","hv_p",p_tank/1e5)*1000; % [kJ/kgK] -> [J/kgK]
    hGL_tank = -hL_tank+hG_tank;
    cl = NFP("N2O","cpl_p",p_tank/1e5) * 1000;   % [kJ/kgK] -> [J/kgK]
    if abs(psat-p_tank)/psat < 1e-4  %SATURATED (èil nostro caso) or TWO PHASE INLET CONDITION 
        omega = 1*niGL_tank/ni_tank + cl*t_tank*(p_tank-p_loss)/ni_tank*(niGL_tank/hGL_tank)^2; % parametro omega
        hasToBeZero = @(eta_critic) eta_critic.^2 + (omega.^2-2*omega).*(1-eta_critic).^2+2*omega.^2.*(log(eta_critic)+(1-eta_critic)); % per poter calcolare eta_critic
        eta_c = fzero(hasToBeZero,[1e-6,1]);
        pc = eta_c*p_tank;
        if p_cc > pc % no  choking
            p = p_cc;
        else
            p = pc; % choking
        end
        mDot_omega = (sqrt(((p_tank-p_loss)/ni_tank)) * sqrt(-2 * (omega * log(p/(p_tank-p_loss)) + (omega-1)*(1 - (p/(p_tank-p_loss))))))/(omega * (((p_tank-p_loss)/p) - 1) + 1)* rocket.inj.A_holes;

    else  % INITIALLY SUBCOOLED LIQUID INLET CONDITION (il ciclo nel nostro caso non dovrebbe mai entrare qui)
        omega_sat = cl*t_tank*(p_tank-p_loss)/ni_tank*(niGL_tank/hGL_tank)^2;
        eta_st = 2*omega_sat/(1+2*omega_sat);
        eta_sat = psat/(p_tank-p_loss);

        if psat > eta_st*p_tank % Low subcooling region
            hasToBeZero = @(eta_critic) (omega_sat+(1/omega_sat)-2)*eta_critic^2/(2*eta_sat)-2*(omega_sat-1)*eta_critic+omega_sat*eta_sat*log(eta_critic/eta_sat)+3/2*omega_sat*eta_sat-1;
            eta_c = fzero(hasToBeZero,[1e-6,1]);
            pc = eta_c*p_tank;
            if p_cc > pc % no  choking
                p = p_cc;
            else
                p = pc; % choking
            end
            eta = p/p_tank;
            mDot_omega = (sqrt((p_tank-p_loss)/niL_tank)*sqrt(2*(1-eta_sat)+2*(omega_sat*eta_sat*log(eta_sat/eta)-(omega_sat-1)*(eta_sat-eta))))/(omega_sat*(eta_sat/eta-1)+1) * rocket.inj.A_holes;

        else % High subcooling region teoricamente non dovrebbe entrare qui
            mDot_omega = sqrt(rho_out*(p_tank-p_loss))*sqrt(2*rho_out*(p_tank-p_loss-psat)) * rocket.inj.A_holes;
        end
    end
    end

% ----------> Cd
    function c_d = C_D(rho_out,p_cc,p_tank) % !!!!!! da vedere perché non è una buona approssimazione (da vedere bene con i test)
    %calcolo della velocità
    K = 1.7;                % head loss coefficient (per bordi dritti)
    v = sqrt((2*(p_tank - p_cc - p_loss))/(K*rho_out));

    %calcolo del numero di Reynolds(v)
    d_fori = 1.5e-3;
    d_piastra = 30e-3;
    mu = NFP("N2O","mul_p",p_tank/1e5);
    Re = (rho_out * v * d_fori)/mu;

    %calcolo del c_d
    z = 0.75;
    beta = d_fori/d_piastra;
    b = 91.71 * beta^2.5;
    c_d_infty = 0.5959 + 0.0312 * beta^2.1 - 0.184 * beta^6;
    c_d = c_d_infty + b/Re^z;
end
end