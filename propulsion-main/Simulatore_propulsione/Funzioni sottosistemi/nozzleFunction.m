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
  This code solves an isentropic chemically-frozen nozzle of the rocket
  calculating the thrust and its condition. It implements every operative 
  case except for the supersonic inlet.

Changelog:
  > version: 1.1 - 04/12/2022 - Alessandro Rampazzo
    - changed the inputs to account for a separated combustion function
    which will perform accurate combustion calculations

  > version: 1.0 - 27/11/2022 - Alessandro Rampazzo
    - added heading, description, and changelog

  > version: 0.1 - 15/03/2022 - Alessandro Rampazzo
    - created
%}

function out = nozzleFunction(p_cc,A_t,eps,gasProp,rocket,env)

% p_cc = combustion chamber pressure [Pa]

% combustion
% total temperature [K]
T_0 = gasProp.T;
% adiabatic expansion coefficient
gamma = gasProp.gamma;
% gas constant [J/kgK]
R = gasProp.R;
% atmospheric pressure [Pa]
p_a = env.p_a;

% exit area
A_e = A_t*eps;

% function correlating expansion ratio and pressure ratio
hasToBeZero = @(p_e) ((gamma+1)/2)^(1/(gamma-1)) .* (p_e./p_cc).^(1/gamma) .*...
    sqrt((gamma+1)/(gamma-1).*(1-(p_e./p_cc).^((gamma-1)/gamma))) - A_t/A_e;

% calculate the critical pressure and find the supersonic outlet pressure
% and the subsonic one
p_crit = p_cc / ((gamma+1)/2)^(gamma/(gamma-1));
p_e_sup = fzero(hasToBeZero,[1e-6,p_crit-1e-6]);
p_e_sub = fzero(hasToBeZero,[p_crit+1e-6,p_cc-1e-6]);

% default shock position (no shock -> A_shock == -1)
% A_shock = -1;


if  p_a < p_cc && p_a >= p_e_sub
    %% subsonic
    cond = "subsonic";

    % adapted exit
    p_e = p_a;
    % exit Mach
    M_e = sqrt(2/(gamma-1)*((p_cc/p_e)^((gamma-1)/gamma)-1));
    % mass flow rate calculated at the exit 
    mDot = A_e * sqrt(gamma/R) * p_cc/sqrt(T_0) * M_e/(1+(gamma-1)/2*M_e^2)^((gamma+1)/(2*(gamma-1)));
    % exit temperature
    T_e = T_0 * (1 + (gamma-1)/2*M_e^2);
    % function that has to be zero for the throat mach
%     hasToBeZero = @(M) sqrt(k/(R*T_0)) * p_cc * M./(1+(k-1)/2*M^2).^((k+1)/(2*(k-1))) - mDot/A_t;
%     M_t = fzero(hasToBeZero,[0,1]);

elseif p_a < p_e_sub && p_a > p_e_sup
    %% overexpanded
%     M_t = 1;

    % some gasdynamics
    mDot = A_t * sqrt(gamma/R) * p_cc/sqrt(T_0) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
    M_e = sqrt(2/(gamma-1)*((p_cc/p_e_sup)^((gamma-1)/gamma)-1));

    % exit pressure if a normal shock is present
    p_e_shock = p_e_sup * (1 + 2*gamma/(gamma+1)*(M_e^2-1));

    if p_a <= p_e_shock
        % no shocks inside the nozzle
        cond = "overexpanded without shock";

        p_e = p_e_sup;
        T_e = T_0 / (1 + (gamma-1)/2*M_e^2);
        
    else
        % shock in the nozzle
        cond = "overexpanded with shock";

        p_e = p_a;
        hasToBeZero = @(M) sqrt(gamma/(R*T_0)) * p_e * M * sqrt(1+(gamma-1)/2*M) - mDot/A_e;
        M_e = fzero(hasToBeZero,[0,1]);
        T_e = T_0 / (1 + (gamma-1)/2*M_e^2);
        % shock formulas
%         p_0_2 = p_e * (1+(k-1)/2*M_e^2)^(k/(k-1));
%         M_2 = @(M_1) sqrt((1+(k-1)/2*M_1^2) / (k*M_1^2-(k-1)/2));
%         hasToBeZero = @(M_1) (1+(k-1)/2*M_1^2)^(-k/(k-1)) *...
%             (1+(k-1)/2*M_2(M_1)^2)^(k/(k-1)) *...
%             (1+2*k/(k+1)*(M_1^2-1))...
%             - p_0_2/p_cc;
        % find shock Mach
%         M_1 = fzero(hasToBeZero,[1,100]);
        % find area at which normal shock occurs
%         A_shock = A_t/M_1 * sqrt(2/(k+1)*(1+(k-1)/2*M_1^2))^((k+1)/(k-1));      
    end
    
elseif p_a <= p_e_sup
    %% underexpanded or perfectly adapted
%     M_t = 1;
    mDot = A_t * sqrt(gamma/R) * p_cc/sqrt(T_0) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
    p_e = p_e_sup;
    M_e = sqrt(2/(gamma-1)*((p_cc/p_e)^((gamma-1)/gamma)-1));
    T_e = T_0 / (1 + (gamma-1)/2*M_e^2);

    if abs(p_a - p_e_sup)/p_e_sup < 10^-4
        cond = "perfectly adapted";
    else
        cond = "underexpanded";
    end
elseif 0 < p_cc && p_cc <= p_a
    % no flow
    cond = "no flow";
%     M_t = 0;
    mDot = 0;
    p_e = 0;
    M_e = 0;
    T_e = T_0;
else
    error("error while solving nozzle")
end

% calculate thrust (S)
if cond == "overexpanded with shock" || cond == "subsonic"
    v = sqrt(gamma*R*T_e) * M_e;
    S = rocket.nozzle.eta_Cf*(mDot*v +  A_e * (p_e - p_a));
else
    Gamma = gamma*(2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
    Cf = Gamma * sqrt(2/(gamma-1) * (1-(p_e/p_cc)^((gamma-1)/gamma))) + A_e/A_t*(p_e - p_a)/p_cc;
    S = rocket.nozzle.eta_Cf * rocket.cc.cstar * Cf * mDot;
end

% packing everything in out struct
out.S = S;
out.mDot_tot = mDot;

% assigning a code to every state
switch cond
    case "no flow"
        out.cond = 0;
    case "subsonic"
        out.cond = 1;
    case "overexpanded with shock"
        out.cond = 2;
    case "overexpanded without shock"
        out.cond = 3;
    case "perfectly adapted"
        out.cond = 4;
    case "underexpanded"
        out.cond = 5;
end
end