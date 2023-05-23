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
  This code implements the mathematical model for the emptying of a tank
  full of N2O from underneath. This indeed implements both the liquid
  discharge and the gas one (tailoff) 

Changelog:
  > version: 1.0 - 04/12/2022 - Alessandro Rampazzo
    - added heading, description, and changelog
    - refactoring
    - polishing some calculations
    - joining this function with the old Modello_Svuotamento_Bifase

  > version: 0.1 - 15/03/2022 - Luca Brovedani, Alberto Scomparin
    - created
%}

function [m_out,p_out,rho_out,T_out,x_out] = tankFunction(T0,m0_ox,V_tank,k)
%% Funzioni Proprietà N2O e altre Funzioni utili

p_func = @(T) tempPressure_N2O(T,"t_p");
rhol_func = @(T) densityl_N2O(T);
rhov_func = @(T) densityv_N2O(T);
sl_func = @(T) entropyl_N2O(T);
sv_func = @(T) entropyv_N2O(T);

%% Calcolo condizioni iniziali

if p_func(T0) > 70e5
    warning('Pressione massima superiore ai 70 bar!')
end

rho_v = rhov_func(T0); %[kg/m^3]
rho_l = rhol_func(T0); %[kg/m^3]

x0 = (V_tank*rho_l*rho_v/m0_ox - rho_v)/(rho_l-rho_v); % Titolo di vapore iniziale

if x0 > 1 
    error('L''ossidante NON è in condizioni di saturazione!')
end

if x0 < 0
    error("non c'è abbastanza volume per alloggiare la massa richiesta")
end


%% Simulazione svuotamento ossidante liquido
% [p,T,m,x] = Modello_Svuotamento_Bifase(T0,f); % Script tankSatBiphase riadattato

t_n = 1; % tempo di estrazione
dt = 5e-4;
n = t_n/dt;

T(1) = T0;
m(1) = m0_ox; % massa iniziale ossidante
dmdt(1) = m(1)/t_n;

p(1) = p_func(T(1));

vl_f = @(T) 1/rhol_func(T);
vv_f = @(T) 1/rhov_func(T);

x_f = @(v, T) (v - vl_f(T))/(vv_f(T) - vl_f(T));
s_f = @(x, T) (1 - x) * sl_func(T) + x * sv_func(T);
v_f = @(x, T) (1 - x) * vl_f(T) + x * vv_f(T);

v(1) = V_tank/m(1);

x(1) = x_f(v(1), T(1));
s(1) = s_f(x(1), T(1));
S(1) = m(1) * s(1);

i = 1; eps = 1e-3; c = 0.5/eps;

while m(i) > 0 && i < 2*n

    dm = dmdt(i) * dt;
    m(i+1) = m(i) - dm;
    S(i+1) = S(i) - dm * sl_func(T(i)); % Ipotizziamo stia uscendo SOLO liquido

    s(i+1) = S(i+1)/m(i+1);
    v(i+1) = V_tank/m(i+1);

    res = 1; x_j = x(i); T_j = T(i); nn = 1;
    while res > eps && nn < 2000
        dvdx = vv_f(T_j) - vl_f(T_j);
        dvdT = c * (v_f(x_j,T_j + eps) - v_f(x_j,T_j - eps));
        dsdx = sv_func(T_j) - sl_func(T_j);
        dsdT = c * (s_f(x_j,T_j + eps) - s_f(x_j,T_j - eps));

        dv_ds = [ v(i+1) - v_f(x_j, T_j); ...
            s(i+1) - s_f(x_j, T_j) ];

        Jinv = 1/(dvdx*dsdT-dsdx*dvdT)* ...
            [dsdT, -dvdT; ...
            -dsdx, dvdx];

        dx_dT = Jinv * dv_ds;

        dx = dx_dT(1);
        dT = dx_dT(2);

        res = sqrt((dT/T_j)^2+(dx/x_j)^2);


        if x_j+dx > 1 %Riscalo dx per stare sotto a x=1
            dx_dT = 0.9*((1-(x_j-dx))/dx)*dx_dT;
            dx = dx_dT(1);
            dT = dx_dT(2);
        end

        x_j = x_j + dx;
        T_j = T_j + dT;

        if abs(1-x_j) < eps && x_j < 1
            break
        end

        if nn == 1999
            fprintf('x_j = %f\n',x_j); %a cosa serve %f??
        end

        nn = nn + 1;
    end

    x(i+1) = x_j;
    T(i+1) = T_j;

    p(i+1) = p_func(T(i+1));
    dmdt(i+1) = p(i+1)/p(1) * dmdt(1);
    i = i+1;

    if abs(1-x_j) < eps && x_j < 1 %Inserisco lo stesso if in modo da uscire da entrambi i cicli
        break
    end


end

%% Tailoff

m2 = linspace(m(end),0);
% svuotamento adiabatico del serbatoio
T_tailoff = T(end) * (m2/m2(1)).^(k-1);
p_tailoff = p(end) * (m2/m2(1)).^(k);
x_tailoff = ones(size(m2));

%% Unione

T_cutoff = 183; % [K] fusion point of N2O (N2O solidifies at temperatures that are under this one)
I = find(T_tailoff > T_cutoff);

T_out = [T,T_tailoff(I(2:end))];
m_out = [m,m2(I(2:end))];
p_out = [p,p_tailoff(I(2:end))];
x_out = [x,x_tailoff(I(2:end))];
rho_out = zeros(size(x_out));
for i = 1:length(x_out)
    if x_out(i) == 1
        rho_out(i) = rhov_func(T_out(i));
    else
        rho_out(i) = rhol_func(T_out(i));
    end
end

end