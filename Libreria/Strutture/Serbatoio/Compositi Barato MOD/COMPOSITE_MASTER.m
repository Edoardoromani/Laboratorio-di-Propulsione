%   COMPOSITE MASTER
%   ----------------
%   Federico Bortolato 2023
%   -----------------------
%   One code to rule them all
%   -------------------------
%   Input : fiber properties, 
%           matrix properties
%           ply stacking
%           loads data
%   Output: matrice di rigidezza
%           tensioni e deformazioni
%           Strenght Ratio - modalità First Ply Failure - Tsai-Hill failure criterion
%           Carichi massimi applicabili con la stessa modalità di carico 
%   --------------------------------------------------------------------
%   Updated from previous code of Francesco Barato 2012
%   ---------------------------------------------------

%%  INITIALIZATION

clear all
close all
% clc

format bank
format compact
 
%% DATA - input - stacking of laminas -> creation of laminate

n = 18;                                             %   number of plies       
t_ply= 0.153*ones(1,n);                            %   plies' thickness [mm]
% teta_ply = [90 90 90 90 90 0 0 0 90 90 0 0 0 90 90 90 90 90]
teta_ply= 52.55*(-ones(1,n)).^(0:n-1);              %   angle [degrees]
material= ones(n,1);                    %   material 


%   material = 1 -> carbon fibre composite
%   material = 2 -> rohacell rigid foam

if length(t_ply) ~= n
    error("n e t_ply devono essere della stessa misura!")
elseif length(teta_ply) ~= n
    error("n e teta_ply devono essere della stessa misura!")
elseif length(material) ~= n
    error("n e material devono essere della stessa misura!")
end

%%  DATA - input - moduli of lamina - carbon/epoxy GENERAL mix

% E_f = 230e9;                %   Young2s modulus of fiber [Pa]
% G_f = 22e9;                 %   Shear Modulus of fiber [Pa]
% v_f = 0.27;                 %   Poissons's ratio of fiber [] - REFERENCED
% 
% V_f = 0.6;                   %   Fiber Volume Fraction []
% 
% E_m = 3.9e9;                %   Young s modulus of matrix [Pa]
% G_m = 1.4e9;                %   Shear Modulus of matrix [Pa]
% v_m = 0.35;                 %   Poissons's ratio of matrix []
% 
% %%  DATA - input - strengths of lamina - carbon/epoxy GENERAL mix
% 
% sigma_f_t_ult = 2067e6;           %   Ultimate tensile strength of fiber [Pa]
% tau_f_ult = 36e6;                 %   Ultimate shear strength of fiber [Pa]
% 
% sigma_m_t_ult = 72e6;             %   Ultimate tensile strength of matrix [Pa]
% sigma_m_c_ult = 102e6;            %   Ultimate compressive strength of matrix [Pa]
% tau_m_ult = 34e6;                 %   Ultimate shear strength of matrix [Pa]
% %   probably equal values tensile/compressive
% 
% %%  PROPERTIES CALCULATOR - carbon/epoxy GENERAL mix
% [U1_m(1),U2_m(1),U3_m(1),U4_m(1),sigma1t_ult(1),sigma1c_ult(1),sigma2t_ult(1),sigma2c_ult(1),tau12_ult(1)] ...
%     = COMPOSITE_PROPERTIES(E_f,G_f,v_f,V_f,E_m,G_m,v_m,sigma_f_t_ult,sigma_m_t_ult,sigma_m_c_ult,tau_m_ult,tau_f_ult);

%%   Ultimate strenghts for the 0° T300/5208

% sigma1t_ult(1)=1500*10^6;
% sigma1c_ult(1)=1500*10^6;
% sigma2t_ult(1)=40*10^6;
% sigma2c_ult(1)=246*10^6;
% tau12_ult(1)=68*10^6;
% %   reduced stiffness matrix for the 0° graphite/epoxy ply
% E1=181.0*10^9;
% E2=10.30*10^9;
% G12=7.17*10^9;
% v12=0.28;
% v21=E2*v12/E1;
% Q11=E1/(1-v12*v21); 
% Q12=v12*E2/(1-v12*v21);  
% Q22=E2/(1-v12*v21); 
% Q66=G12;
% %   four invariants
% U1_m(1) = (3*Q11+3*Q22+2*Q12+4*Q66)/8;
% U2_m(1) = (Q11-Q22)/2;
% U3_m(1) = (Q11+Q22-2*Q12-4*Q66)/8;
% U4_m(1) = (Q11+Q22+6*Q12-4*Q66)/8;

%%   Ultimate strenghts for the 0° H-IM6/epoxy
sigma1t_ult(1)=3500*10^6;
sigma1c_ult(1)=1540*10^6;
sigma2t_ult(1)=56*10^6;
sigma2c_ult(1)=150*10^6;
tau12_ult(1)=98*10^6;
%   reduced stiffness matrix for the 0° graphite/epoxy ply
E1=203*10^9;
E2=11.2*10^9;
G12=8.4*10^9;
v12=0.28;
v21=E2*v12/E1;
Q11=E1/(1-v12*v21); 
Q12=v12*E2/(1-v12*v21);  
Q22=E2/(1-v12*v21); 
Q66=G12;
%   four invariants
U1_m(1) = (3*Q11+3*Q22+2*Q12+4*Q66)/8;
U2_m(1) = (Q11-Q22)/2;
U3_m(1) = (Q11+Q22-2*Q12-4*Q66)/8;
U4_m(1) = (Q11+Q22+6*Q12-4*Q66)/8;

%%  DATA - input - moduli of lamina - rohacell 31 IG-F
% 
% E_f = 36e6;               %   Younges modulus of fiber
% G_f = 13e6;               %   Shear Modulus of fiber
% v = E_f/2/G_f-1;        %   Poissons's ratio ISOTROPE MATERIAl
% v_f = v;                %   Poissons's ratio of fiber -> around 0.38
% 
% V_f = 0;                %   Fiber Volume Fraction
% 
% E_m = 36e6;               %   Younges modulus of matrix
% G_m = 13e6;               %   Shear Modulus of matrix
% v_m = v;                %   Poissons's ratio of matrix
% 
% %%  DATA - input - strengths of lamina - rohacell 31 IG-F
% 
% sigma_f_t_ult = 1e6;        %   Ultimate tensile strength of fiber
% tau_f_ult = 0.4e6;            %   Ultimate shear strength of fiber
% 
% sigma_m_t_ult = 1.0e6;        %   Ultimate tensile strength of matrix
% sigma_m_c_ult = 0.4e6;        %   Ultimate compressive strength of matrix
% tau_m_ult = 0.4e6;            %   Ultimate shear strength of matrix
% 
% %%  PROPERTIES CALCULATOR - rohacell 31 IG-F
% [U1_m(2),U2_m(2),U3_m(2),U4_m(2),sigma1t_ult(2),sigma1c_ult(2),sigma2t_ult(2),sigma2c_ult(2),tau12_ult(2)] ...
%     = COMPOSITE_PROPERTIES(E_f,G_f,v_f,V_f,E_m,G_m,v_m,sigma_f_t_ult,sigma_m_t_ult,sigma_m_c_ult,tau_m_ult,tau_f_ult);

%% DATA - input - Loads acting on laminate

%   Formulas for pressure vessel

p = 6e6;        %   [Pa]    Pressure
FoS = 3;        %   [  ]    Safety Factor
p = p*FoS;      %   [Pa]    New Pressure
r = 0.075-sum(t_ply)/1000;      %   [m]     Internal Radius
Nhoop = p*r;    %   [N/m]   Hoop Stresses
Naxx = p*r/2;   %   [N/m]   Axial Stresses
Nextra = 1.123337859934738e+05;   %   [N/m]   Axial Stresses

Nx = Naxx;  %   [N/m]   Laminate Stress direction x (0°)
Nx = Nx + Nextra;
Ny = Nhoop;             %   [N/m]   Laminate Stress direction y (90°)

Loads = [ Nx Ny 0 0 0 0 ]; %  [N/m]   Loads vector

% fprintf("\n\nNx = %8.4f [N/m]    Ny = %8.4f [N/m]\n\n",Naxx,Nhoop)

%%  PROPERTIES CALCULATOR
[Nxmax,Nymay] = COMPOSITE_CALCULATOR(U1_m,U2_m,U3_m,U4_m,sigma1t_ult,sigma1c_ult,sigma2t_ult,sigma2c_ult,tau12_ult,n,t_ply,teta_ply,material,Loads);
pmax = Nymay/r;
%%  SCONTRINO PRESSIONI
fprintf('----------------------\n')
fprintf('~~~~~~~~~~~~~~~~~~~~~~\n')
fprintf('RISULTATI DEMMERDA\n')
fprintf('~~~~~~~~~~~~~~~~~~~~~~\n')
fprintf('----------------------\n')
fprintf('Design Pressure:\n')
fprintf('Preal = %4.1f [Bar]\n',p/1e5)
fprintf('----------------------\n')
fprintf('Number of plies:\n')
fprintf('n     = %4.0f\n',n)
fprintf('----------------------\n')
fprintf('Max Pressure Sustainable:\n')
fprintf('Pmax  = %4.1f [Bar]\n',pmax/1e5)
fprintf('----------------------\n')
fprintf("Tot thickness = %4.1f\n",sum(t_ply))
fprintf('----------------------\n')
