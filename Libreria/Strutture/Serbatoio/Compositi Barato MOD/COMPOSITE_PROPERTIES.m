%   PROPERTIES CALCULATOR
% -----------------------
%   Adapted from : 
%   AUTAR K. KAW, Mechanics of composite materials, chapter 3
%   ---------------------------------------------------------
%   Code retrieves ply parameters from matrix 
%   and fiber individual properties 

function [U1_m,U2_m,U3_m,U4_m,sigma1t_ult,sigma1c_ult,sigma2t_ult,sigma2c_ult,tau12_ult] = ...
    COMPOSITE_PROPERTIES(E_f,G_f,v_f,V_f,E_m,G_m,v_m,sigma_f_t_ult,sigma_m_t_ult,sigma_m_c_ult,tau_m_ult,tau_f_ult)

%%  Evaluation of the four Elastic Moduli

E1 = E_f*V_f + E_m*(1-V_f);
%   Longitudinal Young's Modulus
E2 = E_f*E_m/(V_f*E_m+(1-V_f)*E_f);
%   Transverse Young's Modulus
v12 = v_f*V_f + v_m*(1-V_f);
%   Major Poisson’s Ratio
G12 = G_f*G_m/(V_f*G_m+(1-V_f)*G_f);
%    In-Plane Shear Modulus

%%  Calculation of reduced stiffness matrix for the 0° ply

v21 = E2*v12/E1;    
%   Minor Poisson's Ratio
Q11 = E1/(1-v12*v21); 
Q12 = v12*E2/(1-v12*v21);  
Q22 = E2/(1-v12*v21); 
Q66 = G12;

%%  Calculation of four invariants

U1_m = (3*Q11+3*Q22+2*Q12+4*Q66)/8;
U2_m = (Q11-Q22)/2;
U3_m = (Q11+Q22-2*Q12-4*Q66)/8;
U4_m = (Q11+Q22+6*Q12-4*Q66)/8;

%%  Longitudinal Tensile Strength - CHAPTER 3.4.1
eps_f_t_ult = sigma_f_t_ult/E_f;       
%   ε_f_t_ult : Ultimate failure tensile strain of the fiber
eps_m_t_ult = sigma_m_t_ult/E_m;       
%   ε_m_t_ult : Ultimate failure tensile strain of the matrix

sigma1t_ult = sigma_f_t_ult*V_f + eps_f_t_ult*E_m*(1-V_f);  
%   Longitudinal Tensile Strength

%%  Longitudinal Compressive Strength - CHAPTER 3.4.2

d_over_s = (4*V_f/pi)^(1/2);    
%   Ratio of the diameter, d, to fiber spacing, s

%   THREE MODES: choose one
MODE = 1;  %   2,3
%   criteria to choose highest, middle or lowest strength from calculations
%-------------------------------------------------------------------
switch MODE
    case 1
        %   INIT MODE 1: --------------------   HIGHEST
        %   Fracture of matrix and/or 
        %   fiber–matrix bond due to tensile strains
        %   in the matrix and/or bond
        
        % eps_2_t_ult = eps_m_t_ult*(1-V_f^(1/3));              
        %   Ultimate transverse tensile strain Versione 1: empirical formula
        eps_2_t_ult = eps_m_t_ult*(d_over_s*(E_m/E_f-1)+1); 
        %   Ultimate transverse tensile strain Versione 2: mechanics of material formula 
        
        sigma1c_ult = E1*eps_2_t_ult/v12;   
        %   Longitudinal Compressive strength
        %   END MODE 1: ----------------------------------------------------
        %-------------------------------------------------------------------
    case 2 
        %   INIT MODE 2: --------------------   MID
        %   Microbuckling of fibers in shear or extensional mode
        %   Usually S1c>S2c; S2c only in low fiber volume fraction composites
        
        S1c = 2*(V_f+(1-V_f)*E_m/E_f)*sqrt(V_f*E_m*E_f/3/(1-V_f));
        %   Extensional mode buckling stress
        S2c = G_m*(1-V_f);
        %   Shear mode buckling stress
        sigma1c_ult = min(S1c,S2c); 
        %   Longitudinal Compressive strength
        %   END MODE 2: ----------------------------------------------------
        %-------------------------------------------------------------------
    case 3
        %   INIT MODE 2: --------------------  LOW 
        %   Shear failure of fibers
        
        sigma1c_ult = 2*(tau_f_ult*V_f + tau_m_ult*(1-V_f)); 
        %   Longitudinal Compressive strength
        %   END MODE 3: ----------------------------------------------------
        %-------------------------------------------------------------------
end

%%  Transverse Tensile Strength - CHAPTER 3.4.3
eps_2_t_ult = eps_m_t_ult*(d_over_s*(E_m/E_f-1)+1); 
%   Ultimate transverse tensile strain Versione 2: mechanics of material formula 
sigma2t_ult = E2*eps_2_t_ult;    
%   Transverse tensile strength

%%  Transverse Compressive Strength - CHAPTER 3.4.4
eps_m_c_ult = sigma_m_c_ult/E_m;       
%   ε_m_c_ult : Ultimate failure compressive strain of the matrix
eps_2_c_ult = eps_m_c_ult*(d_over_s*(E_m/E_f-1)+1); 
%   Ultimate transverse tensile strain Versione 2: mechanics of material formula 
sigma2c_ult = E2*eps_2_c_ult;  
%   Transverse compressive strength

%%  In-Plane Shear Strength - CHAPTER 3.4.5
gamma12_m_ult = tau_m_ult/G_m;
% Ultimate shearing strain of the matrix
gamma12_ult = (d_over_s*G_m/G_f + (1-d_over_s))*gamma12_m_ult;
%   In-plane shearing strains in the composite
tau12_ult = G12*gamma12_ult;
%   In-plane shear strength

end
