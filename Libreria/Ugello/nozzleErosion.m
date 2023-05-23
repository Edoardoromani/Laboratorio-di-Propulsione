function rDot_t = nozzleErosion(t,p_cc,r_t,mDot_tot,O_F,rocket,model)

% throat erosion model
if t >= rocket.nozzle.onsetTime
    
    A_t = pi * r_t^2; % throat area
    
    if model == model_1

        % from paper: Erosion Rate Investigation of Various Nozzle Materials in Hybrid Rocket Motors
        rDot_t = rocket.nozzle.a * O_F^rocket.nozzle.m * (mDot_tot/A_t)^rocket.nozzle.n;

    elseif model == model_2

        % from paper: Numerical and Experimental Investigation of Nozzle Thermochemical Erosion in Hybrid Rockets
        OF_paper = 3.42;
        alfa = 0.0179/OF_paper; % mm/s, with cc pressure in bar, alfa is a function of the propellenat combination and nozzle shape
        rDot_t = alfa*O_F*p_cc^0.8;

    end
    
else  
    
    rDot_t = 0;
    
end
end