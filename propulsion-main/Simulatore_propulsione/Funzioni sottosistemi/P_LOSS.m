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
  This code implements the fluidic losses with the oxidizer flow max
  updating at each iteration (thanks to the fluidicFunction)
%}

function p_loss = P_LOSS(m_tank,mDot_ox,rocket)
    
    modello = "VAR";
    % lista modelli: COST - VAR
    switch modello
        case "COST"
            p_loss = 2.5e5;
        case "VAR"
            % geometria della linea fluidica
            d1 = 0.0306;                                   %[m]   diametro linea fluidica (di ingresso del convergente)
            d2 = 0.0165;                                   %[m]   diametro utile della piastra di iniezione
            l_conv = 0.05;                                 %[m]   lunghezza tratto convergente
            theta_mez = atan((d1-d2)/(2*l_conv));          %[rad] angolo di apertura tratto convergente
            l_tot = 0.375;                                 %[m]   lunghezza totale della linea fluidica    
    
            % definizione dei coefficienti di perdita di carico (già calcolati a
            % priori)
            lambda_imb = 0.5;                                   % lambda imbocco
            lambda_distr = 0.073;                               % lamdba distribuito
            lambda_gom = 2;                                     % lambda gomito
            lambda_conv = 0.8*sin(theta_mez)*(1-(d2/d1)^2)^2;   % lambda convergente
            rho_out = rocket.tank.rho_out_func(m_tank);
            p_loss = rho_out*((mDot_ox/ (rho_out * pi*d1^2/4))^2)/2 * (lambda_imb + lambda_distr * (l_tot/d1) + 2*lambda_gom + lambda_conv);
    
    end
end