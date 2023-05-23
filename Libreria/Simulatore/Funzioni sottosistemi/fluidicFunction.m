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
  This code implements the iteration process between injectorFunction and
  P_LOSS. It also adds some elements to the "out" structure
%}

function out = fluidicFunction(p_cc,p_tank,m_tank,rocket)

    iter = 0;
    itMax = 100;
    tol = 1e-6;
    err = 1e50;
    p_loss_old = 0;
    p_loss = p_loss_old;
    mDot_ox = 0;
    while err>tol && iter<itMax
        iter = iter + 1;
        mDot_ox_old = mDot_ox;
        [mDot_ox,c_d,v,m_vap,m_liq] = injectorFunction(p_cc,p_tank,m_tank,p_loss,rocket);
        p_loss_old = p_loss;
        p_loss = P_LOSS(m_tank,mDot_ox,rocket);
        err = max(abs((p_loss_old - p_loss)/(p_loss)),abs((mDot_ox_old - mDot_ox)/(mDot_ox)));
    end

    % Valori calcolati
    out.mDot_ox = mDot_ox;
    out.p_loss = p_loss;
    out.c_d = c_d;
    out.v = v;
    out.m_vap = m_vap;
    out.m_liq = m_liq;
   
end