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
  This code simulates the whole propulsion system of a hybrid rocket based
  on N2O and paraffin as propellants. It is equivalent to a static fire
  test. It solves the following set of differential equations:

    { d(m_tank)/dt = -mDot_ox
    { d(r_port)/dt = rDot = a * G^n
    { d(r_throat)/dt = <throat erosion model>
    { d(p_cc)/dt = R*T_cc/V_cc * (mDot_ox + mDot_f - p_cc/(R*T_cc)*d(V_cc)/dt - <nozzle func>)
    { d(V_cc)/dt = 2*pi*r_port*L_grain*rDot

  Where G = mDot_ox/(pi*r_port^2), with m_tank the remaining oxidizer mass
  in the tank, r_port the port radius, r_t the throat radius, L_grain the
  grain length, V_cc the combustion chamber volume, T_cc the combustion
  chamber mean temperature and R the specific gas constant. The nozzle
  function is most of the times the following expression: -p_cc*A_t/c_star


Changelog:
  > version: 1.5 - 14/03/2023 - Alessandro Rampazzo & Giulia Quinzi & Riccardo Dei Rossi & Luca
    % Dituri
    - added elements to the out structure

  > version: 1.4 - 29/01/2023 - Alessandro Rampazzo
    - added printing variables while itegrating
    - updated model in the header

  > version: 1.3 - 27/01/2023 - Alessandro Rampazzo
    - added x_tank as an output variable
    - some bugfix and refactoring
    - changed problem to account for unsteady effects in combustion
    chamber, added p_cc and V_cc integration

  > version: 1.2 - 16/12/2022 - Alessandro Rampazzo
    - fixed solveOde and tankEmpty event when changed method to input
    options in solveOde. removed also dt and runEvent = "before" since it's
    deprecated

  > version: 1.1 - 04/12/2022 - Alessandro Rampazzo
    - added a new combustionFunction that takes the O/F and outputs the gas
    properties in a struct called gasProp with fields T, k and R
    - modified the set of equation for p_cc to account for the new equation
    - added the deltaP_inj field in out struct
    - changed intergration method to 'EXP', being better for big and
    complex systems

  > version: 1.0 - 27/11/2022 - Alessandro Rampazzo
    - added heading, description, and changelog

  > version: 0.1 - 15/03/2022 - Alessandro Rampazzo
    - created
%}

function out = propulsionSimulator2(rocket,env,tSpan)

    iter = 0;
    % Initial conditions [m_ox_tank, port_radius, throat_radius]
    Y0 = [rocket.tank.m_vec(1), rocket.cc.D_port0/2,...
        rocket.nozzle.D_throat0/2, 3e5, rocket.cc.V_tot0];

    % terminal event for integration
    breakEvent = @(t,Y) tankEmpty(t,Y,rocket,env);
    
    % solve the problem (for big and complicated systems explicit method is
    % better)

    [t,y] = solveOde(@(t,Y) propulsiveSystem(t,Y,rocket,env),'RK45_A',tSpan,Y0,...
        "breakevent", breakEvent,'abstol',inf,'reltol',1e-6);
%     op.RelTol = 1e-6;
%     op.Events = breakEvent;
%     [t,y] = ode45(@(t,Y) propulsiveSystem(t,Y,rocket,env),tSpan,Y0,op);
    
    % extract the solution undersampling to reduce data weight
    Idx = 1:5:length(t);
    out.t = t(Idx);
    out.m_tank = y(Idx,1);
    out.r_port = y(Idx,2);
    out.r_throat = y(Idx,3);
    out.p_cc = y(Idx,4);
    out.V_cc = y(Idx,5);
    if out.r_port(end) > rocket.cc.D_ext
        warning("burned more fuel than there were in the combustion chamber!")
    end
    
    out.A_throat = pi * out.r_throat.^2;
    out.eps = rocket.nozzle.A_exit./out.A_throat;
    
    % initialize some vectors
    out.rDot_port = zeros(size(out.t));
    out.p_tank = zeros(size(out.t));
    out.T_tank = zeros(size(out.t));
    out.x_tank = zeros(size(out.t));
    out.mDot_ox = zeros(size(out.t));
    out.mDot_fuel = zeros(size(out.t));
    out.mDot_tot = zeros(size(out.t));
    out.S = zeros(size(out.t));
    out.Isp = zeros(size(out.t));
    out.deltaP_inj = zeros(size(out.t));
    out.rDot_throat = zeros(size(out.t));
    out.O_F = zeros(size(out.t));
    out.cond = zeros(size(out.t));
    out.rho_out = zeros(size(out.t));
    out.p_loss = zeros(size(out.t));
    out.c_d = zeros(size(out.t));
    out.v = zeros(size(out.t));
    out.m_vap = zeros(size(out.t));
    out.m_liq = zeros(size(out.t));

    % extract the other variables redoing the calculations (see
    % propulsiveSystem nested function)
    for i = 1:length(out.t)
        out.p_tank(i) = rocket.tank.p_func(out.m_tank(i));
        out.rho_out(i) = rocket.tank.rho_out_func(out.m_tank(i));
        out.T_tank(i) = rocket.tank.T_func(out.m_tank(i));
        out.x_tank(i) = rocket.tank.x_func(out.m_tank(i));
    
        out.mDot_ox(i) = rocket.inj.fluidicFunction(out.p_cc(i),out.p_tank(i),out.m_tank(i)).mDot_ox;
        out.p_loss(i) = rocket.inj.fluidicFunction(out.p_cc(i),out.p_tank(i),out.m_tank(i)).p_loss;
        out.c_d(i) = rocket.inj.fluidicFunction(out.p_cc(i),out.p_tank(i),out.m_tank(i)).c_d;
        out.v(i) = rocket.inj.fluidicFunction(out.p_cc(i),out.p_tank(i),out.m_tank(i)).v; % velocità nell'iniettore
        out.m_vap(i) = rocket.inj.fluidicFunction(out.p_cc(i),out.p_tank(i),out.m_tank(i)).m_vap;
        out.m_liq(i) = rocket.inj.fluidicFunction(out.p_cc(i),out.p_tank(i),out.m_tank(i)).m_liq;
        grainOut2 = rocket.cc.grainFunction(out.mDot_ox(i),out.r_port(i));
        out.mDot_fuel(i) = grainOut2.mDot_fuel;
        out.rDot_port(i) = grainOut2.rDot_port;
        out.O_F(i) = out.mDot_ox(i)/out.mDot_fuel(i);

        out.mDot_tot(i) = out.mDot_fuel(i) + out.mDot_ox(i);
        out.rDot_throat(i) = rocket.nozzle.nozzleErosion(out.t(i),out.p_cc(i),out.A_throat(i),out.mDot_tot(i),out.O_F(i));
        out.deltaP_inj(i) = out.p_tank(i) - out.p_cc(i) - rocket.inj.p_loss;
        out.gasProp(i) = rocket.cc.combustionFunction(out.O_F(i),out.T_tank(i),rocket.cc.T_fuel);
        nozzleOut = rocket.nozzle.nozzleFunction(out.p_cc(i),out.A_throat(i),out.eps(i),out.gasProp(i));

        out.S(i) = nozzleOut.S;
        out.Isp(i) = out.S(i)/(env.g*nozzleOut.mDot_tot);
        out.cond(i) = nozzleOut.cond;
    end

    clear iter

    %% nested functions

    function dY = propulsiveSystem(t,Y,rocket,env)
        % main ode system

        % initialize
        dY = zeros(size(Y));

        % extract the variables from Y vector
        m_tank = Y(1);
        r_port = Y(2);
        r_throat = Y(3);
        p_cc = Y(4);
        V_cc = Y(5);

        if abs(p_cc - 0.331081467706986e6) < 1e-4
            keyboard
        end

        % tank variables using interpolation functions
        p_tank = rocket.tank.p_func(m_tank);
        rho_out = rocket.tank.rho_out_func(m_tank);
        T_tank = rocket.tank.T_func(m_tank);

%         if rho_out < 200
%             disp("lol")
%         end
        
        % to find breakevent
        if isnan(p_tank)
            return
        end

        % throat area
        A_throat = pi * r_throat^2;

        % find p_cc solving the algebraic system of equation cited above
%         [p_cc,mDot_ox,mDot_fuel,rDot_port] = findPcc(p_tank,rho_out,T_tank,A_throat,r_port,rocket,env);
        mDot_ox = rocket.inj.fluidicFunction(p_cc,p_tank,m_tank).mDot_ox;

        % to prevent errors if the time step becomes too big around
        % discontinuities (like in the rho_out vector). pressure in
        % combustion chamber is greater than the tank pressure and mDot_ox
        % is complex
        if ~isreal(mDot_ox)
            return
        end
        grainOut = rocket.cc.grainFunction(mDot_ox,r_port);
        O_F = mDot_ox/grainOut.mDot_fuel;
        gasProp = rocket.cc.combustionFunction(O_F,T_tank,rocket.cc.T_fuel);
        eps = rocket.nozzle.A_exit/A_throat;

        % derivative of oxidizer mass flow rate
        dY(1) = -mDot_ox;
    
        % derivative of port radius
        dY(2) = grainOut.rDot_port;
    
        % derivative of nozzle throat radius     
        dY(3) = rocket.nozzle.nozzleErosion(t,p_cc,A_throat,grainOut.mDot_fuel + mDot_ox,O_F);

        % derivative of combustion chamber volume
        dY(5) = 2*pi*grainOut.rDot_port * r_port * rocket.cc.L_grain;

        % derivative of combustion chamber pressure
        dY(4) = gasProp.R * gasProp.T / V_cc * (mDot_ox + grainOut.mDot_fuel -...
            p_cc/(gasProp.R * gasProp.T) * dY(5) - ...
            rocket.nozzle.nozzleFunction(p_cc,A_throat,eps,gasProp).mDot_tot);
        

    end

    function [position,isterminal,direction] = tankEmpty(t,Y,rocket,env)
        % uncomment this line to print the time on the command window
        if mod(iter,200) == 0
            fprintf("t = %6.3f s | m_t = %5.3f kg | r_f = %6.3f mm | r_t = %5.3f mm " + ...
                "| p_cc = %6.3f bar | V_cc = %5.3f L\n",t,Y(1),Y(2)*1e3,Y(3)*1e3,Y(4)/1e5,Y(5)*1000);
        end
        iter = iter + 1;
        % same thing as above but...
        m_tank = Y(1);
%         r_f = Y(2);
%         r_t = Y(3);
        p_cc = Y(4);
%         V_cc = Y(5);
        
        p_tank = rocket.tank.p_func(m_tank);
%         rho_out = rocket.tank.rho_out_func(m_tank);
%         T_tank = rocket.tank.T_func(m_tank);

        if isnan(p_tank)
            position = 0; % The value that we want to be zero        
            isterminal = true;  % Halt integration ?
            direction = false;   % The zero can be approached from either direction ?
            return
        end

        % I want to block the integration if the remaining oxidizer mass is
        % out of the boundaries of the model of the tank or when the
        % combustion chamber pressure is lower than atmospheric

        % is position = 0 the event is raised
        position = (m_tank > rocket.tank.m_vec(end)+1e-6) && p_cc > 101325 + 1; % The value that we want to be zero        
        isterminal = true;  % Halt integration ?
        direction = false;   % The zero can be approached from either direction ?
    end
end