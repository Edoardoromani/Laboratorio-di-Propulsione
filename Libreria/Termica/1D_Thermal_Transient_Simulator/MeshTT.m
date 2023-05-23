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
  This code is the mesh object with which perform a thermal analysis. It
  stores the geometric, material, temperature and other variables to be
  used later in plots or other analyses

Changelog:
  > version: 1.0 - 13/11/2022 - Alessandro Rampazzo
    - created
%}

classdef MeshTT
    properties (SetAccess = private)
        % vector of node centroids
        x
        % matrix of temperature on centroids for every time step
        T
        % vector of material for every cells
        mat
        % vector of ablation front for every time step
        x_abl
        % time vector
        t
        % structure representing boundary conditions (must have a sx and a
        % dx field)
        BC
    end
    
    methods
        function obj = MeshTT(meshParams,T0)
            %% constructor
            % - meshParams: cell defining various sheets of materials as
            %     following 
            %     meshParams = {
            %       % thickness [m]     ID     N nodes
            %            5e-3,           2,       50;
            %            1e-3,           5,       20;
            %            7e-3,           1,       50};
            % - T0: initialization temperature, given as a value or a
            %     vector of the same size of x

            % parse meshParams and build mesh
            for i = 1:size(meshParams,1)

                thickness = meshParams{i,1};
                materialID = meshParams{i,2};
                N = meshParams{i,3};

                xx = linspace(0,thickness,N);
                if i == 1
                    obj.x = xx;
                    obj.mat = materialID * ones(1,N);
                else
                    obj.x = [obj.x, obj.x(end) + xx(2:end)];
                    obj.mat = [obj.mat, materialID * ones(1,N-1)];
                end
            end
            
            % initialize temperature, x_abl and time
            obj.x_abl = 0;
            obj.t = 0;
            if length(T0) == 1
                obj.T = ones(size(obj.x)) * T0;
            else
                obj.T = T0(:)';
            end
        end

        function obj = setBoundaryConditions(obj,BC)
            %% set boundary conditions
            obj.BC = BC;
            
            % find regression index
            I = find(obj.x < obj.x_abl);
            if isempty(I)
                regIdx = 0;
            else
                regIdx = I(end);
            end
            
            % set mesh temperature values to mach boudnary conditions
            if BC.sx.type == "temp" || BC.sx.type == "temperature"
                obj.T(1:max(regIdx,1)) = BC.sx.value;
            end

            if BC.dx.type == "temp" || BC.sx.type == "temperature"
                obj.T(end) = BC.dx.value;
            end
        end

        function CFL = calcCFLCondition(obj)
            %% calc Courant-Friedrichs-Lewy condition

            % find maximum heat transfer coefficient
            alpha_max = 0;
            for i = 1:length(obj.x)
                lambda = getLambda(obj.mat(i),obj.T(i));
                Cp = getCp(obj.mat(i),obj.T(i));
                rho = getRho(obj.mat(i),obj.T(i));
                alpha_max = max(alpha_max, lambda/(Cp*rho));
            end
            % dt_min = 1/2 * min(deltaX^2/alpha)
            CFL = 0.5 * (min(obj.x(1:end-1)-obj.x(2:end)))^2 / alpha_max;
        end

        function obj = propagate(obj,tspan,CFLFactor,enableAbl)
            %% propagate the simulation for a fixed amount of time
            % - tspan: time span, given as a number of second to propagate
            % - CFLFactor: factor to be multiplied to the CFL condition.
            %     low CFLFactors result in low dt
            % - enableAbl: enable ablation if the first left material 
            %     can ablate
            
            % default value if not given
            if ~exist("enableAbl","var")
                enableAbl = false;
            end
            % calc CFL conditions
            CFL = obj.calcCFLCondition();
            
            % print info
            fprintf("CFL condition is: dt < % 1.3e s\n",CFL)
            fprintf("using dt = % 1.3e s\n",CFL*CFLFactor)
            fprintf("\nstarting integration...\n\n")
            
            % create time vector
            tVec = 0:CFLFactor*CFL:tspan;
            
            % define print time function 
            tSkip = 2000;
            opts.Events = @(t,Y,dt) obj.printTime(t,Y,[tVec(1:tSkip:end-1),tVec(end)],dt);
            opts.RunEvents = "before";

            % initial conditions [T_0, T_1, T_2, ..., x_abl]
            y0 = [obj.T(end,:),obj.x_abl(end)];
            
            % solver, explicit intergrator has been used for stability
            S = solveOde(@(t,Y) FVM2(t,Y,obj.x,obj.mat,obj.BC,enableAbl),"EXP",tVec,y0,opts);
            
            % force values to correct some numerical errors due to ablation
            if enableAbl
                for i = 1:length(S.t)
                    xabl = S.y(i,end);
                    S.y(i,obj.x < xabl) = obj.BC.sx.value;
                end
            else
                for i = 1:length(S.t)
                    xabl = S.y(i,end);
                    ii = find(obj.x > xabl);
                    y = interp1(obj.x(ii:ii+1),S.y(i,ii:ii+1),xabl,'linear','extrap');
                    S.y(i,obj.x < xabl) = y;
                end
            end
            
            % add to the existing data the new ones
            obj.t = [obj.t;obj.t(end) + S.t(2:end)];
            obj.T = [obj.T;S.y(2:end,1:end-1)];
            obj.x_abl = [obj.x_abl;S.y(2:end,end)];
        end

        function [position,isterminal,direction] = printTime(~,t,Y,tVec,dt)
            %% event function to print time during integration

            % measure time
            if t == 0
                tic
            end
            tStep = toc;
            tic
            tspan = tVec(end);

            % print time only when t is equal to a value of tVec (a t
            % subset)
            I = find(tVec == t);

            % calculate seconds, minutes, hours and print them
            if mod(I,1) == 0
                ETA = (tspan - t)*tStep/dt;
                ETA_h = floor(ETA/3600);
                ETA_m = floor((ETA - ETA_h*3600)/60);
                ETA_s = floor(ETA - ETA_h*3600 - ETA_m*60);
                fprintf("t: %1.2f s / %1.1f s | ETA: %1.0f h %1.0f m %1.0f s\n",t,tspan,ETA_h,ETA_m,ETA_s)
            end

            % to detect divergence
            if any(isnan(Y)) || any(Y<0)
                warning("divergence detected: try reducing time step to meet CFL conditon")
                position = 0; % The value that we want to be zero
            else
                position = 1; % The value that we want to be zero
            end
            isterminal = 1;  % Halt integration
            direction = 0;   % The zero can be approached from either direction
        end
    end
end



