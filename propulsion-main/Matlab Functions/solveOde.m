function varargout = solveOde(odeFun,scheme,tSpan,y0,varargin)
%SOLVEODE numerically solves an ode or a set of ode 
%
% Property of THRUST, unauthorized distribution is not allowed
% version: 1.4
% DO NOT USE THE "FIND EVENT" METHOD WITH IMPLICIT METHODS (IMP and CN)
%
%   See also ODE23, ODE45

%{
Changelog:
  > version: 1.5 - 27/01/2023 - Alessandro Rampazzo
    - changed minimum allowable dt to 16*eps (as in ode45)

  > version: 1.4 - 16/12/2022 - Alessandro Rampazzo
    - added header
    - fixed bugs regarding events
    - removed the option runEvent = "before"
    - changed the way the options are given in input, now they are no
    longer given with a struct but by fields and specifications

  > version: 1.3 - 05/06/2022 - Alessandro Rampazzo
    - fixed some bugs
    - added events to stop integration

  > version: 1.2 - 01/04/2022 - Alessandro Rampazzo
    - implemented explicit method, implicit method, crank nicholson and
    runge kutta 45 adaptive

  > version: 1.1 - 02/04/2022 - Alessandro Rampazzo
    - fixed some bugs

  > version: 1.0 - 01/06/2022 - Alessandro Rampazzo
    - created (only runge kutta 45)
%}

%% options
if mod(length(varargin),2) ~= 0
    error("fields and specifications should always come in pair")
end

relTol = 1e-5;
absTol = 1e-5;
initialDt = 1e-5;
itMax = 100;
maxDTFactor = 0.05;
breakEvent = @(t,Y) 1;
% runEvent is deprecated
runEvent = "after";
includeEvent = false;
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case "reltol"
            relTol = varargin{i+1};
        case "abstol"
            absTol = varargin{i+1};
        case "initialdt"
            initialDt = varargin{i+1};
        case "itmax"
            itMax = varargin{i+1};
        case "maxdtfactor"
            maxDTFactor = varargin{i+1};
        case "breakevent"
            breakEvent = varargin{i+1};
        case "runevent"
            % runEvent = varargin{i+1};
        case "includeevent"
            includeEvent = varargin{i+1};
        otherwise
            error("option not recognized")
    end
end
% old way of defining options through struct
% if exist('options','var')
%     if isfield(options,'RelTol')
%         relTol = options.RelTol;
%     end
%     if isfield(options,'AbsTol')
%         absTol = options.AbsTol;
%     end
%     if isfield(options,'MaxIter')
%         maxIter = options.MaxIter;
%     end
%     if isfield(options,'Event')
%         breakEvent = options.Event;
%     end
%     if isfield(options,'RunEvent')
% %         runEvent = lower(options.RunEvent);
%         if runEvent ~= "after" && runEvent ~= "before"
%             error("RunEvents option don't recognized")
%         end
%     end
%     if isfield(options,"IncludeEvent")
%         includeEvent = options.IncludeEvent;
%     end
%     if isfield(options,"InitialDt")
%         initialDt = options.InitialDt;
%     end
%     if isfield(options,"maxDTFactor")
%         maxDTFactor = options.maxDTFactor;
%     end
% end

%% check consistency of intergration variable
tSpan = tSpan(:)';
if any(sign(diff(tSpan)) ~= 1) && any(sign(diff(tSpan)) ~= -1)
    error("integration variable must monotonically increase or decrease!")
end

%% initialization
if scheme == "RK45_A"
    y = zeros(1,length(y0));
    y(1,:) = y0;
    t = tSpan(1);
    t = t(:);
else
    y = zeros(length(tSpan(:)),length(y0));
    y(1,:) = y0;
    t = tSpan(:);
end

maxDT = maxDTFactor*abs(tSpan(1)-tSpan(end));

initialDt = initialDt * (1 - 2*(tSpan(end) < tSpan(1)));
sdt = sign(initialDt);

includeEventBreakFlag = false;

%% solver
switch lower(scheme)
    case "exp"
        %% Explicit
        for i = 2:length(t)
            dt = t(i) - t(i-1);
%             [doBreak,includeEventBreakFlag] = runEventsBefore(runEvent,includeEvent,breakEvent,t(i-1),y(i-1,:));
%             if doBreak
%                 break
%             end

            y(i,:) = y(i-1,:) + odeFun(t(i-1),y(i-1,:))*dt;
            doBreak = runEventsAfter(runEvent,includeEventBreakFlag,breakEvent,t(i),y(i,:));
            if doBreak
                break
            end
        end
    case "imp"
        %% Implicit
        for i = 2:length(t)
            dt = t(i) - t(i-1);

%             [doBreak,includeEventBreakFlag] = runEventsBefore(runEvent,includeEvent,breakEvent,t(i-1),y(i-1,:));
%             if doBreak
%                 break
%             end
            
            relErr = 1e+50;
            absErr = 1e+50;
            iter = 0;
            y(i,:) = y(i-1,:);
            while (relErr > relTol || absErr > absTol) && iter < itMax
                iter = iter + 1;
                yOld = y(i,:);
                y(i,:) = y(i-1,:) + odeFun(t(i),y(i,:))*dt;
                absErr = max(abs(y(i,:) - yOld));
                yOldCorrected = yOld;
                yOldCorrected(yOld == 0) = 1e+50;
                relErr = max(abs((y(i,:) - yOld)./yOldCorrected));
            end
            doBreak = runEventsAfter(runEvent,includeEventBreakFlag,breakEvent,t(i),y(i,:));
            if doBreak 
                break
            end
        end
    case "cn"
        %% Crank-Nicolson
        for i = 2:length(t)
            dt = t(i) - t(i-1);

%             [doBreak,includeEventBreakFlag] = runEventsBefore(runEvent,includeEvent,breakEvent,t(i-1),y(i-1,:));
%             if doBreak
%                 break
%             end
            
            relErr = 1e+50;
            absErr = 1e+50;
            iter = 0;
            y(i,:) = y(i-1,:);
            while (relErr > relTol || absErr > absTol) && iter < itMax
                iter = iter + 1;
                yOld = y(i,:);
                y(i,:) = y(i-1,:) + 0.5*dt*(odeFun(t(i),y(i,:)) + odeFun(t(i-1),y(i-1,:)));
                absErr = max(abs(y(i,:) - yOld));
                yOldCorrected = yOld;
                yOldCorrected(yOld == 0) = 1e+50;
                relErr = max(abs((y(i,:) - yOld)./yOldCorrected));
            end

            doBreak = runEventsAfter(runEvent,includeEventBreakFlag,breakEvent,t(i),y(i,:));
            if doBreak
                break
            end
        end
    case "rk45"
        %% Runge-Kutta 45
        a = [0,0,0,0,0,0;
            1/4,0,0,0,0,0;
            3/32,9/32,0,0,0,0;
            1932/2197,-7200/2197,7296/2197,0,0,0;
            439/216,-8,3680/513,-845/4104,0,0;
            -8/27,2,-3544/2565,1859/4104,-11/40,0];

        b = [16/135,0,6656/12825,28561/56430,-9/50,2/55];
        c = [0,1/4,3/8,12/13,1,1/2];
        k = zeros(size(a,1),length(y0));

        for i = 2:length(t)
            dt = t(i) - t(i-1);

%             [doBreak,includeEventBreakFlag] = runEventsBefore(runEvent,includeEvent,breakEvent,t(i-1),y(i-1,:));
%             if doBreak
%                 break
%             end

            for j = 1:size(k,1)
                k(j,:) = odeFun(t(i-1) + c(j)*dt,y(i-1,:) + dt*a(j,:)*k);
            end

            y(i,:) = y(i-1,:) + dt*b*k;
            doBreak = runEventsAfter(runEvent,includeEventBreakFlag,breakEvent,t(i),y(i,:));
            if doBreak
                break
            end
        end
    case "rk45_a"
        %% Runge-Kutta 45 adaptive
%       Runge–Kutta–Fehlberg
        a = [0,0,0,0,0,0;
            1/4,0,0,0,0,0;
            3/32,9/32,0,0,0,0;
            1932/2197,-7200/2197,7296/2197,0,0,0;
            439/216,-8,3680/513,-845/4104,0,0;
            -8/27,2,-3544/2565,1859/4104,-11/40,0];

        b = [16/135,0,6656/12825,28561/56430,-9/50,2/55];
        c = [0,1/4,3/8,12/13,1,1/2];
        bt = [1/360,0,-128/4275,-2197/75240,1/50,2/55];

%        Dormand-Prince
%         a = [0             0            0            0          0                0         0
%              1/5           0            0            0          0                0         0
%              3/40          9/40         0            0          0                0         0
%              44/45        -56/15        32/9         0          0                0         0
%              19372/6561   -25360/2187   64448/6561  -212/729    0                0         0
%              9017/3168    -355/33       46732/5247   49/176     -5103/18656      0         0
%              35/384        0            500/1113     125/192    -2187/6784       11/84     0  
%             ];
%         b =  [35/384        0            500/1113     125/192    -2187/6784       11/84    0];
%         c = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
%         bt = [71/57600, 0, -71/16695, 71/1920, -17253/339200, 22/525, -1/40];
        
        k = zeros(size(a,1),length(y0));

        dt = initialDt;
        alreadyWarned = false;
        relDt = 0;
        absDt = 0;
        while t(end)*sdt < tSpan(end)*sdt

%             [doBreak,includeEventBreakFlag] = runEventsBefore(runEvent,includeEvent,breakEvent,t(end),y(end,:));
%             if doBreak
%                 break
%             end

            for j = 1:size(k,1)
                k(j,:) = odeFun(t(end) + c(j)*dt,y(end,:) + dt*a(j,:)*k);
            end

%             yy = y(end,:) + dt*b*k;

            % absolute error
            absErr = max(abs(bt * k));

            % relative error
            yend = abs(y(end,:));
            yend(yend == 0) = 1e+50;
            relErr = max(abs(bt * k./yend));

            if relErr < relTol && absErr < absTol
                t(end+1) = t(end) + dt;
                y(end+1,:) = y(end,:) + dt*b*k;
                doBreak = runEventsAfter(runEvent,includeEventBreakFlag,breakEvent,t(end),y(end,:));
%                 fprintf("relDt = %1.3e | absDt = %1.3e | dt = %1.3e | relErr = %1.3e | absErr = %1.3e\n",relDt,absDt,dt,relErr,absErr)
%                 disp(yend)
%                 disp(k)
                if doBreak
                    break
                end
            end
           
            relDt = 0.8 * abs(dt) * (relTol/relErr)^0.2;
            absDt = 0.8 * abs(dt) * (absTol/absErr)^0.2;
            dt = min([relDt,absDt,maxDT,abs(t(end) - tSpan(end))]);
            dt = sdt * max(dt,eps*16);
           
% old way of calculating dt         
%             dt = sdt * min([max(0.9 * abs(dt) * (tol/err)^0.2, 1e-12), ...
%             maxDT, abs(t(end) - tSpan(end))]);
            if abs(t(end) - tSpan(end)) == 0
                break
            end
            if dt == eps*16 && ~alreadyWarned
                alreadyWarned = true;
                warning("smallest possible time step reached, results may be inaccurate")
            end                
        end
    otherwise
        error('scheme not recognized')
end

% fixing output variables
if scheme == "RK45_A"
    if doBreak && ~includeEvent
        t = t(1:end-1);
        y = y(1:end-1,:);
    end

%     if doBreak && runEvent == "before"
%         t = t(1:end-1);
%         y = y(1:end-1,:);
%     end

    if length(tSpan) > 2
        S.t = tSpan(tSpan*sdt <= t(end)*sdt);
        S.y = zeros(length(S.t),size(y,2));
        for i = 1:size(S.y,2)
            S.y(:,i) = interp1(t,y(:,i),S.t)';            
        end
    else
        S.y = y;
        S.t = t(:);
    end
else
    if ~includeEvent
        S.t = t(1:i-1);
        S.y = y(1:i-1,:);
    else
        S.t = t(1:i);
        S.y = y(1:i,:);
    end
end

if nargout == 1
    varargout = {S};
elseif nargout == 2
    varargout = {S.t,S.y};
else
    error("too many output variables")
end

%% Nested functions
    function [doBreak,includeEventBreakFlag] = runEventsBefore(runEvent,includeEvent,breakEvent,t,y)
    doBreak = false;
    includeEventBreakFlag = false;
    if runEvent == "before"
        cont = breakEvent(t,y);
        if cont <= 0
            if includeEvent
                includeEventBreakFlag = true;
            else
                doBreak = true;
            end
        end
    end
end

function doBreak = runEventsAfter(runEvent,includeEventBreakFlag,breakEvent,t,y)
    doBreak = false;
    if includeEventBreakFlag
        doBreak = true;
        return
    end

    if runEvent == "after"
        cont = breakEvent(t,y);
        if cont <= 0
            doBreak = true;
        end
    end
end

end