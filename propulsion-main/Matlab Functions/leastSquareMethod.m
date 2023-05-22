function [c,f_pred_out,R2,err] = leastSquareMethod(X,p,f_pred_in,vars,c0,varargin)
%LEASTSQUAREMETHOD finds the best fit curve for a given set of data using a
%predictor function
%
% Property of THRUST, unauthorized distribution is not allowed
% version: 1.3
%
%   C = LEASTSQUAREMETHOD(X,P,F_PRED,VARS,C0) finds the coefficient C 
%   that best approximate the given set of data P with input X using a
%   predictor function. P must be a vector and X something to be
%   substituted into the predictor function or a cell of them (if they are
%   more than one). F_PRED is a symbolic function and VARS is the list of
%   all symbolic variables that appear in the function in order and they
%   are the coefficient which will be optimized. F_PRED must include 
%   length(X) symbolic variable named x1, x2 x3... (or x if only one is 
%   present) which must not be included in the VARS vector. C0 is the
%   initial guesses for the coefficients
%
%   [C,F_PRED_OUT] = LEASTSQUAREMETHOD(X,P,F_PRED,VARS,C0) returns the 
%   matlab function F_PRED_OUT which is the optimized function for the
%   given data
%
%   [C,F_PRED_OUT,R2] = LEASTSQUAREMETHOD(X,P,F_PRED,VARS,C0) returns
%   coefficient of determination of the interpolaton. The closer to 1 it is
%   the better is the interpolation
%
%   [C,F_PRED_OUT,R2,ERR] = LEASTSQUAREMETHOD(X,P,F_PRED,VARS,C0) returns
%   the error of the final iteration
%   
%   [C,F_PRED_OUT,R2,ERR] = LEASTSQUAREMETHOD(X,P,F_PRED,VARS,C0,param,
%   value) modify internal parameters for the calculation. Available 
%   parameters are: 
%     - tol: convergence tolerance. Default = 1e-6
%     - itMax: maximum number of iterations. Default = 100
%     - URF: Under Relaxation Factor. Default = 0.5
%     - uncertainty: vector of variances of every data point to allow the
%     wheighted fit. Default = ones(size(X))
%     - verbosity = type of output. Available options: ["none","end","full"]
%     Default = "end"
% 
%   Example: simple linear regression
%
%     Data = [1,1;..
%           2,3.5;..
%            ...];
%     syms a b x
%     f_pred = a*x + b;
%     [c,f_pred] = LeastSquareMethod(Data(:,1),Data(:,2),f_pred,[a,b],[1,1]);
%
%   See also POLYFIT, POLYVAL

%{
Changelog:
  > version: 1.5 - 06/01/2023 - Alessandro Rampazzo
    - added the coefficient of determination R2

  > version: 1.4 - 03/01/2023 - Alessandro Rampazzo
    - normalization of hessian matrix

  > version: 1.3 - 13/11/2022 - Alessandro Rampazzo
    - implemented the multidimensional optimization by allowing the
    predictor function to have multiple input variables (x1, x2, etc)

  > version: 1.2 - 01/09/2022 - Alessandro Rampazzo
    - added the automatic jacobian calculation and optimized function is
    now given as an output

  > version: 1.1 - 02/06/2022 - Alessandro Rampazzo
    - fixed some bugs

  > version: 1.0 - 01/06/2022 - Alessandro Rampazzo
    - created
%}

if mod(length(varargin),2) ~= 0
    error("fields and specifications should always come in pair")
end

% parse options
tol = 1e-6;
itMax = 100;
URF = 0.5;
uncertainty = ones(size(p));
verbosity = "end";
verbosityOptions = ["none","end","full"];
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case "urf"
            URF = varargin{i+1};
        case "itmax"
            itMax = varargin{i+1};
        case "tol"
            tol = varargin{i+1};
        case "uncertainty"
            uncertainty = varargin{i+1}/max(varargin{i+1});
        case "verbosity"
            verbosity = lower(varargin{i+1});
            if all(verbosity ~= verbosityOptions)
                error("verbosity option ""%s"" not recognized",verbosity)
            end
        otherwise
            error("option not recognized")
    end
end

% assert that the predictor function and vars are consistent
svar = symvar(f_pred_in);

if ~iscell(X) || (iscell(X) && size(X,2) == 1)
    if ~any("x" == string(svar))
        error("predictor function must have a symbolic variable named 'x'")
    else
        subsvar = svar("x" == string(svar));
    end
    if iscell(X)
        X{1} = X{1}(:);
    else
        X = X(:);
    end
else
    cc = 1;
    for i=1:size(X,2)
        if ~any("x"+i == string(svar))
            error("predictor function must have a symbolic variable named 'x'")
        else
            subsvar(cc) = svar("x"+i == string(svar));
            cc = cc + 1;
        end
        if iscell(X)
            X{i} = X{i}(:);
        end
    end
end

for i = 1:length(vars)
    if ~any(string(vars(i)) == string(svar))
        error("number of input variables does not coincide with numeber of" + ...
            " variables present in predictor function")
    end
end

if length(vars) ~= length(c0)
    error("number of symbolic variables does not match with number of initial guesses")
end

% substitute the inputs into the predictor function
f_pred = subs(f_pred_in,subsvar,X);


% generate Jacobian matrix
J_func = [];
for i = 1:length(vars)
    J_func = [J_func,diff(f_pred,vars(i))];
end

% transform f_pred e J_func into matlab functions
J_func = matlabFunction(J_func,'Vars',vars);
f_pred = matlabFunction(f_pred,'Vars',vars);
str_f = "f_pred = @(c) f_pred(";
str_J = "J_func = @(c) J_func(";
for i = 1:length(vars)
    if i > 1
        str_f = str_f + ",";
        str_J = str_J + ",";
    end
    str_f = str_f + "c(" + i + ")";
    str_J = str_J + "c(" + i + ")";
end
str_f = str_f + ");";
str_J = str_J + ");";

eval(str_f)
eval(str_J)

% if det() == 0
%     error("jacobian matrix is not invertible, remove the excess redundant samples")
% end
% % assert that shapes are correct
% if size(J_func(c0),1) ~= length(p)
%     error("shape mismatch, transpose the f_pred function")
% end

% make everything columnwise
c0 = c0(:);
p = p(:);
uncertainty = uncertainty(:);

% intialize coefficients
c = c0;

% main loop
err = 1e+50;
iter = 0;
while err > tol && iter < itMax
    iter = iter + 1;
    % calculate Jacobian
    J = J_func(c);
    % normalize derivatives
    maxs = max(abs(J));
    J = J./maxs;
    % calculate residuals
    y = (p - f_pred(c))./(uncertainty);
    % calculate step
    delta_c = URF * ((J.'*J)\J.'*y);   
    % denormalize delta_c
    delta_c = delta_c ./ maxs(:);
    % calculate error
    err = norm(delta_c./c0);
    % output infos
    if verbosity == "full"
        fprintf("iter %i | error on data: %1.3e | delta x: %1.3e\n",iter,norm(y./p)/sqrt(length(y)),err)
    end   
    % add step
    c = c + delta_c;
end

if any(or(isnan(c),isinf(c)))
    error("nan or inf encoutered, try modifying the predictor function " + ...
        "and make sure that every parameter is indipendent")
end

% calculate the coefficient of determination
R2 = 1 - sum((p - f_pred(c)).^2)/sum((p - mean(p)).^2);

% post main loop checks
if iter == itMax
    warning("maximum iterations reached")
elseif verbosity == "full" || verbosity == "end"
    fprintf("\nbest fit curve found:\niter %i | R2 = %1.4f | error on data: %1.3e | delta x: %1.3e\n\n",iter,R2,norm(y./p)/sqrt(length(y)),err)
end

% tranform f_pred in a matlabFunction with the optimized coefficients
f_pred_out = f_pred_in;

for i = 1:length(vars)
    f_pred_out = subs(f_pred_out,vars(i),c(i));
end

f_pred_out = matlabFunction(f_pred_out,'Vars',symvar(f_pred_out));

end
