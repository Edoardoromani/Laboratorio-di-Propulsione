function geometry = BellNozzle(D_t,A_ratio,l_percent,n,R_link,entrant_angle)
% function to create a Bell nozzle
%
% D_t: throat diameter
% A_ratio: area ratio (A_e/A_t)
% l_percent: percentile length with respect to a cone nozzle with the same
% exit area and semi divergence angle = 15°
% n: mesh density (3 section with the same density, default = 100)
% R_link: linking section between throat and bell nozzle, normalized with
% respect to throat radius (default 0.5)
% entrant_angle: initial convergent angle 
%
% the function can be called with only two arguments: geometry = BellNozzle(0.01,6)

if nargin < 3
    l_percent = 80;
end
if nargin < 4
    n = 100;
end
if nargin < 5
    R_link = 0.5;
end
if nargin < 6
    entrant_angle = 45;
end

A_t = pi/4 * D_t^2;

Rt = sqrt(A_t / pi);
% find wall angles (theta_n, theta_e) for given aratio (ar)
[~,theta_n,theta_e] = find_wall_angles(A_ratio, Rt, l_percent);
%     angles = [theta_n,theta_e];
% entrant functions
ea_radian = deg2rad(entrant_angle);
ea_start = -pi/2 - ea_radian;
ea_end = -pi/2;
angle_list = linspace(ea_start, ea_end, n);
xe = zeros(size(angle_list));
ye = xe;
for i = 1:length(angle_list)
    xe(i) = 1.5 * Rt * cos(angle_list(i));
    ye(i) = 1.5 * Rt * sin(angle_list(i)) + 2.5 * Rt;
end

% exit section
ea_start = -pi/2;
ea_end = theta_n - pi/2;
angle_list = linspace(ea_start, ea_end, n);
xe2 = zeros(size(angle_list));
ye2 = xe2;
R_link = Rt * R_link;
for i = 1:length(angle_list)
    xe2(i) = R_link * cos(angle_list(i));
    ye2(i) = R_link * sin(angle_list(i)) + R_link + Rt;
end

% bell section
% Nx, Ny-N is defined by [Eqn. 5] setting the angle to (θn – 90)
Nx = R_link * cos(theta_n - pi / 2);
Ny = R_link * sin(theta_n - pi / 2) + R_link + Rt;
% Ex - [Eqn. 3], and coordinate Ey - [Eqn. 2]
Ex = l_percent/100 * ((sqrt(A_ratio) - 1) * Rt) / tan(deg2rad(15));
Ey = sqrt(A_ratio) * Rt;
% gradient m1,m2 - [Eqn. 8]
m1 = tan(theta_n);
m2 = tan(theta_e);
% intercept - [Eqn. 9]
C1 = Ny - m1 * Nx;
C2 = Ey - m2 * Ex;
% intersection of these two lines (at point Q)-[Eqn.10]
Qx = (C2 - C1) / (m1 - m2);
Qy = (m1 * C2 - m2 * C1) / (m1 - m2);

% Selecting equally spaced divisions between 0 and 1 produces
% the points described earlier in the graphical method
% The bell is a quadratic Bézier curve, which has equations:
% x(t) = (1 − t)^2 * Nx + 2(1 − t)t * Qx + t^2 * Ex, 0≤t≤1
% y(t) = (1 − t)^2 * Ny + 2(1 − t)t * Qy + t^2 * Ey, 0≤t≤1 [Eqn. 6]
int_list = linspace(0, 1, n);
xbell = zeros(size(int_list));
ybell = xbell;
for i = 1:length(int_list)
    t = int_list(i);
    xbell(i) = ((1 - t) ^ 2) * Nx + 2 * (1 - t) * t * Qx + (t ^ 2) * Ex;
    ybell(i) = ((1 - t) ^ 2) * Ny + 2 * (1 - t) * t * Qy + (t ^ 2) * Ey;
end
% create negative values for the other half of nozzle
%     data = {xe(:),ye(:),xe2(:),ye2(:),xbell(:),ybell(:)};
geometry.x = [xe';xe2(2:end)';xbell(2:end)'] - min(xe);
geometry.r = [ye';ye2(2:end)';ybell(2:end)'];
geometry.A = pi*geometry.r.^2;


% find wall angles (theta_n, theta_e) in radians for given aratio (ar)
function [Ln, theta_n, theta_e] = find_wall_angles(ar, Rt, l_percent)
    % wall-angle empirical data
    if nargin < 3
        l_percent = 80;
    end
    Aratio = [4, 5, 10, 20, 30, 40, 50, 100];
    theta_n_60 = [20.5, 20.5, 16.0, 14.5, 14.0, 13.5, 13.0, 11.2];
    theta_n_80 = [21.5, 23.0, 26.3, 28.8, 30.0, 31.0, 31.5, 33.5];
    theta_n_90 = [20.0, 21.0, 24.0, 27.0, 28.5, 29.5, 30.2, 32.0];
    theta_e_60 = [26.5, 28.0, 32.0, 35.0, 36.2, 37.1, 35.0, 40.0];
    theta_e_80 = [14.0, 13.0, 11.0, 9.0, 8.5, 8.0, 7.5, 7.0];
    theta_e_90 = [11.5, 10.5, 8.0, 7.0, 6.5, 6.0, 6.0, 6.0];

    % nozzle length
    f1 = ((sqrt(ar) - 1) * Rt) / tan(deg2rad(15));

    if l_percent == 60
        theta_n = theta_n_60;
        theta_e = theta_e_60;
        Ln = 0.8 * f1;
    elseif l_percent == 80
        theta_n = theta_n_80;
        theta_e = theta_e_80;
        Ln = 0.8 * f1;
    elseif l_percent == 90
        theta_n = theta_n_90;
        theta_e = theta_e_90;
        Ln = 0.9 * f1;
    else
        theta_n = theta_n_80;
        theta_e = theta_e_80;
        Ln = 0.8 * f1;
    end
    % find the nearest ar index in the aratio list

    [~,x_index] = min(abs(Aratio - ar));
    % if the value at the index is close to input, return it
    if round(Aratio(x_index), 1) == round(ar, 1)
        theta_n = deg2rad(theta_n(x_index));
        theta_e = deg2rad(theta_e(x_index));
        return
    end
    % check where the index lies, and slice accordingly
    if x_index > 3
        % slice couple of middle values for interpolation
        ar_slice = Aratio(x_index - 1:x_index + 2);
        tn_slice = theta_n(x_index - 1:x_index + 2);
        te_slice = theta_e(x_index - 1:x_index + 2);
        % find the tn_val for given ar
        tn_val = interp1(ar_slice, tn_slice, ar);
        te_val = interp1(ar_slice, te_slice, ar);
    elseif (length(Aratio) - x_index) <= 2
        % slice couple of values initial for interpolation
        ar_slice = Aratio(x_index - 1:len(x_index));
        tn_slice = theta_n(x_index - 1:len(x_index));
        te_slice = theta_e(x_index - 1:len(x_index));
        % find the tn_val for given ar
        tn_val = interp1(ar_slice, tn_slice, ar);
        te_val = interp1(ar_slice, te_slice, ar);
    else
        % slice couple of end values for interpolation
        ar_slice = Aratio(1:x_index + 2);
        tn_slice = theta_n(1:x_index + 2);
        te_slice = theta_e(1:x_index + 2);
        % find the tn_val for given ar
        tn_val = interp1(ar_slice, tn_slice, ar);
        te_val = interp1(ar_slice, te_slice, ar);
    end
    theta_n = deg2rad(tn_val);
    theta_e = deg2rad(te_val);
end
end
