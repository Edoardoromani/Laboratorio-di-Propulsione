

%%      VENTING HOLES ACCORDING TO THEORY BY JOHN J. SCIALDONE, NASA, "SPACECRAFT COMPARTMENT VENTING"

%   initialization
clear all
close all
clc

%   data

dp_dt = 0.042e5;    %   [Pa/s]  External pressure change rate - from simulator

V = 0.0037;          %   [m^3]   Internal volume (Upper part) - 0.0037 Lower part

C = 0.65;           %   [ ]     Orifice coefficient

g = 9.81;           %   [m/s^2] Gravity acceleration

R = 29.2;           %   [m/K] (should be 287 J/KgK but he divided it by g = 9.81 m/s^2)

T = 216;            %   [K]     Air temperature

P0 = 0.2e5;         %   [Pa]    Initial pressure at which steep pressure variation occurs

DP = 0.01e5;        %   [Pa]    Max delta P accepted bewteen interior and esterior

%   calculations

A = dp_dt * V / (C*sqrt(2*g*R*T*DP*P0));    %   [m^2]
d = sqrt(A/pi)*2;                           %   [m]
d = d*1000;                                 %   [mm]

fprintf("%4f\n",d)