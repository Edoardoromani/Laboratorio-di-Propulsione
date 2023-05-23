clear all
close all
clc

Ncr = 2500; %[N]
E = 70000; % [Mpa]
L = 1150; % [mm]
Lo = L; % fattore 1: trave con doppia cerniera, conservativo

Imin = Ncr*Lo^2/pi^2/E; % Inerzia minima della sezione

H_quad =  nthroot(12*Imin,4); % Lato quadrato pieno

h = H_quad-2.4;
H_cavo = nthroot(12*Imin + h^4,4);

fprintf("\n\nLato quadrato pieno: L = %2.4f mm",H_quad)
fprintf("\n\nLato quadrato cavo : L' = %2.4f mm..." + ...
          "\nSpessore           : t = %2.4f mm",H_cavo,(H_cavo-h)/2)

A_quad = H_quad^2;
A_cavo = H_cavo^2-h^2;

fprintf("\n\n Area Quadrato Pieno: A = %2.4f mm^2",A_quad)
fprintf("\n\n Area Quadrato Cavo : A' = %2.4f mm^2",A_cavo)