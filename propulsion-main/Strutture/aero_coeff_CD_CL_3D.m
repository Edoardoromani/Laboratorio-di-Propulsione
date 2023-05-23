%   function [cl,cd]=aero_coeff3(Ma,r,alfa)
clear all
close all
clc

Ma = 0:0.025:2;
alfa = -4:0.1:4;
CDs = ones(length(Ma),length(Ma));
CLs = ones(length(Ma),length(Ma));

for i = 1:length(Ma)
    for j = 1:length(alfa)

        % caso motore on (r=on);
        cd0=0.45;
        if Ma(i)<0.8
            cd=cd0;
        elseif Ma(i)>4
            cd=0.85*(2-4^0.4);
        elseif Ma(i)>1.1
            cd=0.85*(2-Ma(i)^0.4);
        elseif Ma(i)<1
            cd=cd0+(Ma(i)-0.8)*(0.75-cd0)/(0.2);
        else
            cd=0.75+(Ma(i)-1)*(0.816968899-0.75)/(0.1);
        end
            
        % caso motore off (r=off);
        %if strcmp(r,'off')==1   % se le due stringe sono uguali
        cd=cd+0.1;
        %end
        
        a0=0.1665;
        
        if Ma(i)<0.6
            cl=a0*alfa(j);
        elseif Ma(i)<0.9
            %cl=a0*alfa/sqrt(1-Ma(i)^2);
            cl=a0*alfa(j)*(1+1.3*(Ma(i)-0.6));
        elseif Ma(i)>1.1
            cl=0.26*alfa(j)*(1-0.1*Ma(i));
        else 
            cl=0.26*alfa(j)*(1-0.1*1.1);
        end

        CLs(i,j) = cl;
        CDs(i,j) = cd;

    end
end

CDX = max(CDs,[],"all");
[CdR,CdC] = find(CDs==CDX);
Cd_Ma_max = Ma(CdR);
Cd_alfa_max = alfa(CdC);
Cd_z = ones(1,length(Cd_alfa_max))*CDX

CLX = max(CLs,[],"all")
[ClR,ClC] = find(CLs==CLX)
Cl_Ma_max = Ma(ClR)
Cl_alfa_max = alfa(ClC)
Cl_z = ones(1,length(Cl_alfa_max))*CLX

figure(1)
surf(alfa,Ma,CDs)
hold on
xlabel("alpha")
ylabel("Mach")
zlabel("Cd")
plot3(Cd_alfa_max,Cd_Ma_max,Cd_z,"or")
hold off

figure(2)
surf(alfa,Ma,CLs)
hold on
xlabel("alpha")
ylabel("Mach")
zlabel("Cl")
plot3(Cl_alfa_max,Cl_Ma_max,Cl_z,"or")
hold off