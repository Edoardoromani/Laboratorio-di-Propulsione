function [cl,cd]=aero_coeff_fun(Ma,r,alfa)


% caso motore on (r=on);
cd0=0.45;
if Ma<0.8
    cd=cd0;
elseif Ma>4
    cd=0.85*(2-4^0.4);
elseif Ma>1.1
    cd=0.85*(2-Ma^0.4);
elseif Ma<1
    cd=cd0+(Ma-0.8)*(0.75-cd0)/(0.2);

else
    cd=0.75+(Ma-1)*(0.816968899-0.75)/(0.1);
end
    
% caso motore off (r=off);
if strcmp(r,'off')==1   % se le due stringe sono uguali
    cd=cd+0.1;
end

a0=0.1665;

if Ma<0.6
    cl=a0*alfa;
elseif Ma<0.9
    %cl=a0*alfa/sqrt(1-Ma(i)^2);
    cl=a0*alfa*(1+1.3*(Ma-0.6));
elseif Ma>1.1
    cl=0.26*alfa*(1-0.1*Ma);
else 
    cl=0.26*alfa*(1-0.1*1.1);
end

