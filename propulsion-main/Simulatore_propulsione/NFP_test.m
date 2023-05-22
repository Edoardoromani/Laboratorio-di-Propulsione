clc
clear
close all

addpath("Funzioni sottosistemi")

rhol_func = @(T) densityl_N2O(T);
rhol_NFP_func = @(T) NFP("N2O","rl_t",T);
N = 10000;

tic
for i = 1:N
    rhol_func(293.15);
end
t1 = toc;
fprintf("normal interpolation took: %1.3e s\n",t1/N)

tic
for i = 1:N
    rhol_NFP_func(293.15);
end
t2 = toc;
fprintf("NFP interpolation took: %1.3e s\n",t2/N)

%%
T0 = 293.15;
m0_ox = 4;
V_tank = 7.25e-3;
k = 1.33;

%%
tic
[m_normal,p_normal,rho_normal,T_normal,x_normal] = tankFunction(T0,m0_ox,V_tank,k);
t4 = toc;
fprintf("Normal tank simulation took: %1.3e s\n",t4)
%%
tic
[m_NFP,p_NFP,rho_NFP,T_NFP,x_NFP] = tankFunctionNFP(T0,m0_ox,V_tank,k);
s_NFP = m_NFP*0;

for i = 1:length(m_NFP)
    if x_NFP(i) < 1
        s_NFP(i) = NFP("N2O","s_px",p_NFP(i)/1e5,x_NFP(i));
    else
        s_NFP(i) = NFP("N2O","s_pt",p_NFP(i)/1e5,T_NFP(i));      
    end
end

t3 = toc;
fprintf("NFP tank simulation took: %1.3e s\n",t3)
%%
load("S")

%%
figure
hold on
plot(m_NFP,T_NFP,'o-',"LineWidth",1)
plot(m_normal,T_normal,'-.',"LineWidth",1)
plot(S.m,S.T,'--',"LineWidth",1)
legend("NFP","normal","NFP old","location","best")
title("Temperature")

figure
hold on
plot(m_NFP,p_NFP,'o-',"LineWidth",1)
plot(m_normal,p_normal,'-.',"LineWidth",1)
plot(S.m,S.p,'--',"LineWidth",1)
legend("NFP","normal","NFP old","location","best")
title("Pressure")

figure
hold on
% mm = linspace(m_NFP(1),m_NFP(end),5000);
% plot(mm,interp1(m_NFP,rho_NFP,mm,'pchip'),'o-',"LineWidth",1)
plot(m_NFP,rho_NFP,'o-',"LineWidth",1)
plot(m_normal,rho_normal,'-.',"LineWidth",1)
plot(S.m,S.rho_out,'--',"LineWidth",1)
legend("NFP","normal","NFP old","location","best")
title("Density")

figure
hold on
plot(m_NFP,x_NFP,'o-',"LineWidth",1)
plot(m_normal,x_normal,'-.',"LineWidth",1)
plot(S.m,S.x,'--',"LineWidth",1)
legend("NFP","normal","NFP old","location","best")
title("Vapor quality")



figure
hold on
plot_bell_curve('r')
plot_isobars([2,5,10,25,50],'g')
grid on
box on
xlabel("s [kJ/kgK]")
ylabel("T [K]")
plot(s_NFP,T_NFP,'b')
plot(s_NFP(1),T_NFP(1),'bo')
text(s_NFP(1),T_NFP(1)+5,'start')
plot(s_NFP(end),T_NFP(end),'bo')
text(s_NFP(end),T_NFP(end)+5,'end')

function plot_bell_curve(plot_arg)
    N = 100;
    s_crit = NFP("N2O","s_px",NFP("N2O","pcrit"),0);
    s_val = [linspace(-0.0064,s_crit-1e-6,N),linspace(s_crit+1e-6,2.032,N)];
    T_val = s_val*0;
    for i = 1:length(s_val)
        if s_val(i) < s_crit
            T_val(i) = fzero(@(T) NFP("N2O","sl_t",T) - s_val(i),[183,NFP("N2O","tcrit")]);
        else
            T_val(i) = fzero(@(T) NFP("N2O","sv_t",T) - s_val(i),[183,NFP("N2O","tcrit")]);
        end
    end
    plot(s_val,T_val,plot_arg)
end

function plot_isobars(p_val,plot_arg)
    N = 100;
    s_val = zeros(1,N);
    T_val = linspace(185,308,N);

    for j = 1:length(p_val)
        offset = 10;
        for i = 1:length(T_val)
            s_val(i) = NFP("N2O","s_pt",p_val(j),T_val(i));
        end
        plot(s_val,T_val,plot_arg)
        ss = s_val(end-offset);
        TT = T_val(end-offset);
        dar = get(gca,"DataAspectRatio");
        ang = atan2d(dar(1)* (T_val(end-offset+1) - T_val(end-offset-1)),...
            dar(2) * (s_val(end-offset+1) - s_val(end-offset-1)));
        text(ss,TT,sprintf("%1.1f bar",p_val(j)),'rotation',ang,...
            "FontSize",10)
    end
end


