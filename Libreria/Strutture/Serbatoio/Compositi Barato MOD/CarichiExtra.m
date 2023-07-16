clear 
close all
clc

dist1 = 1.774;                      %   [m] distanza inizio serbatoio da punta del razzo
dist2 = 3.678;                      %   [m] distanza inizio serbatoio da punta del razzo
% dist1 = 0.35;                       %   [m] distanza inizio serbatoio da punta del razzo
% dist2 = 4.072;                      %   [m] distanza inizio serbatoio da punta del razzo

dist = linspace(dist1,dist2,100);       %   [m] distanza della sezione dalla punta del razzo
Ftot = zeros(1,length(dist));           %   inizializzazione vettore carichi
Fmomenti = zeros(1,length(dist));       %   inizializzazione vettore carichi
Fpres = zeros(1,length(dist));          %   inizializzazione vettore carichi

for i = 1:length(dist)
    x = dist(i)-0.35;               %   [m] braccio della forza agente su cp dell'ogiva
    y = 4.072-dist(i);              %   [m] braccio della forza agente su cp delle alette
    M_og = 650*x;                   %   [Nm] Momento della forza agente su cp ogiva
    M_alette = 650*y;               %   [Nm] Momento della forza agente su cp alette
    J = pi*(0.075^4-0.073^4)/4;     %   [m^4] Momento inerzia sezione circolare
    
    thick = 2e-3;                   %   [m] Spessore carbonio serbatoio
    max_y = 0.075;                  %   [m] Massima distanza dal piano neutro

    sigma_max_momenti = abs(max_y/J*(M_og-M_alette));   %   [N/m^2] massima tensione dovuta ai momenti
    % sigma_avg_momenti = sigma_max_momenti/2;          %   [N/m^2] tensione media dovuta ai momenti
    Fmomenti(i) = thick*sigma_max_momenti;              %   [N/m]   Nx
    Area_cil = pi*(0.075^2-0.073^2);                    %   [m^2]   Area sezione circolare
    Fpres(i) = 648000;                                  %   [N/m]
    Ftot(i) = Fmomenti(i)+Fpres(i);          %   [N/m]
    % fprintf('Distanza = %4.3f  [N]\n',dist(i))
    % fprintf('Nmomenti = %4.0f [N/m]\n',Fmomenti)
    % fprintf('Nspinta  = %4.0f   [N/m]\n',Fspinta)
    % fprintf('Ntot     = %4.0f [N/m]\n',Ftot(i))
    % fprintf('----------------------\n')
end

Nfinale = max(Ftot)-Fpres(i);

%%  mega-plot

figure()
plot(dist,Ftot,"DisplayName","Sforzo totale")
hold on
plot(dist(end),Ftot(end),"or")
plot(dist,Fpres,"DisplayName","Sforzo pressione")
plot(dist,Fmomenti,"DisplayName","Sforzo lift")

[~,numero] = min(Fmomenti);
mid = fix(length(Fmomenti)-numero)/2;

pb = pbaspect;
ra = range(xlim)/pb(1);
rb = range(ylim)/pb(2); 
b = Fmomenti(end) - Fmomenti(numero);
a = dist(end)-dist(numero);

theta = rad2deg(atan(  (b/rb) / (a/ra)  ));

h1 = text(dist(end),Fmomenti(end),"Sforzo lift");
set(h1,'Rotation',theta,'VerticalAlignment',"bottom","HorizontalAlignment","right");
h2 = text(dist(mid),Ftot(mid),"Sforzo totale");
set(h2,'Rotation',theta,'VerticalAlignment',"bottom");
title("Carichi assiali agenti sulla skin del serbatoio")
xlabel("Distanza dalla punta [m]")
ylabel("Carico per unità di spessore [N/m]")
xline(dist1,'--k',"LabelHorizontalAlignment","left","Label","testa superiore")
xline(dist2,'--k',{"testa inferiore"})
text(dist(end),Fpres(end),"Sforzo pressione",HorizontalAlignment="right",VerticalAlignment="bottom")
titolino = sprintf("\\underline{%4.0f kN/m}",Ftot(end)/1000);
text(dist(end)-2/100*dist(end),Ftot(end),titolino,HorizontalAlignment="right",Interpreter="latex",Color="r")
hold off
set(gcf,"Color","w")
set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'PaperPositionMode', 'auto')
