function [Nx_max,Ny_max] = COMPOSITE_CALCULATOR(U1_m,U2_m,U3_m,U4_m,sigma1t_ult_m,sigma1c_ult_m,sigma2t_ult_m,sigma2c_ult_m,tau12_ult_m,n,t_ply,teta_ply,material,Loads)

% four invariants
for k=1:1:n
    U1(k)=U1_m(material(k));
    U2(k)=U2_m(material(k));
    U3(k)=U3_m(material(k));
    U4(k)=U4_m(material(k));
    sigma1t_ult(k)=sigma1t_ult_m(material(k));
    sigma1c_ult(k)=sigma1c_ult_m(material(k));
    sigma2t_ult(k)=sigma2t_ult_m(material(k));
    sigma2c_ult(k)=sigma2c_ult_m(material(k));
    tau12_ult(k)=tau12_ult_m(material(k));
end

hl=0;
%da riferimento locale a globale
for k=1:1:n
    t(k)=t_ply(k)/1000;
    teta=pi*teta_ply(k)/180;
    Q(1,1,k)=U1(k)+U2(k)*cos(2*teta)+U3(k)*cos(4*teta);
    Q(1,2,k)=U4(k)-U3(k)*cos(4*teta);
    Q(2,1,k)=Q(1,2,k);
    Q(2,2,k)=U1(k)-U2(k)*cos(2*teta)+U3(k)*cos(4*teta);
    Q(1,3,k)=0.5*U2(k)*sin(2*teta)+U3(k)*sin(4*teta);
    Q(3,1,k)=Q(1,3,k);
    Q(2,3,k)=0.5*U2(k)*sin(2*teta)-U3(k)*sin(4*teta);
    Q(3,2,k)=Q(2,3,k);
    Q(3,3,k)=0.5*(U1(k)-U4(k))-U3(k)*cos(4*teta);
    s=sin(teta);
    c=cos(teta);
    T(:,:,k)=[c^2 s^2 2*s*c; s^2 c^2 -2*s*c; -s*c s*c c^2-s^2];
    hl=hl+t(k);  %spessore laminato
end

h0=-hl/2; %posizione piano medio
sum=0;
for k=1:1:n
    sum=sum+t(k);
    h(k)=-hl/2+sum;
end

%Stiffness Matrix
for i=1:1:3
    for j=i:1:3
        sumA=Q(i,j,1)*(h(1)-h0);
        sumB=Q(i,j,1)*(h(1)^2-h0^2);
        sumD=Q(i,j,1)*(h(1)^3-h0^3);
        for k=2:1:n
            sumA=sumA+Q(i,j,k)*(h(k)-h(k-1));
            sumB=sumB+Q(i,j,k)*(h(k)^2-h(k-1)^2);
            sumD=sumD+Q(i,j,k)*(h(k)^3-h(k-1)^3);
        end
        A(i,j)=sumA;
        A(j,i)=sumA;
        B(i,j)=sumB/2;
        B(j,i)=sumB/2;
        D(i,j)=sumD/3;
        D(j,i)=sumD/3;
    end
end

% A  %matrice di rigidezza assiale
% B  %matrice di rigidezza di accoppiamento
% D  %matrice di rigidezza flessionale

% Stiffness Matrix
SM = [A B; B D];
% eps=inv(SM)*LOAD';
eps = SM\Loads';

%Global Strains
epsx(1)=eps(1)+h0*eps(4);
epsy(1)=eps(2)+h0*eps(5);
gammaxy(1)=eps(3)+h0*eps(6);
for k=2:1:n+1
epsx(k)=eps(1)+h(k-1)*eps(4);
epsy(k)=eps(2)+h(k-1)*eps(5);
gammaxy(k)=eps(3)+h(k-1)*eps(6);
end

%Global and Local Stresses
for k=1:1:n
sigmatop(k,:)=Q(:,:,k)*[epsx(k) epsy(k) gammaxy(k)]';
sigmabottom(k,:)=Q(:,:,k)*[epsx(k+1) epsy(k+1) gammaxy(k+1)]';
sigmatop_loc(k,:)=T(:,:,k)*sigmatop(k,:)';
sigmabottom_loc(k,:)=T(:,:,k)*sigmabottom(k,:)';

if sigmatop_loc(k,1)>0
    X1=sigma1t_ult(k);
else 
    X1=sigma1c_ult(k);
end
if sigmatop_loc(k,2)>0
    X2=sigma1t_ult(k);
    Y=sigma2t_ult(k);
else
    X2=sigma1c_ult(k);
    Y=sigma2c_ult(k);
end
SRtop(k)=1/sqrt((sigmatop_loc(k,1)/X1)^2+(sigmatop_loc(k,2)/Y)^2+(sigmatop_loc(k,3)/tau12_ult(k))^2-(sigmatop_loc(k,1)/X2*sigmatop_loc(k,2)/X2));
% criterio di rottura Tsai-Hill faccia superiore della lamina

if sigmabottom_loc(k,1)>0
    X1=sigma1t_ult(k);
else 
    X1=sigma1c_ult(k);
end
if sigmabottom_loc(k,2)>0
    X2=sigma1t_ult(k);
    Y=sigma2t_ult(k);
else
    X2=sigma1c_ult(k);
    Y=sigma2c_ult(k);
end
SRbottom(k)=1/sqrt((sigmabottom_loc(k,1)/X1)^2+(sigmabottom_loc(k,2)/Y)^2+(sigmabottom_loc(k,3)/tau12_ult(k))^2-(sigmabottom_loc(k,1)/X2*sigmabottom_loc(k,2)/X2));
% criterio di rottura Tsai-Hill faccia inferiore della lamina
end

[SRmin,n_ply]=min([SRtop SRbottom]);

if (n_ply/2-floor(n_ply/2))==0
    n_ply=n_ply/2;
else
    n_ply=floor(n_ply/2)+1;
end

Nx_max = Loads(1)*SRmin;
Ny_max = Loads(2)*SRmin; 
Nxy_max = Loads(3)*SRmin; 

Mx_max=Loads(4)*SRmin;
% disp('[Nm/m]')
My_max=Loads(5)*SRmin;
% disp('[Nm/m]')
Mxy_max=Loads(6)*SRmin;  
% disp('[Nm/m]')

%%  SCONTRINO RISULTATI

fprintf('First ply to fail: %d\n',n_ply)
fprintf('----------------------\n')
fprintf('Strenght Ratio = %4.3f\n',SRmin)
fprintf('----------------------\n')
fprintf("Active loads:\n")
fprintf('Nxxx  = %10.0f [N/m]\n',Loads(1))
fprintf('Nyyy  = %10.0f [N/m]\n',Loads(2))
fprintf('----------------------\n')
fprintf('Max loads:\n')
fprintf('Nxmax = %10.0f [N/m]\n',Nx_max)
fprintf('Nymax = %10.0f [N/m]\n',Ny_max)
fprintf('----------------------\n')

for i=1:1:n-1
    pos(2*i)=h(i);
    pos(2*i-1)=h(i);
    sigma_x(2*i-1)=sigmabottom(i,1);
    sigma_x(2*i)=sigmatop(i+1,1);
    sigma_y(2*i-1)=sigmabottom(i,2);
    sigma_y(2*i)=sigmatop(i+1,2);
    tau_xy(2*i-1)=sigmabottom(i,3);
    tau_xy(2*i)=sigmatop(i+1,3);    
end

%%  PLOT RISULTATI

pos=-[h0 pos -h0]*1000;
sigma_x=[sigmatop(1,1) sigma_x sigmabottom(n,1)]*10^-6;
sigma_y=[sigmatop(1,2) sigma_y sigmabottom(n,2)]*10^-6;
tau_xy=[sigmatop(1,3) tau_xy sigmabottom(n,3)]*10^-6;

% f1 = figure;
% plot(sigma_x,pos)
% xlabel('tensioni sigma_x  [MPa]')
% ylabel('h  [mm]')
% 
% f2 = figure;
% plot(sigma_y,pos)
% xlabel('tensioni sigma_y  [MPa]')
% ylabel('h  [mm]')
% 
% f3 = figure;
% plot(tau_xy,pos)
% xlabel('tensioni tau  [MPa]')
% ylabel('h  [mm]')
% 
% h=[h0 h];
% 
% f4 = figure;
% plot(epsx,-h*1000)
% xlabel('deformazioni eps_x')
% ylabel('h  [mm]')
% 
% f5 = figure;
% plot(epsy,-h*1000)
% xlabel('deformazioni eps_y')
% ylabel('h  [mm]')
% 
% f6 = figure;
% plot(gammaxy,-h*1000)
% xlabel('deformazioni gamma xy')
% ylabel('h  [mm]')
% 
% %%  spread figs
% 
% movegui(f1,"northwest")
% movegui(f2,"north")
% movegui(f3,"northeast")
% movegui(f4,"southwest")
% movegui(f5,"south")
% movegui(f6,"southeast")

%% PRINT USEFUL MAX RESULTS TO WINDOW

% fig1 = figure();
% width = 500;
% x = 1920/2 - width/2;
% height = 50;
% y = 1080/2 - 2*height;
% 
% fig1.Position = [x y width height];
% textfinal = sprintf("Nx       = %6.0f [N/m] Ny     = %6.0f [N/m]\n" + ...
%                     "Nxmax = %6.0f [N/m] Nymax = %6.0f [N/m]",Loads(1),Loads(2),Nx_max,Ny_max);
% uicontrol('style','text',"FontSize",15,"HorizontalAlignment","left",...
%               'units','norm',...
%               'pos',[0 0 1 1],...
%               'string',textfinal,...
%               'min',0,'max',10);