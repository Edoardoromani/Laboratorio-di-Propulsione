%%  INITIALIZE
% federico bortolato 2023 - laboratorio di propulsione spaziale
% calcolatore sforzi di taglio giunzioni incollaggi

close all
clear all
clc

%%  Choice of theories to be computed
%   comment the rows to be ignored

theories_list = [   
                    % "VOLKERSEN";        %   shear             - elastic - single lap
                    % "REISSNER";         %   shear and bending - elastic - single lap
                    % % "HART_ELASTO_S1";   %   shear and bending - elastic - single lap #sembra sovrastimare
                    % "HART_ELASTO_S2";   %   shear and bending - elastic - single lap
                    % "HART_PLASTO_S";    %   shear and bending - plastic - single lap
                    "HART_PLASTO_D";    %   shear             - plastic - double lap

                ];

%%  DATI

P = 6e6; 		        % [Pa] pressione in serbatoio
r_int = 0.0687; 		% [m] raggio sferico interno testa 
A_int = pi*r_int^2;     % [m^2] area sezione testa sferica, ovvero l'area 
                        %       su cui la pressione agente genera una
                        %       tensione longitudinale
r_ext = 0.0722; 		    % [m] raggio cilindrico esterno testa 
Circ = 2*pi*r_ext; 	    % [m] circonferenza cilindrica esterna testa 
sigma = 28e6; 		    % [Pa] lap shear strength epoxy resin
FoS = 3;              % Fattore di sicurezza #ECSS
t = 0.0028;              % [m] spessore minore (carbonio)
E = 70e9;               % [Pa] Modulo Young minore (alluminio)
ta = 0.0001;            % [m] spessore adesivo
poisson_coeff = 0.35;    % Poisson alluminio
Ga = 1e9;               % [Pa] adhesive shear modulus 

% Eepx = 10e9;            % young's modulus epoxy
% poisson_epx = 0.35;     % poisson epoxy
% G = Eepx/2/(1+poisson_epx); % shear modulus epoxy

%%  AVERAGE STRESS CALCULATIONS

F_long = A_int*P*FoS; 	% [N] Forza "stappante" 
A_req = F_long/sigma; 	% [m^2] Area necessaria per la tenuta dell'incollaggio
H = A_req/Circ; 	    % [m] altezza minima di incollaggio a partire da saldatura testa

L =  0.05; % [m] incollaggio massimo disponibile
N = F_long/Circ; % [N/m]
tau_avg = N/L; % [Pa]

% fprintf("\nAltezza di incollaggio  :       H = %2.4f mm\n", H*1000)
% fprintf("\nAverage Adhesive Stress : tau_avg = %2.4f  MPa\n", tau_avg/1e6)

%%  DATA HANDLING

npoints = 1000;  %   number of points to plot for each theory

Xs = linspace(0,L/2,npoints);

datazero = zeros(1,length(Xs)*2);

for j = 1:length(theories_list)

    data.(theories_list(j)).name = strings;
    data.(theories_list(j)).color = strings;
    data.(theories_list(j)).max = strings;
    data.(theories_list(j)).taus = datazero;
    data.(theories_list(j)).xss = datazero;

end

%%  THEORIES CALCUATIONS

tau_vec = zeros(1,length(Xs));

for j = 1:length(theories_list)
    theory_name = theories_list(j);
    switch theory_name

        case "VOLKERSEN" % shear - elastic - single lap

            tt = t; tb = t; b = Circ;

            for i = 1:length(Xs)
                y = Xs(i);
                omega = sqrt(Ga/E/tt/ta*(1+tt/tb));
                tau_vec(i) = F_long*omega/2/b*cosh(omega*y)/sinh(omega*L/2)+((tt-tb)/(tt+tb))*(omega*L/2)*sinh(omega*y)/cosh(omega*L/2);
            end

            data.VOLKERSEN.xss = [fliplr(-1*Xs) Xs]*1e3;
            data.VOLKERSEN.taus = [fliplr(tau_vec) tau_vec]/1e6;
            data.VOLKERSEN.name = "Volkersen";
            data.VOLKERSEN.color = "#0072BD";
            data.VOLKERSEN.max = sprintf("\nMax Adhesive Stress Volk: tau_max = %2.4f MPa\n", max(tau_vec)/1e6);

        case "REISSNER" % shear and bending - elastic - single lap

            c = L/2; Pmed = N;
            u = sqrt(3*Pmed/t/E*(1-poisson_coeff^2)/2)/t;
            beta = sqrt(8*Ga*t/ta/E);
            k = (cosh(u*c))/(cosh(u*c)+2*sqrt(2)*sinh(u*c));
            
            for i = 1:length(Xs)
                y = Xs(i);
                tau_vec(i) = 1/8*Pmed/c*( beta*c/t*(1+3*k)*(cosh(beta*c*y/t/c))/(sinh(beta*c/t)) + 3*(1-k));
            end
            tau_max = Pmed/E/t*Ga/ta/beta;
            data.REISSNER.xss = [fliplr(-1*Xs) Xs]*1e3;
            data.REISSNER.taus = [fliplr(tau_vec) tau_vec]/1e6;    
            data.REISSNER.name = "Goland-Reissner";
            data.REISSNER.color = "#D95319";
            data.REISSNER.max = sprintf("\nMax Adhesive Stress Reiss: tau_max = %2.4f MPa\n", max(tau_vec)/1e6);

        case "HART_ELASTO_S1" % shear and bending - elastic - single lap        
            %   ECSS FORMULAS: 10.11.2.1 SINGLE LAP SHEAR JOINT SHEAR STRESS

            ni = poisson_coeff^2; 
            
            for i = 1:length(Xs)
                y = Xs(i);
                alpha = 1/(1+2*sqrt(2)*tanh(L/t*((3*(1-ni)*N)/(2*E*t))));
                W = sqrt( (2*(1-ni)*Ga)/(E*t*ta) );
                K = 1/4*(W*L*(1+3*alpha)*cosh(W*y*2)/sinh(L*W)+3*(1-alpha));
                tau_vec(i) = K*tau_avg;
            end
            
            data.HART_ELASTO_S1.xss = [fliplr(-1*Xs) Xs]*1e3;
            data.HART_ELASTO_S1.taus = [fliplr(tau_vec) tau_vec]/1e6;
            data.HART_ELASTO_S1.name = "Hart-elastic-s1";
            data.HART_ELASTO_S1.color = "#4DBEEE";
            data.HART_ELASTO_S1.max = sprintf("\nMax Adhesive Stress Hart: tau_max = %2.4f MPa\n", max(tau_vec)/1e6);
           
        case "HART_ELASTO_S2" % shear and bending - elastic - single lap

            c = L/2; Pmed = N; 
            lamda2 = sqrt(((1+3*(1-poisson_coeff^2))/4)*2*Ga/ta/E/t); 
            D = E*t^3/12/(1-poisson_coeff^2);
            csi = sqrt(Pmed/D);
            M = Pmed*(t+ta)/2/(1+csi*c+csi^2*c^2/6);
            A2 = Ga/ta/E/t*(Pmed+6*M/t*(1-poisson_coeff^2))/(2*lamda2*sinh(2*lamda2*c));
            C2 = (Pmed - A2/lamda2*sinh(2*lamda2*c))/2/c;
            
            for i = 1:length(Xs)
                y = Xs(i);
                tau_vec(i) = A2*cosh(2*lamda2*y) + C2;
            end

            data.HART_ELASTO_S2.xss = [fliplr(-1*Xs) Xs]*1e3;
            data.HART_ELASTO_S2.taus = [fliplr(tau_vec) tau_vec]/1e6;
            data.HART_ELASTO_S2.name = "Hart-elastic";
            data.HART_ELASTO_S2.color = "#7E2F8E";
            data.HART_ELASTO_S2.max = sprintf("\nMax Adhesive Stress Har2: tau_max = %2.4f MPa\n", max(tau_vec)/1e6);

        % case "HART_PLASTO_Sxxx" % shear and bending - plastic - single lap #NONFUNZIA

            % c = L/2; Pmed = N; 
            % D = E*t^3/12/(1-poisson_coeff^2);
            % csi = sqrt(Pmed/D);
            % k = 1/(1+csi*c+csi^2*c^2/6);
            % tau_p = sigma;
            % lamda = sqrt(2*Ga/t/ta/E);
            % lamda2 = sqrt((1+3*(1-poisson_coeff^2))/4)*lamda; 
            % gamma_e = tau_p/Ga;
            % 
            % %   NUMERICAL OPTIMIZATION
            % % f1 = @(d,K,ssr) Pmed/L/tau_p*lamda2*L - 2*lamda2*(L-d)/2 - (1-K)*lamda2*d - K*tanh(lamda2*d);
            % % f2 = @(d,K,ssr) (1+3/(1+csi*c+csi^2*c^2/6)*(1-poisson_coeff^2)*(1+ta/t))*Pmed/tau_p*lamda^2*(L-d)/2 ...
            % % -2*ssr-K*(2*lamda2*(L-d)/2)^2;
            % % f3 = @(d,K,ssr) 2*ssr-K*( (2*lamda2*(L-d)/2+tanh(lamda2*d))^2 -(tanh(lamda2*d))^2);
            % % 
            % % F = @(par) [f1(par(1),par(2),par(3)); f2(par(1),par(2),par(3)); f3(par(1),par(2),par(3));];
            % % [par] = fsolve(F,[0; 0;0]);
            % % 
            % % d = par(1);
            % % K = par(2);
            % % ssr = par(3);
            % 
            % %   SYMBOLIC CALCULATION - probably slower
            % syms d K gamma_p
            % f1 =  tau_avg/tau_p*lamda2*L - 2*lamda2*(L-d)/2 - (1-K)*lamda2*d - K*tanh(lamda2*d);
            % f2 =  (1+3*k*(1-poisson_coeff^2)*(1+ta/t))*Pmed/tau_p*lamda^2*(L-d)/2 ...
            % -2*gamma_p/gamma_e-K*(2*lamda2*(L-d)/2)^2;
            % f3 =  2*gamma_p/gamma_e-K*( (2*lamda2*(L-d)/2+tanh(lamda2*d))^2 -(tanh(lamda2*d))^2);
            % [par] = solve(f1,f2,f3);
            % 
            % d = par.d;
            % K = par.K;
            % gamma_p = par.gamma_p;
            % 
            % A2 = K*tau_p/cosh(lamda2*d);
            % 
            % for i = 1:length(Xs)
            % 
            %     y = Xs(i);
            %     if y < d/2
            %         tau_vec(i) = A2*cosh(2*lamda2*y) + tau_p*(1-K);
            %     else
            %         tau_vec(i) = tau_p;
            %     end
            % end
            % 
            % data.HART_PLASTO_S.xss = [fliplr(-1*Xs) Xs]*1e3;
            % data.HART_PLASTO_S.taus = [fliplr(tau_vec) tau_vec]/1e6;
            % data.HART_PLASTO_S.name = "Hart-plastic-s";
            % data.HART_PLASTO_S.color = "#77AC30";
            % data.HART_PLASTO_S.max = sprintf("\nMax Adhesive Stress Har3: tau_max = %2.4f MPa\n", max(tau_vec)/1e6);

        case "HART_PLASTO_S" % shear and bending - plastic - single lap

            c = L/2; Pmed = N; 
            D = E*t^3/12/(1-poisson_coeff^2);
            csi = sqrt(Pmed/D);
            k = 1/(1+csi*c+csi^2*c^2/6);
            tau_p = sigma;
            lamda = sqrt(2*Ga/t/ta/E);
            lamda2 = sqrt((1+3*(1-poisson_coeff^2))/4)*lamda; 
            gamma_e = tau_p/Ga;

            %   NUMERICAL OPTIMIZATION

            f1 = @(A2,C2,A3,B3,d) A2/Ga*cosh(lamda*d) + C2/Ga - gamma_e;
            f2 = @(A2,C2,A3,B3,d) A2/Ga*2*sinh(lamda*d) - B3;
            f3 = @(A2,C2,A3,B3,d) A2/Ga*4*lamda^2*cosh(lamda*d) - 2*A3;
            f4 = @(A2,C2,A3,B3,d) 2*A3*(L-d)/2 + B3 - Pmed/E/t/ta*(1+3*k*(1-poisson_coeff^2)*(1+ta/t));
            f5 = @(A2,C2,A3,B3,d) A2/2/lamda*sinh(lamda*d) + C2*d/2 + tau_p*(L-d)/2 - Pmed/2;

            F = @(par) [

                f1(par(1),par(2),par(3),par(4),par(5));
                f2(par(1),par(2),par(3),par(4),par(5));
                f3(par(1),par(2),par(3),par(4),par(5));
                f4(par(1),par(2),par(3),par(4),par(5));
                f5(par(1),par(2),par(3),par(4),par(5))

                ];

            options = optimoptions('fsolve','Display','iter');
            options = optimoptions(options,MaxIterations=1e4,MaxFunctionEvaluations=1e5);

            [par] = fsolve(F,[0;0;0;0;0],options);

            d = par(5);
            A2 = par(1);
            C2 = par(2);
            K = A2/tau_p*cosh(lamda2*d);
            A3 = par(3);
            B3 = par(4);

            %   SYMBOLIC CALCULATION - slower

            % syms A2 C2 A3 B3 d 
            % 
            % f1 = A2/Ga*cosh(lamda*d) + C2/Ga - gamma_e;
            % f2 = A2/Ga*2*sinh(lamda*d) - B3;
            % f3 = A2/Ga*4*lamda^2*cosh(lamda*d) - 2*A3;
            % f4 = 2*A3*(L-d)/2 + B3 - Pmed/E/t/ta*(1+3*k*(1-poisson_coeff^2)*(1+ta/t));
            % f5 = A2/2/lamda*sinh(lamda*d) + C2*d/2 + tau_p*(L-d)/2 - Pmed/2;
            % [par] = vpasolve(f1,f2,f3,f4,f5);
            % 
            % d = par.d;
            % A2 = par.A2;
            % C2 = par.C2;

            for i = 1:length(Xs)
                
                y = Xs(i);
                if y < d/2
                    tau_vec(i) = A2*cosh(2*lamda*y) + C2;
                else
                    tau_vec(i) = tau_p;
                end
            end

            data.HART_PLASTO_S.xss = [fliplr(-1*Xs) Xs]*1e3;
            data.HART_PLASTO_S.taus = [fliplr(tau_vec) tau_vec]/1e6;
            data.HART_PLASTO_S.name = "Hart-plastic-single";
            data.HART_PLASTO_S.color = "#77AC30";
            data.HART_PLASTO_S.max = sprintf("\nMax Adhesive Stress Har3: tau_max = %2.4f MPa\n", max(tau_vec)/1e6);

        case "HART_PLASTO_D" % shear - plastic - double lap
            tau_p = sigma;
            lamda = sqrt(2*Ga/t/ta/E);
            lamda2 = sqrt(((1+3*(1-poisson_coeff^2))/4))*lamda; 
            
            %   NUMERICAL OPTIMIZATION 
            f1 =  @(d,ssr) (lamda*(L-d)/2 +tanh(lamda*d/2))^2 - (tanh(lamda*d/2))^2 - 2*ssr;
            f2 =  @(d,ssr) tau_avg/tau_p - (1-d/L) - (tanh(lamda*d/2))/(lamda*L/2); 
            F = @(par) [f1(par(1),par(2));f2(par(1),par(2))];
            [par] = fsolve(F,[0,0]);
            
            d = par(1);
            ssr = par(2);
            
            %   SYMBOLIC CALCULATION - probably slower
            % syms d ssr
            % f1 =  (lamda*(L-d)/2 +tanh(lamda*d/2))^2 - (tanh(lamda*d/2))^2 - 2*ssr;
            % f2 =  tau_avg/tau_p - (1-d/L) - (tanh(lamda*d/2))/(lamda*L/2); 
            % 
            % [par] = solve(f1,f2)
            % 
            % d = par.d;
            % ssr = par.ssr;
            
            A = tau_p/cosh(lamda*d/2);
            
            for i = 1:length(Xs)
                
                y = Xs(i);
                if y < d/2
                    tau_max = A*cosh(lamda*y);
                else
                    tau_max = tau_p;
                end
                tau_vec(i) = tau_max;
            end

            data.HART_PLASTO_D.xss = [fliplr(-1*Xs) Xs]*1e3;
            data.HART_PLASTO_D.taus = [fliplr(tau_vec) tau_vec]/1e6;
            data.HART_PLASTO_D.name = "Hart-plastic-double";
            data.HART_PLASTO_D.color = "#EDB120";
            data.HART_PLASTO_D.max = sprintf("\nMax Adhesive Stress Har4: tau_max = %2.4f MPa\n", tau_max/1e6);
    end
end

%% COMPARATIVE PLOT

figure(Name = "stress comparison plot", NumberTitle="off");

for index = 1:length(theories_list)
    H = plot(data.(theories_list{index}).xss,data.(theories_list{index}).taus, ...
        LineWidth=2,DisplayName=data.(theories_list{index}).name,Color=data.(theories_list{index}).color);
    if index == 1
        hold on
    end
end

extra1 = yline(tau_avg/1e6, Color = 'r', LineWidth=1,LineStyle='--',DisplayName="tau medio");
extra2 = xline(0, Color = 'k', LineWidth=0.5,LineStyle='-',DisplayName="mezzeria");
title("Stress Distribution in Single Lap Bonded Joint - COMPARISON OF THEORIES")
xlabel("length [mm]")
ylabel("shear stress [MPa]")
legend()
hold off
set(gcf,"Color","w")
set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'PaperPositionMode', 'auto')

%%  extra: shear strains single lap hart

% figure(Name = "shear plot", NumberTitle="off");
% 
% gamma_p= @(x) A3.*x.^2 + B3.*x + gamma_e;
% X1 = Xs(Xs<d/2);
% X2 = Xs(Xs>=d/2); X2g = X2 - X2(1);
% gammas_e = tau_vec(Xs<d/2)/Ga;
% gammas_p = gamma_p(X2g);
% plot(X1,gammas_e)
% hold on
% plot(X2,gammas_p)
% title("Shear Strain in Hart-Smith Elasto-Plastic Analysis of Single Lap Joint")
% ylabel("Deformazione a taglio")
% xlabel("Distanza dal centro dell'incollaggio [m]")
% xline(d/2 ,Label="Tratto Elastico",LabelHorizontalAlignment="left",LabelOrientation="horizontal",LineStyle="-.",Color="k")
% xline(d/2,Label="Tratto Plastico",LabelHorizontalAlignment="right",LabelVerticalAlignment="bottom",LineStyle="-.",Color="k")
% % yline(gamma_e,Label="Zona Plastica",LabelHorizontalAlignment="left",LabelVerticalAlignment="top",LineStyle="-.",Color="k")
% % yline(gamma_e,Label="Zona Elastica",LabelHorizontalAlignment="left",LabelVerticalAlignment="bottom",LineStyle="-.",Color="k")
% hold off
% set(gcf,"Color","w")
% set(gcf, 'InvertHardcopy', 'off');
% set(gcf, 'PaperPositionMode', 'auto')

%%  extra: shear strains double lap hart

figure(Name = "shear plot double", NumberTitle="off");
gamma_e = tau_p/Ga;
gamma_p1 = ssr*gamma_e;

gamma_p= @(x) lamda^2/2/Ga*tau_p.*x.^2 + lamda*tau_p/Ga*tanh(lamda*d/2).*x + tau_p/Ga;
X1 = Xs(Xs<d/2);
X2 = Xs(Xs>=d/2);X2g = X2 - X2(1);
gammas_e = tau_vec(Xs<d/2)/Ga;
gammas_p = gamma_p(X2g);
plot(X1,gammas_e)
hold on
plot(X2,gammas_p)
title("Shear Strain in Hart-Smith Elasto-Plastic Analysis of Double Lap Joint")
ylabel("Deformazione a taglio")
xlabel("Distanza dal centro dell'incollaggio [m]")
xline(d/2 ,Label="Tratto Elastico",LabelHorizontalAlignment="left",LabelOrientation="horizontal",LineStyle="-.",Color="k")
xline(d/2,Label="Tratto Plastico",LabelHorizontalAlignment="right",LabelVerticalAlignment="bottom",LineStyle="-.",Color="k")
% yline(gamma_e+gamma_p1,Label="Zona Plastica",LabelHorizontalAlignment="left",LabelVerticalAlignment="bottom",LineStyle="-.",Color="k")
% yline(gamma_e,Label="Zona Elastica",LabelHorizontalAlignment="left",LabelVerticalAlignment="bottom",LineStyle="-.",Color="k")

hold off
set(gcf,"Color","w")
set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'PaperPositionMode', 'auto')
