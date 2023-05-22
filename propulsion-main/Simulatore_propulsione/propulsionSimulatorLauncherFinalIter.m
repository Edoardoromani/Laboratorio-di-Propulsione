%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  _____  __  __  ____  _   _  ____  _____                %
%                 |_   _||  ||  || __ \| | | |/ ___||_   _|               %
%                   | |  |  __  ||    /| |_| |\___ \  | |                 %
%                   |_|  |__||__||_|\_\ \___/ |____/  |_|                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Project: Propulsion Simulator

Property of THRUST, unauthorized distribution is not allowed

Description:
  This code simulates the whole propulsion system of a hybrid rocket based
  on N2O and paraffin as propellants. It is equivalent to a static fire
  test. Specifically it is the wrapper of PropulsionSimulator.m which is
  the actual function that does all the calculations

Changelog:
  > version: 1.5 - 25/03/2023 - Giulia Quinzi & Riccardo Dei Rossi
     - added plot (figure 6) about mDot_ox depending on deltaP_inj
  > version: 1.5 - 14/03/2023 - Alessandro Rampazzo & Giulia Quinzi & Riccardo Dei Rossi & Luca
    % Dituri
    - added plots (figures 4,5) with some useful information about the
    simulation

  > version: 1.3 - 05/02/2023 - Alessandro Rampazzo
    - changed batch execution
    - created function saveDynamicFile that saves the output in a dynamic
    named .mat file

  > version: 1.2 - 27/01/2023 - Alessandro Rampazzo
    - changed tSpan definition
    - bugfix the calculation of the mean with the integral mean
    - added geometry specification
    - refactoring 

  > version: 1.1 - 04/12/2022 - Alessandro Rampazzo
    - added a new combustionFunction to the buildRocket nested function to
    simulate combustion
    - added new info outputs on the command window
    - fixed some plots
    - changed name to the tank function to tankFunction
    - fixed the nominal tb automatic identification
    - added the rocket.tank.k field

  > version: 1.0 - 27/11/2022 - Alessandro Rampazzo
    - added heading, description, and changelog

  > version: 0.1 - 15/03/2022 - Alessandro Rampazzo
    - created
%}

% some standard commands
clc
clear
close all

format compact

% adds the paths to the files needed for the simulation 
addpath("Funzioni sottosistemi");
% addpath("Funzioni sottosistemi\Proprietà_N2O");
addpath("../Chimica/NFP")
% addpath("Matlab Functions")
addpath("../Matlab Functions")


%% options

opts.saveFile = false;             % save the simulation output 
opts.plots = true;                 % activate plots
opts.printInfo = true;            % print some informations regarding the simulation
opts.batchCalculation = false;      % group calculation in batches
opts.plotInfo = true;             % plots other informations like mean, max and min of some variables

opts.folderName = "Engine Data";    % nome della cartella in cui salvare i dati    
opts.fileName = "TH-1000";           % nome iniziale di una serie di dati

% time integration vector
opts.tSpan = [0,30]; % [s]

% propulsion model type - "simplified" or "unsteady"
opts.propulsionModel = "unsteady";

%% environment
env.g = 9.81;   % [m/s^2]
env.h = 0;      % [m] altitudine di test 
[~,~,env.p_a,~] = atmosisa(env.h);  % [Pa] pressione all'altitudine operativa 

%% tank
rocket.tank.V = 0.031;        % [m^3] volume serbatoio
rocket.tank.T0 = 25 + 273.15; % [°C] temperatura iniziale (nominal 25°C)
rocket.tank.m0 = 20;          % [kg] massa iniziale di ossidante da caricare nel tank
                              % (valore oppure "max" per ma massima massa caricabile)
rocket.tank.k = 1.33;         % coefficiente di espansione adiabatica dell'N2O
rocket.tank.f_min = 0.05;     % percentuale di vapore minima in volume

%% injection
rocket.inj.Cd = 0.72;         % discharge coefficient dell'iniettore 
rocket.inj.n_holes = 50;       % numero di fori iniettore 9
rocket.inj.D_holes = 1.2e-3;  % [m] diametro dei fori 1.5
rocket.inj.p_loss = 5e+5;   % [Pa] perdite nei condotti fluidici
rocket.inj.model = "DYER-ISOH0";
% chose between: SPI (non funziona, bisogna usare l'SPI-ISOH0) SPI-ISOH0, DYER-ISOH0 (il modello Dyer è usato per i test), OMEGA-ISOH0

%% combustion chamber

% injector slot geometry (cylinder)
rocket.cc.D_inj_slot = 33.75 * 1e-3; % [m] -> [mm] diametro slot iniettore
rocket.cc.L_inj_slot = 12.5 * 1e-3; % [m] -> [mm] lunghezza slot iniettore
% prechamber geometry (cylinder)
rocket.cc.D_prechamber = 142* 1e-3; % [m] -> [mm] diametro precamera
rocket.cc.L_prechamber = 37.5 * 1e-3; % [m] -> [mm] lunghezza precamera
% grain geometry (cylinder)
rocket.cc.L_grain = 265 * 1e-3; % [mm] -> [m] lunghezza grano
rocket.cc.D_port0 = 51 * 1e-3;  % [mm] -> [m] diametro iniziale interno del grano
rocket.cc.D_ext = 142 * 1e-3; % [mm] -> [m] diametro esterno del grano
% diaphragm geometry (cylinder)
rocket.cc.D_diaphragm = 22.5 * 1e-3; % [m] -> [mm] diametro diaframma
rocket.cc.L_diaphragm = 15 * 1e-3; % [m] -> [mm] lunghezza diagramma
% postchamber geometry (cylinder)
rocket.cc.D_postchamber = 142 * 1e-3; % [m] -> [mm] diametro postcamera
rocket.cc.L_postchamber = 37.5 * 1e-3; % [m] -> [mm] lunghezza postcamera

rocket.cc.a = 0.132e-3;    % [m/s] coefficiente della regression rate
rocket.cc.n = 0.555;       % esponente della regression rate
rocket.cc.rho_fuel = 924;  % [kg/m^3] densità fuel
rocket.cc.cstar_th = 1600; % [m/s] star teorico
rocket.cc.eta_cstar = 0.9; % efficienza sul cstar
rocket.cc.T_fuel = 1000; % [K] temperatura superficiale del fuel
rocket.cc.cstar = rocket.cc.cstar_th*rocket.cc.eta_cstar;   % [m/s] cstar reale
rocket.cc.T = 2750;       % [K] temperatura di combustione con 90% di 
                          % efficienza sul cstar
rocket.cc.gamma = 1.22;   % costante di espansione adiabatica (Cp/Cv) a 90%
                          % di efficienza
rocket.cc.MM = 25.680;    % [g/mol] massa molecolare dei prodotti di 
                          % combustione a 90% di efficienza


%% nozzle
% nozzle inhibitor geometry (cylinder)
rocket.nozzle.L_inh = 5.5 * 1e-3; % [m] -> [mm] spessore inibitore ugello
rocket.nozzle.D_inh = 43.75 * 1e-3; % [m] -> [mm] diametro inibitore ugello
% nozzle convergent geometry (truncated cone)
rocket.nozzle.D_convergent = 43.75 * 1e-3;% [mm] -> [m] diametro convergente
rocket.nozzle.L_convergent = 22.5 * 1e-3;% [mm] -> [m] lunghezza convergente
% nozzle throat geometry (cylinder)
rocket.nozzle.D_throat0 = 25 * 1e-3; % [mm] -> [m] diametro iniziale di gola
rocket.nozzle.L_throat = 9.75 * 1e-3; % [mm] -> [m] lunghezza di gola
% nozzle divergent geometry (truncated cone)
rocket.nozzle.eps0 = 6.65; % rapporto di espansione iniziale
rocket.nozzle.L_divergent = 40 * 1e-3; % [mm] -> [m] lunghezza divergente

rocket.nozzle.eta_Cf = 0.97;        % efficienza sul Cf
rocket.nozzle.rho_graphite = 1820;  % [kg/m^3] densità grafite dell'ugello
rocket.nozzle.onsetTime = 0;        % [s] tempo che la gola ci mette a
                                    % raggiungere la temperatura di ablazione
rocket.nozzle.erosionModel = -1;    % modello di erosione da usare (-1 per
                                    % non considerare l'erosione)

if opts.batchCalculation
    % calcolo a gruppi

    % variabili da testare - batches - usare solo celle
    
    batches.nozzle_L_inh = {5.5e-3};
    batches.tank_m0 = num2cell([4,4.5]);
    batches.tank_m0{end+1} = "max";
    batches.tank_T0 = num2cell((20:5:25) + 273.15);

    
    batch_fields = fieldnames(batches);
    batch = cell(size(batch_fields));
    for i = 1:length(batch_fields)
        batch{i} = batches.(batch_fields{i});
    end

    % main loop
    for i = 1:length(batch{1})
        for j = 1:length(batch{2})
            for k = 1:length(batch{3})
                
                % stampo informazioni sul calcolo
                fprintf("running batch %i/%i\n",...
                    1 + (k-1) + (j-1)*length(batch{3}) + (i-1)*length(batch{3})*length(batch{2}), ...
                    length(batch{1})*length(batch{2})*length(batch{3}));
                
                % assegno le variabili alla struttura principale: "rocket"
                eval("rocket." + replace(batch_fields{1},"_",".")+"= batch{1}{i};")
                eval("rocket." + replace(batch_fields{2},"_",".")+"= batch{2}{j};")
                eval("rocket." + replace(batch_fields{3},"_",".")+"= batch{3}{k};")
                
                
                % "costruisco" il razzo (calcolo variabili secondarie tipo
                % area della gola, volume del tank, costante del gas R,
                % definisco funzioni etc.)
                % P.S. DA FAR PARTIRE SEMPRE PRIMA DI RUNNARE LA SIMULAZIONE
                try
                    rocket = buildRocket(rocket,env);
                    % runno la simulazione
                    out = singleSimulation(rocket,env,opts);
                catch ME
                    if ME.message == "L'ossidante NON è in condizioni di saturazione!" || ...
                            ME.message == "Non c'è abbastanza volume per alloggiare la massa richiesta alla temperatura data"

                        warning(ME.identifier,"skipping batch because this error occured: %s",ME.message)

                        continue
                    else
                        rethrow(ME)
                    end
                end

                if opts.saveFile
                    saveDynamicFile(out,opts,batch_fields,batch,{batch{1}{i},batch{2}{j},batch{3}{k}})
                end
            end
        end
    end
    
    fprintf("finished\n")
else
    % calcolo singolo

    % "costruisco" il razzo (calcolo variabili secondarie tipo
    % area della gola, volume del tank, costante del gas R,
    % definisco funzioni etc.)
    % P.S. DA FAR PARTIRE SEMPRE PRIMA DI RUNNARE LA SIMULAZIONE
    rocket = buildRocket(rocket,env);

    % runno la simulazione
    out = singleSimulation(rocket,env,opts);

    if opts.saveFile
        % creo la cartella se non esiste
        if ~exist(opts.folderName,"dir")
            mkdir(opts.folderName)
        end
        % costruisco il nome del file da salvare
        name = opts.folderName+"\"+opts.fileName + ".mat";
        % salvo
        save(name,"out");
    end
end


%% functions

function out = singleSimulation(rocket,env,opts)
    
    % runno la simulazione
    if opts.propulsionModel == "simplified"
        out = propulsionSimulator(rocket,env,opts.tSpan);
    elseif opts.propulsionModel == "unsteady"
        out = propulsionSimulator2(rocket,env,opts.tSpan);
    else
        error("propulsion model not understood")
    end

    % calcolo qualche altra variabile

    % flusso di massa [kg/sm^2]
    out.G = out.mDot_ox./(pi * out.r_port.^2);
    
    out.nominal_tb = out.t(find(out.x_tank == 1,1));

    % tempo di burning totale
    out.tb_tot = out.t(end); % [s]

    S_func = @(x) interp1(out.t,out.S,x,'linear');
    out.I_tot = integral(S_func,0,out.tb_tot); % [Ns]
    out.I_tot_nominal = integral(S_func,0,out.nominal_tb); % [Ns]
     
    % stampa qualche informazione
    out.Idx_nominal = find(out.t < out.nominal_tb);

    if opts.printInfo
        fprintf("Info: \n")
        fprintf(" - Nominal tb: %1.2f s\n",out.nominal_tb)
        fprintf(" - Total burning time: %1.2f s\n",out.tb_tot)
        fprintf(" - Initial G: %1.2f kg/sm^2\n",out.G(1))
        fprintf(" - Pre-tailoff G: %1.2f kg/sm^2\n",out.G(out.Idx_nominal(end)))
        fprintf(" - Total nominal impulse is: %1.1f Ns\n",out.I_tot_nominal)
        fprintf(" - Total impulse is: %1.1f Ns\n",out.I_tot)
        fprintf(" - Average nominal thrust is: %1.1f N\n",integral_mean(out.t(out.Idx_nominal),out.S(out.Idx_nominal)))
        fprintf(" - Average nominal O/F: %1.3f\n",integral_mean(out.t(out.Idx_nominal),out.O_F(out.Idx_nominal)))
        fprintf(" - Average nominal regression rate: %1.2f mm/s\n",1e3*integral_mean(out.t(out.Idx_nominal),out.rDot_port(out.Idx_nominal)))
        fprintf(" - Average nominal ox mass flow rate: %1.3f kg/s\n",integral_mean(out.t(out.Idx_nominal),out.mDot_ox(out.Idx_nominal)))
        fprintf(" - Average nominal fuel mass flow rate: %1.3f kg/s\n",integral_mean(out.t(out.Idx_nominal),out.mDot_fuel(out.Idx_nominal)))
        fprintf(" - Average nominal chamber pressure: %1.2f bar\n",1e-5*integral_mean(out.t(out.Idx_nominal),out.p_cc(out.Idx_nominal)))
        fprintf(" - Starting tank pressure: %1.2f bar\n",1e-5*out.p_tank(1))
        fprintf(" - Average nominal tank pressure: %1.2f bar\n",1e-5*integral_mean(out.t(out.Idx_nominal),out.p_tank(out.Idx_nominal)))
        fprintf(" - Starting tank temperature: %1.1f °C\n",out.T_tank(1)-273.15)
        fprintf(" - Average nominal tank temperature: %1.1f °C\n",integral_mean(out.t(out.Idx_nominal),out.T_tank(out.Idx_nominal))-273.15)
        fprintf(" - Average nominal injector pressure drop: %1.2f bar\n",1e-5*integral_mean(out.t(out.Idx_nominal),out.deltaP_inj(out.Idx_nominal)))
        fprintf(" - Remaining thermal protection thickness is: %1.2f mm\n",(rocket.cc.D_ext/2 - out.r_port(end))*1000)
        fprintf(" - Average throat diameter: %1.2f mm\n",1e3*integral_mean(out.t(out.Idx_nominal),out.r_throat(out.Idx_nominal))*2)
        fprintf(" - Average expansion ratio: %1.2f\n",integral_mean(out.t(out.Idx_nominal),out.eps(out.Idx_nominal)))
    end 

    % plotto i grafici
    if opts.plots
        plotGraphs(out,rocket,opts)
    end
end


function saveDynamicFile(out,opts,batch_fields,batch,values)
    % nome dinamico del file
    
    for i = 1:length(values)
        if length(batch{i}) == 1
            continue
        end
        if isnumeric(values{i})

            writeField(i) = string(batch_fields{i});
            if values{i} < 0.1
                exponent = floor(log10(values{i}));
                writeValues(i) = sprintf("%1.2fe%i",values{i}/10^exponent,exponent);
            else
                writeValues(i) = sprintf("%1.2f",values{i});
            end
        else
            writeField(i) = string(batch_fields{i});

            writeValues(i) = string(values{i});
        end
    end
    writeField = rmmissing(writeField);
    writeValues = rmmissing(writeValues);
    fullFileNameCommand = "fullFileName = opts.fileName";
    for i = 1:length(writeValues)
        fullFileNameCommand = fullFileNameCommand + "+'_' + writeField("+i+") + '=' + writeValues("+i+")";
    end
    eval(fullFileNameCommand+";")
    % creo la cartella se non esiste
    if ~exist(opts.folderName,"dir")
        mkdir(opts.folderName)
    end
    % costruisco il nome del file da salvare
    name = opts.folderName+"\"+fullFileName + ".mat";
    % salvo
    save(name,"out");
end

function plotGraphs(out,rocket,opts)
    
    linewidth = 1;

    % definisco qualche funzione utile
    plot_mean = @(m,color) yline(m,"--","mean: " + sprintf("%1.2f",m),...
        'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right',...
        "Color",color);
    set_axis = @(m,M) ylim([m,M]);

    % creo la prima figura e la splitto in una matrice 2x3 di sottofigure
    figure(1)
    m = 2;
    n = 3;
    
    % massa di ossidante rimasta
    subplot(m,n,1)
    plot(out.t,out.m_tank,"LineWidth",linewidth)
    hold on
    title('oxidizer mass vs time')
    grid on
    ylabel('m [kg]')
    xlabel('t [s]')
    set_axis(0,max(out.m_tank))
    
    % pressioni
    subplot(m,n,2)
    p1 = plot(out.t,out.p_tank/1e+5,"LineWidth",linewidth);
    hold on
    p2 = plot(out.t,out.p_cc/1e+5,"LineWidth",linewidth);
    title('pressures vs time')
    grid on
    ylabel('p [bar]')
    xlabel('t [s]')
    set_axis(0,max(out.p_tank)/1e+5)
    if opts.plotInfo
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.p_tank(out.Idx_nominal))/1e+5,p1.Color)
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.p_cc(out.Idx_nominal))/1e+5,p2.Color)
    end
    legend([p1,p2],["p_{tank}","p_{cc}"],'location','southwest')
    
    % portate di massa
    subplot(m,n,3)
    p1 = plot(out.t,out.mDot_ox,"LineWidth",linewidth);
    hold on
    p2 = plot(out.t,out.mDot_fuel,"LineWidth",linewidth);
    p3 = plot(out.t,out.mDot_ox + out.mDot_fuel,"LineWidth",linewidth);
    title('mass flow rates vs time')
    grid on
    ylabel('mDot [kg/s]')
    xlabel('t [s]')
    set_axis(0,max(out.mDot_ox + out.mDot_fuel))
    if opts.plotInfo
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.mDot_ox(out.Idx_nominal)),p1.Color)
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.mDot_fuel(out.Idx_nominal)),p2.Color)
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.mDot_tot(out.Idx_nominal)),p3.Color)
    end
    legend([p1,p2,p3], ["mDot_{ox}","mDot_f","mDot_{tot}"],'location','best')

    % spinta e condizione dell'ugello
    subplot(m,n,4)
    p1 = plot(out.t,out.S,"LineWidth",linewidth);
    hold on   
    title('thrust vs time')
    grid on
    ylabel('S [N]')
    xlabel('t [s]')
    set_axis(0,max(out.S))
    if opts.plotInfo
        p2 = plot(out.t,(out.cond+1)*max(out.S)/7,"LineWidth",linewidth);
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.S(out.Idx_nominal)),"k")             
        yline(max(out.S)/7*1,':','no flow',"Color",p2.Color)
        yline(max(out.S)/7*2,':','sub.',"Color",p2.Color)
        yline(max(out.S)/7*3,':','shock',"Color",p2.Color)
        yline(max(out.S)/7*4,':','o.e.',"Color",p2.Color)
        yline(max(out.S)/7*5,':','p.a.',"Color",p2.Color)
        yline(max(out.S)/7*6,':','u.e.',"Color",p2.Color)
        legend([p1,p2],["thrust","nozzle condition"],'location','southwest')
    end
    

    % Isp
    subplot(m,n,5)
    plot(out.t,out.Isp,"LineWidth",linewidth)
    hold on
    title('I_{sp} vs time')
    grid on
    ylabel('I_{sp} [s]')
    xlabel('t [s]')
    set_axis(0,max(out.Isp))
    if opts.plotInfo
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.Isp(out.Idx_nominal)),"k");
    end

    % O/F
    subplot(m,n,6)
    plot(out.t,out.O_F,"LineWidth",linewidth)
    hold on
    title('O/F vs time')
    grid on
    ylabel('O/F')
    xlabel('t [s]')
    set_axis(0,max(out.O_F))
    if opts.plotInfo
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.O_F(out.Idx_nominal)),"k");
    end

    % nuova figura
    figure(2)

    % area ratio
    subplot(m,n,1)
    plot(out.t,out.eps,"LineWidth",linewidth)
    hold on
    title('area ratio vs time')
    grid on
    ylabel('\epsilon')
    xlabel('t [s]')
    set_axis(0,max(out.eps))
    if opts.plotInfo
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.eps(out.Idx_nominal)),"k");
    end
    
    % erosion rate della gola
    subplot(m,n,2)
    plot(out.t,out.rDot_throat*1000,"LineWidth",linewidth)
    hold on
    title('throat erosion rate vs time')
    grid on
    ylabel('rDot_{throat} [mm/s]')
    xlabel('t [s]')
    set_axis(0,max(out.rDot_throat)*1000+1e-6)
    if opts.plotInfo
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.rDot_throat(out.Idx_nominal))*1000,"k");
    end

    % raggio interno del grano
    subplot(m,n,3)
    plot(out.t,out.r_port*1000,"LineWidth",linewidth)
    title('port radius vs time')
    hold on
    grid on
    ylabel('D [mm]')
    xlabel('t [s]')
    set_axis(0,rocket.cc.D_ext * 1.1*1000/2)
    if opts.plotInfo
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.r_port(out.Idx_nominal))*1000,"k")
        yline(rocket.cc.D_ext*1000/2,'k','case','LabelVerticalAlignment','bottom')
        yline(rocket.cc.D_port0*1000/2,'k--',"min: " + rocket.cc.D_port0*1000/2,...
            'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right')
        yline(out.r_port(end)*1000,'k--',"max: " + out.r_port(end)*1000,...
            'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right')
    end
    
    % erosion rate del grano
    subplot(m,n,4)
    plot(out.t,out.rDot_port*1000,"LineWidth",linewidth)
    hold on
    title('regression rate vs time')
    grid on
    ylabel('r [mm/s]')
    xlabel('t [s]')
    set_axis(0,max(out.rDot_port)*1000)
    if opts.plotInfo
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.rDot_port(out.Idx_nominal))*1000,"k")
    end

    % delta p attraverso l'iniettore
    deltaP = out.deltaP_inj;
    subplot(m,n,5)
    plot(out.t,deltaP/1e+5,"LineWidth",linewidth)
    hold on
    title('\DeltaP vs time')
    grid on
    ylabel('\DeltaP')
    xlabel('t [s]')
    set_axis(0,max(deltaP)/1e+5)
    if opts.plotInfo
        plot_mean(integral_mean(out.t(out.Idx_nominal),deltaP(out.Idx_nominal))/1e+5,"k");
    end
    
    % diametro di gola
    subplot(m,n,6)
    plot(out.t,out.r_throat * 2 * 1000,"LineWidth",linewidth)
    title('throat diameter vs time')
    hold on
    grid on
    ylabel('D_{throat} [mm]')
    xlabel('t [s]')
    set_axis(0,rocket.nozzle.D_exit/2*2 * 1.1 * 1000)
    if opts.plotInfo
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.r_throat(out.Idx_nominal))*2 * 1000,"k")
        yline(rocket.nozzle.D_exit*1000,'k','exit diameter','LabelVerticalAlignment','bottom')
    end

    
    figure(3)
    m = 1;
    n = 3;
    % flusso di massa
    subplot(m,n,1)
    plot(out.t,out.G,"LineWidth",linewidth)
    hold on
    title("G vs time")
    grid on
    xlabel('t [s]')
    ylabel("G [kg/sm^2]")
    if opts.plotInfo
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.G(out.Idx_nominal)),"k")
    end
    
    % temperatura nel serbatoio
    subplot(m,n,2)
    plot(out.t,out.T_tank-273.15,"LineWidth",linewidth)
    hold on
    title("Tank temperature vs time")
    grid on
    xlabel('t [s]')
    ylabel("T [°C]")
    if opts.plotInfo
        plot_mean(integral_mean(out.t(out.Idx_nominal),out.T_tank(out.Idx_nominal))-273.15,"k")
    end

    % titolo di vapore
    subplot(m,n,3)
    plot(out.t,out.x_tank,"LineWidth",linewidth)
    hold on
    title("Vapor quality vs time")
    grid on
    xlabel('t [s]')
    ylabel("x")
    drawnow

    figure(4)
    m = 1;
    n = 3;

    % perdite fluidiche
    subplot(m,n,1)
    plot(out.t,out.p_loss/1e5,"LineWidth",linewidth)
    hold on
    title("Fluidic losses vs time")
    grid on
    xlabel('t [s]')
    ylabel("p_loss [bar]")
    drawnow

    % Discharge coefficient Cd
    subplot(m,n,2)
    plot(out.t,out.c_d,"LineWidth",linewidth)
    hold on
    title("Cd vs time")
    grid on
    xlabel('t [s]')
    ylabel("Cd")
    drawnow

    % Velocità nell'iniettore
    subplot(m,n,3)
    plot(out.t,out.v,"LineWidth",linewidth)
    hold on
    title("Velocità nell'iniettore vs time")
    grid on
    xlabel('t [s]')
    ylabel("v [m/s]")
    drawnow

    figure(5)
    m = 1;
    n = 2;

    % massa di liquido all'interno del serbatoio
    subplot(m,n,1)
    plot(out.t,out.m_liq,"LineWidth",linewidth)
    hold on
    title("Liquid mass in the tank vs time")
    grid on
    xlabel('t [s]')
    ylabel("m_liq [kg]")
    drawnow

    % massa di vapore all'interno del serbatoio
    subplot(m,n,2)
    plot(out.t,out.m_vap,"LineWidth",linewidth)
    hold on
    title("Vapour mass in the tank vs time")
    grid on
    xlabel('t [s]')
    ylabel("m_vap [kg]")
    drawnow
 
    figure(6)
    m = 1;
    n = 2;

    % portata massica di ossidante in funzione del DeltaP a cavallo
    % dell'iniettore
    subplot(m,n,1)
    plot(out.deltaP_inj*1e-5,out.mDot_ox,"LineWidth",linewidth)
    hold on
    title("Portata massica di ossidante in funzione del deltaP")
    grid on
    xlabel('DeltaP [bar]')
    ylabel("mDotOX [kg/s]")
    drawnow

    figure(7)
    deltaP = out.deltaP_inj;
    deltaPperc=deltaP/(out.p_cc/100);
    plot(out.t,deltaPperc,"LineWidth",linewidth);
    hold on
    title("DeltaP percentuale nel tempo")
    grid on
    xlabel('t [s]')
    ylabel("DeltaP/p_c_c [%]")
    %char="Diametro di gola =" + rocket.nozzle.D_throat0*1000 + "mm";
    %legend([p],{char})
end

function rocket = buildRocket(rocket,env)
    % calcola tutte le variabili secondarie e le funzioni necessarie a far
    % funzionare la simulazione

    rocket.cc.V_inj_slot = rocket.cc.L_inj_slot * pi*(rocket.cc.D_inj_slot^2/4); % [m^3]
    rocket.cc.V_prechamber = rocket.cc.L_prechamber * pi*(rocket.cc.D_prechamber^2/4); % [m^3]
    rocket.cc.V_diaphragm = rocket.cc.L_diaphragm * pi*(rocket.cc.D_diaphragm^2/4); % [m^3]
    rocket.cc.V_grain = rocket.cc.L_grain * pi*(rocket.cc.D_port0^2/4); % [m^3]
    rocket.cc.V_postchamber = rocket.cc.L_postchamber * pi*(rocket.cc.D_postchamber^2/4); % [m^3]
    rocket.nozzle.V_inh = rocket.nozzle.L_inh * pi*(rocket.nozzle.D_inh^2/4); % [m^3]
    rocket.nozzle.V_convergent = rocket.nozzle.L_convergent/3 * pi * ...
        (rocket.nozzle.D_convergent^2/4 + rocket.nozzle.D_throat0^2/4 + ...
        rocket.nozzle.D_convergent * rocket.nozzle.D_throat0/4); % [m^3]
    rocket.nozzle.V_throat0 = rocket.nozzle.L_throat * pi*(rocket.nozzle.D_throat0^2/4); % [m^3]

    rocket.cc.V_tot0 = rocket.cc.V_inj_slot + rocket.cc.V_prechamber + ...
        rocket.cc.V_grain + rocket.cc.V_diaphragm + ...
        rocket.cc.V_postchamber + rocket.nozzle.V_inh + ...
        rocket.nozzle.V_convergent + rocket.nozzle.V_throat0; % [m^3]
    
    rocket.inj.A_holes = rocket.inj.n_holes * rocket.inj.D_holes^2/4 * pi; % [m^2] area totale dell'iniettore

    rocket.cc.R = 8314/rocket.cc.MM; % [kJ/kgK] constante del gas

    rocket.nozzle.A_throat0 = pi * rocket.nozzle.D_throat0^2/4;   % [m^2] area iniziale di gola
    rocket.nozzle.A_exit = rocket.nozzle.A_throat0 * rocket.nozzle.eps0; % [m^2] area di uscita
    rocket.nozzle.D_exit = sqrt(rocket.nozzle.A_exit/pi)*2;   % [m] diametro di uscita

    rhol_ox = NFP("N2O","rl_t",rocket.tank.T0);
    rhov_ox = NFP("N2O","rv_t",rocket.tank.T0);

    V_vap = rocket.tank.V * rocket.tank.f_min;
    m0_max = V_vap * rhov_ox + (rocket.tank.V - V_vap) * rhol_ox;

    % risolvi il serbatoio
    if isstring(rocket.tank.m0)
        if lower(rocket.tank.m0) == "max"
            % caso di massima massa inseribile
            m0 = m0_max;            
        else
            error("string not recognized")
        end
    else
        % caso classico di massa data
        if rocket.tank.m0 > m0_max
            warning("la frazione di vapore è inferiore a quella minima")
        end
        m0 = rocket.tank.m0;
    end

    [rocket.tank.m_vec,rocket.tank.p_vec,rocket.tank.rho_out_vec,...
                rocket.tank.T_vec,rocket.tank.x_vec] = ...
                tankFunctionNFP(rocket.tank.T0,m0,rocket.tank.V,rocket.tank.k);
    
    % taglio la lunghezza del vettore perchè 2000 elementi sono troppi, lo
    % rendo lungo 100 interpolando in mezzo
%     n = 200;
%     t_val = linspace(1,length(rocket.tank.m_vec),n);
%     rocket.tank.m_vec = interp1(1:length(rocket.tank.m_vec),rocket.tank.m_vec,t_val);
%     rocket.tank.p_vec = interp1(1:length(rocket.tank.p_vec),rocket.tank.p_vec,t_val);
%     rocket.tank.rho_out_vec = interp1(1:length(rocket.tank.rho_out_vec),rocket.tank.rho_out_vec,t_val);
%     rocket.tank.T_vec = interp1(1:length(rocket.tank.T_vec),rocket.tank.T_vec,t_val);
%     rocket.tank.x_vec = interp1(1:length(rocket.tank.x_vec),rocket.tank.x_vec,t_val);
    
    % costruisco le funzioni 
    rocket.tank.p_func = @(m) interp1(rocket.tank.m_vec,rocket.tank.p_vec,m,'pchip');
    rocket.tank.rho_out_func = @(m) interp1(rocket.tank.m_vec,rocket.tank.rho_out_vec,m,'pchip');
    rocket.tank.T_func = @(m) interp1(rocket.tank.m_vec,rocket.tank.T_vec,m,'pchip');
    rocket.tank.x_func = @(m) interp1(rocket.tank.m_vec,rocket.tank.x_vec,m,'pchip');

    rocket.inj.fluidicFunction = @(p_cc,p_tank,m_tank) fluidicFunction(p_cc,p_tank,m_tank,rocket);  
    rocket.cc.grainFunction = @(mDot_ox,r_f) grainFunction(mDot_ox,r_f,rocket);     
    rocket.nozzle.nozzleErosion = @(t,p_cc,A_t,mDot_tot,O_F) nozzleErosion(t,p_cc,A_t,mDot_tot,O_F,rocket,env);
    rocket.nozzle.nozzleFunction = @(p_cc,A_t,eps,gasProp) nozzleFunction(p_cc,A_t,eps,gasProp,rocket,env);
    rocket.cc.combustionFunction = @(O_F,T_ox,T_f) combustionFunction(O_F,T_ox,T_f,rocket);
end

function mean = integral_mean(x,y)
   mean = integral(@(xx) interp1(x,y,xx),min(x),max(x))/(max(x) - min(x));
end