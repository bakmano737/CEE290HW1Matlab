% --------- Written by Jasper A. Vrugt ---------
% ------ University of California Irvine ------- 
% --------- CEE - 290 : Models & Data ----------

clear
close all
clc



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % %  define measurements and constants
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Define the measured head
headObs = [0.55 0.47 0.30 0.22 0.17 0.14];
% Define the time
timeObs = [5.0 10.0 20.0 30.0 40.0 50.0];
% Define the distance from injection
d = 10; 
% Define amount of injected water
Q = 50;

% default values for the parameters:
S = 0.15;
T = 0.4;


% make a big figure
bigfigure(1,1,1,'figureNumber',1001)
% prepare some arrays for the visualization 
hu = repmat(NaN,[1,6]);
hl = repmat(NaN,[1,6]);

for iSubplot=1:6
    if iSubplot>1
        strS = ['(Press Enter to accept previous value of ',num2str(S),')'];
        strT = ['(Press Enter to accept previous value of ',num2str(T),')'];
    else
        strS = ['(Press Enter to accept default value of ',num2str(S),')'];
        strT = ['(Press Enter to accept default value of ',num2str(T),')'];
    end

    tmp = input(['Enter a value for ''S'' ',strS,'. >> ']);
    if ~isempty(tmp)
        S=tmp;
    end
    tmp = input(['Enter a value for ''T'' ',strT,'. >> ']);
    if ~isempty(tmp)
        T=tmp;
    end

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % % % %  make prediction
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % concatentate the parameters into a 1xN vector
    pars = [S,T];
    % define for which the forward model should give a prediction
    timeSim = timeObs;
    % run the forward model with the parameters etc.
    headSim = slugmodel(pars,timeSim,Q,d);


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % % % %  visualization
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % re-run the slugmodel function to create a smooth interpolation 
    timeSimSmooth = linspace(timeObs(1),timeObs(end),200);
    % run the forward model with the parameters etc.
    headSimSmooth = slugmodel(pars,timeSimSmooth,Q,d);

    % do this in a convoluted way to deal with octave shizzle
    UWH=[0.27,0.35];
    LWH=[0.27,0.08];
    L = [0.05,0.38,0.72,0.05,0.38,0.72];
    UB = [0.64,0.64,0.64,0.14,0.14,0.14];
    LB = [0.55,0.55,0.55,0.05,0.05,0.05];

    if isnan(hu(iSubplot))
        hu(iSubplot)=axes('position',[L(iSubplot),UB(iSubplot),UWH]);
    else
        axes(hu(iSubplot))
        cla
    end
    plot(timeObs,headObs,'om',...
         timeSim,headSim,'b.',...
         timeSimSmooth,headSimSmooth,'-b')
    legend('obs','sim','Location','NorthEast')
    set(gca,'ylim',[0,0.7],'xticklabel',[])
    text(30,0.6,['S = ',num2str(S),char(10),'T = ',num2str(T)],...
         'fontsize',9,'horizontalalignment','center')
    ylabel('head [m]')
    if isnan(hl(iSubplot))
        hl(iSubplot)=axes('position',[L(iSubplot),LB(iSubplot),LWH]);
    else
        axes(hl(iSubplot))
        cla
    end
    stem(timeObs,headSim-headObs,'-r.')
    set(gca,'ylim',[-1,1]*0.2)
    xlabel('time')
    ylabel('residual')

    disp(['Results visualized in subplot ',num2str(iSubplot),char(10)])

end