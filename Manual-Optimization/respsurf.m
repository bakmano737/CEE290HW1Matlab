% --------- Written by Jasper A. Vrugt ---------
% ------ University of California Irvine ------- 
% --------- CEE - 290 : Models & Data ----------

clear all; close all;

% The observations:

% Define the measured head
h = [0.55 0.47 0.30 0.22 0.17 0.14];
% Define the time
t = [5.0 10.0 20.0 30.0 40.0 50.0];

% Define the distance from injection
d = 10;

% Define amount of injected water
Q = 50;


% Calculate the response surface as a function of S and T
nS = 250;
nT = 500;
S = linspace(0.002,0.6,nS);
T = linspace(0.02,1.1,nT);
for iS=1:nS
    for iT=1:nT
        % Run the slug model
        h_pred = slugmodel([S(iS),T(iT)],t,Q,d);

        % Calculate the error residuals
        residuals = (h - h_pred);

        % Calculate the objective function value
        resp_surface(iT,iS) = sum(residuals.^2);
    end
end

figure
imagesc(S,T,resp_surface,[0,0.5])
hold on
plot(0.106 +[-1,-1]*0.0048,[0,2],'--',...
    0.106 + [1,1]*0.0048,[0,2],'--','linewidth',1,'color',0.9*[1,1,1]);
plot([0,0.6],0.507+[1,1]*0.0352,'--',...
    [0,0.6],0.507+[-1,-1]*0.0352,'--','linewidth',1,'color',0.9*[1,1,1]);
set(gca,'ydir','normal')
colorbar
title('response surface for slug injection (SSR)')
xlabel('S')
ylabel('T')


set(gcf,'paperpositionmode','auto','inverthardcopy','off')
print('respsurf-sluginj.eps','-depsc2','-r300','-loose')
print('respsurf-sluginj.png','-dpng','-r300')