function [parset] = Run_Optimization ( parsIni , n_iter , opt_method );

% --------- Written by Jasper A. Vrugt ---------
% ------ University of California Irvine ------- 
% --------- CEE - 290 : Models & Data ----------

% Set default values if insufficient number of input arguments
if nargin == 0,
    % define the starting point:
    parsIni = [0.40,0.60]; n_iter = 10; opt_method = 1;
end;
if nargin == 1,
    n_iter = 10; opt_method = 1;
end;
if nargin == 2,
    opt_method = 1;
end;

% Define number of iterations
nParsetsLM = n_iter; nParsetsGN = n_iter;

% define whether to show the response surface
showRespSurf = true;
% define whether to print the starting point's Jacobian
printJacobian = true;

switch opt_method 
    case 1,
        % Gauss-Newton
        includeGauNew = true; includeLevMar = false;
    case 2
        % Levenberg-Marquardt
        includeLevMar = true; includeGauNew = false;
end;

% define the lambda factor to use in Levenberg-Marquardt (higher value 
% (e.g. ~100) results in more steepest-descent-like behavior, lower 
% value (e.g. ~0.001) is more Gauss-Newton)
lambdaFactor = 2.0;
% scale factor to decrease the step size in Gauss-Newton for improved convergence
alphaGauNew = 1.0;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %                                                                   % % 
% %            NO NEED TO CHANGE ANYTHING BELOW THIS POINT            % % 
% %                                                                   % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% Define the measured head
h = [0.55 0.47 0.30 0.22 0.17 0.14];
% Define the time
t = [5.0 10.0 20.0 30.0 40.0 50.0];
% Define the distance from injection
d = 10; 
% Define amount of injected water
Q = 50;


% Calculate the response surface as a function of S and T
nS = 100; 
nT = 200;
S = linspace(0.002,0.6,nS);
T = linspace(0.02,1.1,nT);
if showRespSurf
    for iS=1:nS
        for iT=1:nT
            % Run the slug model
            h_pred = slugmodel([S(iS),T(iT)],t,Q,d);

            % Calculate the error residuals
            residuals = (h - h_pred);

            % Calculate the objective function value
            OF(iT,iS) = sum(residuals.^2);
        end
    end
end

% prepare a figure
bigfigure
clf
ax121=subplot(1,2,1);
set(gca,'Fontsize',14)
ax122=subplot(1,2,2);
set(gca,'Fontsize',14)

axes(ax121)
tSmooth = linspace(t(1),t(end),200);
hSimIni = slugmodel(parsIni,t,Q,d);
hSimIniSmooth = slugmodel(parsIni,tSmooth,Q,d);
handles121 = plot(t,h,'om','markerfacecolor','m');
legendEntries121 = {'obs'};
hold on
tmphandles = plot(t,hSimIni,'ko',...
                  tSmooth,hSimIniSmooth,'--k','linewidth',2);
handles121 = [handles121;tmphandles(2)];
legendEntries121 = [legendEntries121;'sim-ini'];
clear tmphandles

if showRespSurf
    axes(ax122)
    imagesc(S,T,log10(OF))
    colormap(flipud(gray))
    hold on
    contour(S,T,log10(OF),[-4:1:2],'-k','showText','on')
    g=colorbar;
    set(get(g,'title'),'string','log_{10}(OF)')
end
handles122=[];
legendEntries={};


if includeGauNew
    % Gauss-Newton method
    pars = parsIni;
    % Iteratively improve the initial starting parameter values
    for iParsetGN = 1:nParsetsGN

        % Run the model with the current parameter values
        h_pred = slugmodel(pars,t,Q,d); 
        % Calculate the residuals
        residuals = (h - h_pred)';
        % Store the objective function value
        parsetsGN(iParsetGN,1:3) = [pars,sum(residuals.^2)];
        % Calculate Jacobian
        J = jacSlugModel(pars,t,Q,d);
        % Calculate the update
        update = alphaGauNew * - ( inv(J'*J) ) * J' * residuals;
        % Update parameter values
        pars = pars + update';

    end

	% run the forward model with the parameters etc (for visualization purposes only).
    hSimEnd = slugmodel(parsetsGN(nParsetsGN,:),t,Q,d);
    hSimEndSmooth = slugmodel(parsetsGN(nParsetsGN,:),tSmooth,Q,d);
    axes(ax121)
    hold on
    tmphandles = plot(t,hSimEnd,'bs',...
                      tSmooth,hSimEndSmooth,'-b','linewidth',2);
    handles121 = [handles121;tmphandles(2)];
    legendEntries121 = [legendEntries121;'GN-sim-final'];
    clear tmphandles
 
    axes(ax122)
    hold on
    plot(parsetsGN(:,1),parsetsGN(:,2),'-bs','linewidth',2)
    % Return parset
    parset = parsetsGN;
end

if includeLevMar
	% Levenberg-Marquardt 
    lambdaFactorMin = 1e-12;
    lambdaFactorMax = 1e16;
    
	nPars = numel(parsIni);
    % Run the model with the initial parameter values
    h_pred = slugmodel(parsIni,t,Q,d);
    % Calculate the residuals
    residuals = (h - h_pred)';
    parsOld = parsIni;
    scoreOld = sum(residuals.^2);
    % start the record of points in parameter space:
    parsetsLM(1,1:nPars+1) = [parsOld,scoreOld];
    % keep a record of rejected points for visualization:
    parsetsLMRej = repmat(NaN,[nParsetsLM,nPars+1]);
    % Iteratively improve the initial starting parameter values
    for iParsetLM=2:nParsetsLM 
            % calculate the Jacobian at point parsOld
            J = jacSlugModel(parsOld,t,Q,d);
            % calculate the Levenberg Marquardt term
            levMarTerm = lambdaFactor * eye(nPars);
            % calculate the size of the update
            update = - inv(J'*J+levMarTerm)*J' * residuals;
            % Calculate the new point in the parameter space
            parsNew = parsOld + update';
            % Run the model with the new parameter values
            h_pred = slugmodel(parsNew,t,Q,d);
            % Calculate the residuals
            residuals = (h - h_pred)';
            % Calculate the new objective score
            scoreNew = sum(residuals.^2);
            if scoreNew<scoreOld
                % make lambdaFactor smaller:
                lambdaFactor = lambdaFactor/2;
                % accept the parsNew point in the list of points
                parsetsLM(iParsetLM,1:nPars+1) = [parsNew,scoreNew];
                % move to the new point:
                parsOld = parsNew;
                scoreOld = scoreNew; 
            else
                % make lambdaFactor larger:
                lambdaFactor = lambdaFactor*2;
                % accept the old point in the list of points
                parsetsLM(iParsetLM,1:nPars+1) = [parsOld,scoreOld];
                parsetsLMRej(iParsetLM,1:nPars+1) = [parsNew,scoreNew];
            end
            if lambdaFactor > lambdaFactorMax
                % make sure lambdaFactor doesn't grow larger
                % than 'lambdaFactorMax'
                lambdaFactor = lambdaFactorMax;
            end
            if lambdaFactor < lambdaFactorMin
                % make sure lambdaFactor doesn't become smaller
                % than 'lambdaFactorMin'
                lambdaFactor = lambdaFactorMin;
            end
    end

	% run the forward model with the parameters etc (for visualization purposes only).
	hSimEnd = slugmodel(parsetsLM(nParsetsLM,:),t,Q,d);
	hSimEndSmooth = slugmodel(parsetsLM(nParsetsLM,:),tSmooth,Q,d);
    axes(ax121)
    hold on
    tmphandles = plot(t,hSimEnd,'co',...
                      tSmooth,hSimEndSmooth,'-c','linewidth',2);
	handles121 = [handles121;tmphandles(2)];
    legendEntries121 = [legendEntries121;'LM-sim-final'];
    clear tmphandles
	axes(ax122)
    hold on
	plot(parsetsLM(:,1),parsetsLM(:,2),'-co','linewidth',2)
    % also include the rejected points:
    % plot(parsetsLMRej(:,1),parsetsLMRej(:,2),'co');
    % Return parset
    parset = parsetsLM;   
end

axes(ax121)
set(gca,'xlim',[0,50],'ylim',[0,0.7],'xtick',t)
xlabel('time')
ylabel('head [m]')
legend(handles121,legendEntries121)

axes(ax122)
axis image
set(gca,'xlim',[S(1),S(nS)],'ylim',[T(1),T(nT)],'ydir','normal')
xlabel('S')
ylabel('T')
if showRespSurf
    title(['response surface with contour',char(10),'lines of SSR in logspace(-4,2,7)'])
else
    title('[S,T] parameter space')
end
if printJacobian
	disp(sprintf('The Jacobian at point [S,T]=[%.3f,%.3f] is:',parsIni))
	J = jacSlugModel(parsIni,t,Q,d);
	disp(sprintf('%14.5g %14.5g \n',J'))
end
if includeGauNew
    disp('Gauss-Newton sampled these parameter sets:')
    disp('iteration              S              T             OF')
    disp(sprintf('%9d %14.5g %14.5g %14.5g \n',[(1:nParsetsGN)',parsetsGN]'))
end
if includeLevMar
    disp('Levenberg-Marquardt sampled these parameter sets:')
    disp('iteration              S              T             OF')
    disp(sprintf('%9d %14.5g %14.5g %14.5g \n',[(1:nParsetsLM)',parsetsLM]'))
end
% set(gcf,'paperpositionmode','auto','inverthardcopy','off')
% print(gcf,'gauss-newton-slugmodel.eps','-loose','-depsc','-r300')
