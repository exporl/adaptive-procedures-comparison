function convergence = calculateConvergence(procedure,values,reversals,nStaircasePlots)
% Calculate the convergence points for simulations of an adaptive
% procedure.
%
% INPUTS:
%     procedure             struct that defines the procedure parameters, 
%                           with required fields:
%      procedure.nLast      number of last trials/reversals to take into 
%                           account for the calculation of the convergence
%                           point
%      procedure.nLastType  either 'trials' or 'reversals'
%     values                matrix with a staircase for each simulation of 
%                           the respective procedure (the amount of rows 
%                           corresponds to the amount of simulations that 
%                           were done)
%     reversals             boolean matrix that specifies for each value of
%                           the starcaises whether a reversal occured.
%     nStaircasePlots       (optional): number of desired staircase plots 
%                           to give an idea of the adaptive track
%
% OUTPUTS: 
%     convergence           array with all estimated convergence points for
%                           each simulation (the length corresponds to the
%                           amount of simulations)
%
% Author: Benjamin Dieudonné, KU Leuven, Department of Neurosciences, ExpORL
% Correspondence: tom.francart@med.kuleuven.be

if nargin<4
    nStaircasePlots = 0;
end

nSimulations = size(values,1);
convergence = nan(nSimulations,1);
experimentsDropped = 0;

switch procedure.nLastType
    case 'trials'
        for i=1:nSimulations
            rowValues = values(i,~isnan(values(i,:)));
            convergence(i) = mean(rowValues(end-procedure.nLast+1:end));
        end
        
        if nStaircasePlots>0 % plot some staircases
            plotIDs = randi([1 nSimulations],nStaircasePlots);
            trialIndices = 1:size(values,2);

            hold on;
            for iPlotID = 1:length(plotIDs)
                plotID = plotIDs(iPlotID);
                h = plot(trialIndices,values(plotID,:),'LineWidth',1.5);
                plot(trialIndices(end-procedure.nLast+1:end),values(plotID,end-procedure.nLast+1:end),'.','color',get(h,'color'),'MarkerSize',10);
            end
            hold off;
            xlabel('Trial index');
        end
        
    case 'reversals'
        for i=1:nSimulations
            rowReversalValues = values(i,reversals(i,:));
            rowReversalValues = rowReversalValues(~isnan(rowReversalValues));
            if length(rowReversalValues)-procedure.nLast+1<1
                warning('Experiment %i of %i is dropped because it does not have enough reversals.',i,nSimulations);
                convergence(i) = nan;
                experimentsDropped = 1;
            else
                convergence(i) = mean(rowReversalValues(end-procedure.nLast+1:end));
            end
        end
        
        if nStaircasePlots>0 % plot some staircases
            plotIDs = randi([1 nSimulations],nStaircasePlots);
            trialIndices = 1:size(values,2);

            hold on;
            for iPlotID = 1:length(plotIDs)
                plotID = plotIDs(iPlotID);
                reversalIndices = trialIndices(reversals(plotID,:));
                rowReversalValues = values(plotID,reversals(plotID,:));

                h = plot(trialIndices,values(plotID,:),'LineWidth',1.5);
                plot(reversalIndices(end-procedure.nLast+1:end),rowReversalValues(end-procedure.nLast+1:end),'.','color',get(h,'color'),'MarkerSize',10);
            end
            hold off;
            xlabel('Trial index');
        end
        
    otherwise
        error('This type of procedure.nLastType is not supported');
end

convergence = convergence(~isnan(convergence));

if experimentsDropped
    warning('%i of %i (%0.1f%%) experiments dropped because they did not have enough reversals.',nSimulations-length(convergence),nSimulations,100*(nSimulations-length(convergence))/nSimulations);
end

end