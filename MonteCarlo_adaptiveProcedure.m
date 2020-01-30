function [values,reversals] = MonteCarlo_adaptiveProcedure(procedure,psychCurve,nSimulations)
% Simulate an adaptive procedure for an imposed psychometric curve for a
% large amount of times.
%
% INPUTS:
%     procedure                 struct that defines the adaptive procedure,
%                               with required fields:
%      procedure.startValue     start value of the stimulus presentation
%      procedure.nUp            amount of trials with false response before
%                               the stimulus value increases (i.e., before
%                               the task becomes easier)
%      procedure.nDown          amount of trials with correct response 
%                               before the stimulus value increases (i.e., 
%                               before the task becomes easier)
%      procedure.stopAfter      number of trials/reversals before the
%                               a run ends
%      procedure.stopAfterType	either 'trials' or 'reversals'
%      procedure.stepSizes      step sizes that are used throughout the
%                               procedure
%      procedure.stepChangeAfters
%                               number of trials/reverals before the step
%                               size is changed
%      procedure.stepChangeAfterType
%                               either 'trials' or 'reversals'
%      procedure.stepSize       function that calculates the step size for
%                               each trial
% 
%     psychCurve                struct that defines the imposed
%                               psychometric curve, with required fields:      
%      psychCurve.function      function that gives the estimated
%                               performance (i.e., the chance of a correct
%                               response) as a function of stimulus value
%                               (mostly a sigmoid)
% 
%     nSimulations              amount of simulations
%
% OUTPUTS: 
%     values            matrix with a staircase for each simulation of the 
%                       respective procedure (the amount of rows 
%                       corresponds to the amount of simulations that were 
%                       done)
%     reversals         boolean matrix that specifies for each value of the
%                       starcaises whether a reversal occured.
%
% Author: Benjamin Dieudonné, KU Leuven, Department of Neurosciences, ExpORL
% Correspondence: tom.francart@med.kuleuven.be




switch procedure.stopAfterType
    case 'trials'
        procedure.maxTrials = procedure.stopAfter;
        
    case 'reversals'
        procedure.maxTrials = round(procedure.stopAfter*5*(procedure.nUp+procedure.nDown)/2); % "arbitrary" choice
        
    otherwise
        error('This type of procedure.stopAfterType is not supported');
        
end

values = nan(nSimulations, procedure.maxTrials + 1);
reversals = zeros(nSimulations, procedure.maxTrials + 1);
values(:,1) = procedure.startValue;

response = @(current) 100*rand(size(current)) < psychCurve.function(current); % [ 0 / 1 ]

nCorrect = zeros(nSimulations, 1); % number of correct responses in a row after change of value
nWrong = zeros(nSimulations, 1); % number of wrong responses in a row after change of value

nDown = procedure.nDown*ones(nSimulations, 1);
nUp = procedure.nUp*ones(nSimulations, 1);

isUp = zeros(nSimulations, 1);
isDown = zeros(nSimulations, 1);

stopIndices = nan(nSimulations, 1);

% do maxTrials number of trials
for trialindex = 1:procedure.maxTrials
    currentValue = values(:,trialindex);
    currentResponse = response(currentValue);
    
    nCorrect = ( nCorrect + currentResponse ) .* currentResponse; % +1 if current=correct, 0 if current=wrong
    nWrong = ( nWrong + 1 ) .* (1-currentResponse);  % +1 if current=wrong, 0 if current=correct
    
    step = (-1)*( nCorrect == nDown ) + (+1)*( nWrong == nUp ); % make next condition harder if nCorrect == nDown, easier if nWrong == nUp
    
    nCorrect = nCorrect.*(1 - ( nCorrect == nDown )); % reset after change of value
    nWrong = nWrong.*(1 - ( nWrong == nUp )); % reset after change of value
    
    nTrials = trialindex-1; % number of finished trials up to now (finished: when the next trial is completely defined)
    isUpPrev = isUp;
    isDownPrev = isDown;
    isUp = (step==0).*isUpPrev + (step > 0); % check if it's going up (can also be plateau)
    isDown = (step==0).*isDownPrev + (step < 0);
    reversals(:, trialindex) = (isUpPrev~=isUp).*(isDownPrev~=isDown);
    nReversals = sum(reversals,2); % number of reversals up to now
    stopIndices(nReversals == procedure.stopAfter & isnan(stopIndices)) = trialindex; % note: if stopAfterType == 'trials', then ALWAYS nReversals < stopAfter
    
    stepSize = procedure.stepSize(procedure,nTrials,nReversals);
    nextValue = currentValue + step.*stepSize;
    nextValueBounded = min(max(nextValue,procedure.minValue),procedure.maxValue);
    if ( sum( nextValue ~= nextValueBounded ) )
        percentageSaturated = 100*sum( nextValue ~= nextValueBounded )/nSimulations;
        nextValue = nextValueBounded;
%      warning('Procedure saturated (trial %i, %0.5f%% of experiments)', trialindex, percentageSaturated);
    end
    
    values(:, trialindex+1) = nextValue;
end

% truncate according to stopping rule
for iSimulation = 1:nSimulations
    stopIndex = stopIndices(iSimulation);
    if ~isnan(stopIndex)
        values(iSimulation, stopIndex + 1:end) = nan; % note: if stopAfterType == 'trials', then ALWAYS nReversals <= stopAfter
        reversals(iSimulation, stopIndex + 1:end) = 0;
    end
end

while strcmp(procedure.stopAfterType,'reversals') && sum(nReversals<procedure.stopAfter)
    nUndone = sum(nReversals<procedure.stopAfter);
    percentageUndone = 100*nUndone/nSimulations;
    warning('Some simulations did not end (%i of %i, %0.1f%% of simulations). Removed those and replaced with random others', nUndone,nSimulations, percentageUndone);
    selectRandom = randi([1 nSimulations],nUndone,1);
    values(nReversals<procedure.stopAfter,:) = values(selectRandom,:);
    reversals(nReversals<procedure.stopAfter,:) = reversals(selectRandom,:);
    nReversals(nReversals<procedure.stopAfter) = nReversals(selectRandom);
end
reversals = logical(reversals);

end