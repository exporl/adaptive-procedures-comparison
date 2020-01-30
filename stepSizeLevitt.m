function stepSize = stepSizeLevitt(procedure,nTrials,nReversals)
% Calculate a step size (for each simulation) according to:
%     Levitt, H. C. C. H. (1971). Transformed up?down methods in 
%     psychoacoustics. The Journal of the Acoustical society of America, 
%     49(2B), 467-477.
%
% INPUTS:
%     procedure             struct that defines the procedure parameters, 
%                           with required fields:
%      procedure.stepSizes      step sizes that are used throughout the
%                               procedure
%      procedure.stepChangeAfters
%                               number of trials/reverals before the step
%                               size is changed
%      procedure.stepChangeAfterType
%                               either 'trials' or 'reversals'
%     nTrials               number of passed trials (for each simulation)
%     nReversals            number of passed reversals (for each
%                           simulation)
%
% OUTPUTS: 
%     stepSize              the calculated step size
%
% Author: Benjamin Dieudonné, KU Leuven, Department of Neurosciences, ExpORL
% Correspondence: tom.francart@med.kuleuven.be

procedure.stepSizes = procedure.stepSizes(:).'; % make row-vector
procedure.stepChangeAfters = procedure.stepChangeAfters(:).'; % make row-vector
nSimulations = size(nReversals,1);

switch procedure.stepChangeAfterType
    case 'trials'
        nTrials = nTrials*ones(size(procedure.stepChangeAfters));
        stepSize = ones(nSimulations,1)*procedure.stepSizes(sum(nTrials >= procedure.stepChangeAfters,2));
    
    case 'reversals'
        nReversals = nReversals*ones(1,size(procedure.stepChangeAfters,2));
        procedure.stepChangeAfters = ones(nSimulations,1)*procedure.stepChangeAfters;
        procedure.stepSizes = ones(nSimulations,1)*procedure.stepSizes;

        stepSizeIndices = sum(nReversals >= procedure.stepChangeAfters,2).';
        stepSizeIndices = sub2ind(size(procedure.stepSizes), 1:size(procedure.stepSizes, 1), stepSizeIndices);
        stepSize = procedure.stepSizes(stepSizeIndices).';
    
    otherwise
        error('This type of p.stepChangeAfterType is not supported');
        
end

end
