% This file compares different two-down one-up adaptive procedures in terms
% of power and duration, by simulating different procedures for
% different imposed psychometric functions. 
%   - The psychometric functions were estimated with a constant procedure 
%     for a discrimination task of interaural time differences, with 
%     different types of stimuli. It is then investigated how well the 
%     different procedures can estimate the differences in discrimination 
%     performance for the different types of stimuli. 
%   - The procedures are defined by the following parameters:
%       - step size
%       - when and how to change the step size (after a number of trials or
%       reversals)
%       - the start value of each run
%       - when a run stops (after a number of trials or reversals)
%       - how to calculate the convergence point of a run
%       - the number of runs per condition
%       - the number of participants
%
% Author: Benjamin Dieudonné, KU Leuven, Department of Neurosciences, ExpORL
% Correspondence: tom.francart@med.kuleuven.be

%% 
clear all
close all
clc

%% Define parameters for desired outcomes
alpha = 0.05; % significance level to calculate power
nStaircasePlots = 5; % plot some (randomly picked) staircase plots for each procedure

%% Define all parameters for the simulations
% Functions to switch between logarithmic and linear scale
ITDlin = @(ITDdB) 10*10.^((ITDdB)/10); % re 10 µs
ITDdB = @(ITDlin) 10*log10(ITDlin/10); % re 10 µs

p = struct;

% MONTE CARLO PARAMETERS
p.nSimulations = 10000;

% PSYCHOMETRIC CURVES PARAMETERS
nPsychCurves = 2;
sigmoid = @(x, p) p.chance + (100 - p.chance) ./ (1 + exp(-(x - p.mu) ./ p.sigma)); % [% ~ dB]
p.psychCurves{1} = struct;
p.psychCurves{1}.chance = 50; % [%]
p.psychCurves{1}.sigma = 3.8; % [dB]
p.psychCurves{2} = p.psychCurves{1};

p.psychCurves{1}.mu = ITDdB(181/2); % [dB] % VALUE FOR CONDITION WITHOUT QUANTIZATION -- DIVIDED BY TWO TO OBTAIN PRESENTED ITD
p.psychCurves{1}.function = @(x) sigmoid(x,p.psychCurves{1});
p.psychCurves{2}.mu = ITDdB(262/2); % [dB] % VALUE FOR CONDITION WITH QUANTIZATION -- DIVIDED BY TWO TO OBTAIN PRESENTED ITD
p.psychCurves{2}.function = @(x) sigmoid(x,p.psychCurves{2});

% PROCEDURE PARAMETERS
nProcedures = 3;
ps = cell(nProcedures,1);

% GENERAL PROCEDURE PARAMETERS
p.procedure = struct;
p.procedure.startValue = ITDdB(2000); % [dB]
p.procedure.maxValue = ITDdB(2000);
p.procedure.minValue = ITDdB(10);

p.procedure.nUp = 1;
p.procedure.nDown = 2;

p.procedure.stepSize = @(p,nTrials,nReversals) stepSizeLevitt(p,nTrials,nReversals);

% PROCEDURE 1 PARAMETERS (our experiment)
ps{1} = p;
ps{1}.title = 'Our experiment';
ps{1}.nSubjects = 16;
ps{1}.nRuns = 2;
ps{1}.nSimulations = p.nSimulations*ps{1}.nRuns*ps{1}.nSubjects;

ps{1}.procedure.stopAfter = 25;
ps{1}.procedure.stopAfterType = 'trials'; % trials / reversals
    
ps{1}.procedure.stepSizes = [3,2,1];
ps{1}.procedure.stepChangeAfters = [0,8,20];
ps{1}.procedure.stepChangeAfterType = 'trials'; % trials / reversals

ps{1}.procedure.nLast = 12;
ps{1}.procedure.nLastType = 'trials'; % trials / reversals

% PROCEDURE 2 PARAMETERS (Bernstein, 2002)
ps{2} = p;
ps{2}.title = 'Bernstein et al., 2002';
ps{2}.nSubjects = 4;
ps{2}.nRuns = 6;
ps{2}.nSimulations = p.nSimulations*ps{2}.nRuns*ps{2}.nSubjects;

ps{2}.procedure.stopAfter = 12;
ps{2}.procedure.stopAfterType = 'reversals'; % trials / reversals
    
ps{2}.procedure.stepSizes = [2,0.5];
ps{2}.procedure.stepChangeAfters = [0,2];
ps{2}.procedure.stepChangeAfterType = 'reversals'; % trials / reversals

ps{2}.procedure.nLast = 10;
ps{2}.procedure.nLastType = 'reversals'; % trials / reversals

% PROCEDURE 3 PARAMETERS (Ehlers, 2016)
ps{3} = p;
ps{3}.title = 'Ehlers et al., 2016';
ps{3}.nSubjects = 11;
ps{3}.nRuns = 1;
ps{3}.nSimulations = p.nSimulations*ps{3}.nRuns*ps{3}.nSubjects;

ps{3}.procedure.stopAfter = 10;
ps{3}.procedure.stopAfterType = 'reversals'; % trials / reversals
    
ps{3}.procedure.stepSizes = [10*log10(3),10*log10(2),10*log10(sqrt(2))];
ps{3}.procedure.stepChangeAfters = [0,2,4];
ps{3}.procedure.stepChangeAfterType = 'reversals'; % trials / reversals

ps{3}.procedure.nLast = 6;
ps{3}.procedure.nLastType = 'reversals'; % trials / reversals

%% Do simulation and calculate outcome variables
fprintf('Starting simulations... \n');
if nStaircasePlots>0
    figure; 
end
for iProcedure = 1:nProcedures 
    fprintf('\tProcedure %d of %d (%s)... \n',iProcedure,nProcedures,ps{iProcedure}.title);   
    for iPsychCurve = 1:nPsychCurves
        fprintf('\t\tPsychometric curve %d of %d... \n',iPsychCurve,nPsychCurves);
        
        [values,reversals] = MonteCarlo_adaptiveProcedure(ps{iProcedure}.procedure,ps{iProcedure}.psychCurves{iPsychCurve},ps{iProcedure}.nSimulations);
        
        % Calculate convergence points (per run)        
        if nStaircasePlots>0 && iPsychCurve==1
            subplot(nProcedures,1,iProcedure);
            convergence = calculateConvergence(ps{iProcedure}.procedure,values,reversals,nStaircasePlots);
            axis([0 80 ITDdB(10) ITDdB(2000)])
            yticks(ITDdB([10 50 100 500 1000 2000]))
            yticklabels([10 50 100 500 1000 2000])
            title(ps{iProcedure}.title);
        else
            convergence = calculateConvergence(ps{iProcedure}.procedure,values,reversals);
        end
        
        % Calculate statistics: durations
        duration = sum(~isnan(values),2)-1;
        ps{iProcedure}.psychCurves{iPsychCurve}.durationPerRun = duration;
        ps{iProcedure}.psychCurves{iPsychCurve}.durationPerSubject = sum(reshape(duration,ps{iProcedure}.nRuns,p.nSimulations*ps{iProcedure}.nSubjects),1).'; % sum over runs per subject => estimation per SUBJECT
        ps{iProcedure}.psychCurves{iPsychCurve}.totalDuration = sum(reshape(duration,ps{iProcedure}.nRuns*ps{iProcedure}.nSubjects,p.nSimulations),1).'; % estimation per COMPLETE EXPERIMENT = SAMPLE DISTRIBUTION
        
        ps{iProcedure}.psychCurves{iPsychCurve}.durationPerRunMean = mean(ps{iProcedure}.psychCurves{iPsychCurve}.durationPerRun);
        ps{iProcedure}.psychCurves{iPsychCurve}.durationPerSubjectMean = mean(ps{iProcedure}.psychCurves{iPsychCurve}.durationPerSubject);
        ps{iProcedure}.psychCurves{iPsychCurve}.totalDurationMean = mean(ps{iProcedure}.psychCurves{iPsychCurve}.totalDuration);

        % Calculate statistics: convergence points (JND/percentageCorrect)
        ps{iProcedure}.psychCurves{iPsychCurve}.convergence = mean(reshape(convergence,ps{iProcedure}.nRuns*ps{iProcedure}.nSubjects,p.nSimulations)).'; % average over subjects and runs per subject => estimation per COMPLETE EXPERIMENT = SAMPLE DISTRIBUTION
        ps{iProcedure}.psychCurves{iPsychCurve}.convergenceMean = mean(ps{iProcedure}.psychCurves{iPsychCurve}.convergence); % SAMPLE DISTRIBUTION MEAN
        ps{iProcedure}.psychCurves{iPsychCurve}.convergenceStd = std(ps{iProcedure}.psychCurves{iPsychCurve}.convergence); % SAMPLE DISTRIBUTION STD

        ps{iProcedure}.psychCurves{iPsychCurve}.JNDMean = 2*ITDlin(ps{iProcedure}.psychCurves{iPsychCurve}.convergenceMean); % *2 FOR JNDs !
        ps{iProcedure}.psychCurves{iPsychCurve}.JNDMin = 2*ITDlin(ps{iProcedure}.psychCurves{iPsychCurve}.convergenceMean-ps{iProcedure}.psychCurves{iPsychCurve}.convergenceStd); % *2 FOR JNDs !
        ps{iProcedure}.psychCurves{iPsychCurve}.JNDMax = 2*ITDlin(ps{iProcedure}.psychCurves{iPsychCurve}.convergenceMean+ps{iProcedure}.psychCurves{iPsychCurve}.convergenceStd); % *2 FOR JNDs !

        ps{iProcedure}.psychCurves{iPsychCurve}.percentageCorrectMean = p.psychCurves{iPsychCurve}.function(ps{iProcedure}.psychCurves{iPsychCurve}.convergenceMean);    
        ps{iProcedure}.psychCurves{iPsychCurve}.percentageCorrectMin = p.psychCurves{iPsychCurve}.function(ps{iProcedure}.psychCurves{iPsychCurve}.convergenceMean-...
            ps{iProcedure}.psychCurves{iPsychCurve}.convergenceStd);
        ps{iProcedure}.psychCurves{iPsychCurve}.percentageCorrectMax = p.psychCurves{iPsychCurve}.function(ps{iProcedure}.psychCurves{iPsychCurve}.convergenceMean+...
            ps{iProcedure}.psychCurves{iPsychCurve}.convergenceStd);
    end
    
    % Calculate statistics: power and variance (by comparing the two psychCurves)
    ps{iProcedure}.significanceLevel = quantile(ps{iProcedure}.psychCurves{1}.convergence,1-alpha);
    ps{iProcedure}.beta = sum(ps{iProcedure}.psychCurves{2}.convergence<=ps{iProcedure}.significanceLevel)/p.nSimulations;
    ps{iProcedure}.powerPercent = 100*(1-ps{iProcedure}.beta);
    
    ps{iProcedure}.variance = var(ps{iProcedure}.psychCurves{1}.convergence);
end
fprintf('Simulations done! \n');


%% Print estimations
fprintf('Results: \n');
for iProcedure = 1:nProcedures
    fprintf('\t%s:\n',ps{iProcedure}.title);
    fprintf('\t\tMean performance: %0.1f us [%0.1f us, %0.1f us] ~ %0.1f %% [%0.1f %%, %0.1f %%] \n',...
        ps{iProcedure}.psychCurves{1}.JNDMean,ps{iProcedure}.psychCurves{1}.JNDMin,ps{iProcedure}.psychCurves{1}.JNDMax,...
        ps{iProcedure}.psychCurves{1}.percentageCorrectMean,ps{iProcedure}.psychCurves{1}.percentageCorrectMin,ps{iProcedure}.psychCurves{1}.percentageCorrectMax);
    fprintf('\t\tMean duration per run: %0.1f trials\n',...
        ps{iProcedure}.psychCurves{1}.durationPerRunMean);
    fprintf('\t\tMean duration per condition per subject: %0.1f trials \n',...
        ps{iProcedure}.psychCurves{1}.durationPerSubjectMean);
    fprintf('\t\tMean total duration per condition (summed over subjects): %0.1f trials \n',...
        ps{iProcedure}.psychCurves{1}.totalDurationMean);
    fprintf('\t\tPower: %.2f%%\n',ps{iProcedure}.powerPercent);
    fprintf('\t\tVariance: %.2f dB\n',ps{iProcedure}.variance);
end

%% Plot duration boxplots
% Labels etc
labels = cell(nProcedures,1);
for iProcedure = 1:nProcedures
    labels{iProcedure} = ps{iProcedure}.title;
end
width = 0.8/nPsychCurves;
legendEntries = strseq('Psychometric curve ',1:nPsychCurves).';

% Average duration per run
figure;
subplot(3,1,1)
hold on;
colors = get(gca,'colororder');
colors(1,:) = [0 0 1];
for iProcedure = 1:nProcedures
    iPsychCurve = 1;
    
    dataWithDummies = nan(length(ps{iProcedure}.psychCurves{iPsychCurve}.durationPerRun),nProcedures);
    dataWithDummies(:,iProcedure) = ps{iProcedure}.psychCurves{iPsychCurve}.durationPerRun;
    boxplot(dataWithDummies,'labels',labels,'position',(1:nProcedures)-1+(iPsychCurve-1)*width+0.5, 'widths',width,'Color', colors(1+iProcedure,:),'Symbol','.');
    plot(NaN,1,'color', colors(iPsychCurve,:)); % dummy plot for legend
end
for iProcedure = 1:nProcedures
    dataDummies = nan(2,nProcedures);
    boxplot(dataDummies,'labels',labels,'position',(1:nProcedures)-1+0.5, 'widths',width);
end
hold off;
xlim([0 nProcedures])
ylim('auto');
title('Duration per run');
ylabel('Duration [number of trials]');
set(gca,'FontSize',18);

% Average duration per subject
subplot(3,1,2)
hold on;
for iProcedure = 1:nProcedures
    iPsychCurve = 1;
    
    dataWithDummies = nan(length(ps{iProcedure}.psychCurves{iPsychCurve}.durationPerSubject),nProcedures);
    dataWithDummies(:,iProcedure) = ps{iProcedure}.psychCurves{iPsychCurve}.durationPerSubject;
    boxplot(dataWithDummies,'labels',labels,'position',(1:nProcedures)-1+(iPsychCurve-1)*width+0.5, 'widths',width,'Color', colors(1+iProcedure,:),'Symbol','.');
    plot(NaN,1,'color', colors(iPsychCurve,:)); % dummy plot for legend
end
for iProcedure = 1:nProcedures
    dataDummies = nan(2,nProcedures);
    boxplot(dataDummies,'labels',labels,'position',(1:nProcedures)-1+0.5, 'widths',width);
end
hold off;
xlim([0 nProcedures])
ylim('auto');
title('Duration per subject');
ylabel('Duration [number of trials]');
set(gca,'FontSize',18);

% Average duration per condition
subplot(3,1,3)
hold on;
for iProcedure = 1:nProcedures
    iPsychCurve = 1;
    
    dataWithDummies = nan(length(ps{iProcedure}.psychCurves{iPsychCurve}.totalDuration),nProcedures);
    dataWithDummies(:,iProcedure) = ps{iProcedure}.psychCurves{iPsychCurve}.totalDuration;
    boxplot(dataWithDummies,'labels',labels,'position',(1:nProcedures)-1+(iPsychCurve-1)*width+0.5, 'widths',width,'Color', colors(1+iProcedure,:),'Symbol','.');
    plot(NaN,1,'color', colors(iPsychCurve,:)); % dummy plot for legend
end
for iProcedure = 1:nProcedures
    dataDummies = nan(2,nProcedures);
    boxplot(dataDummies,'labels',labels,'position',(1:nProcedures)-1+0.5, 'widths',width);
end
hold off;
xlim([0 nProcedures])
ylim('auto');
title('Total duration per condition');
ylabel('Duration [number of trials]');
set(gca,'FontSize',18);

%% Plot psychometric curve with estimations
figure;
hold on;
for iPsychCurve = 1:nPsychCurves
    plot(2*ITDlin(ITDdB(5):ITDdB(2000)),p.psychCurves{iPsychCurve}.function(ITDdB(5):ITDdB(2000)),'Color',colors(1,:),'LineWidth',2);
    for iProcedure = 1:nProcedures
        errorbar(ps{iProcedure}.psychCurves{iPsychCurve}.JNDMean,ps{iProcedure}.psychCurves{iPsychCurve}.percentageCorrectMean,...
            ps{iProcedure}.psychCurves{iPsychCurve}.JNDMean - ps{iProcedure}.psychCurves{iPsychCurve}.JNDMin, ...
            ps{iProcedure}.psychCurves{iPsychCurve}.JNDMax - ps{iProcedure}.psychCurves{iPsychCurve}.JNDMean,...
            'horizontal','.','Color',colors(1+iProcedure,:),'LineWidth',2,'MarkerSize',30);
    end
end

hold off;
axis([10 2000 50 100])
set(gca, 'XScale', 'log')
xticks([10 50 100 500 1000])
xlabel('ITD [µs]');
ylabel('Percentage correct [%]');
legend({'Psychometric curve' labels{1:end}})
set(gca,'FontSize',18);

%% Plot histogramma (=sampling distributions)
figure; 
for iProcedure = 1:nProcedures
    subplot(nProcedures,1,iProcedure)
    hold on;
    for iPsychCurve = 1:nPsychCurves
        histogram(ps{iProcedure}.psychCurves{iPsychCurve}.convergence,'BinWidth',0.2,'FaceColor',colors(1+iProcedure,:));
    end
    plot([ps{iProcedure}.significanceLevel ps{iProcedure}.significanceLevel],[-100 2000],'-r','LineWidth',2);
    plot(ITDdB([50 300]),[0 0],'-k');
    xlim(ITDdB([50 300]))
    ylim([-100 2000])
    xticks(ITDdB([10 50 75 100 150 200 300 1000]))
    xticklabels(2*[10 50 75 100 150 200 300 1000])
    title(ps{iProcedure}.title);
    xlabel('Estimated JND [µs]');
    ylabel('Count');
    set(gca,'FontSize',18);
end
