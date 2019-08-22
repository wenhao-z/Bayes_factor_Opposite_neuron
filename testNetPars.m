% Test the performance of integration and Bayes factor in congruent and
% opposite neurons under different combinations of parameters.
% The results calculated by neurons are compared with theoretical predictions.

% Wen-Hao Zhang, May 10, 2019
% wenhao.zhang@pitt.edu

setWorkPath;

% Load model parameters
parsProbMdl;

% ------------------------------------------
% Set model parameters

parsMdl.tBin = 1/parsMdl.beta; % The time window of observing spikes

% Grid of stimulus direction
parsMdl.Stim = [zeros(1,10); 0: 10: 90] * pi/ 180; % unit: rad
% parsMdl.Stim = [zeros(1,1); 20] * pi/ 180; % unit: rad

% Grid of peak of tuning curves
R = (5:10:50) * parsMdl.tBin;
% R = 20 * parsMdl.tBin;
[R1, R2] = ndgrid(R,R);
parsMdl.R = [R1(:), R2(:)]';
clear R1 R2

parsMdl.nTrials = 50; % number of trials under each parameter s and R.
%                each trial will generate a feedforward spiking input,
%                corresponding to x and \Lambda

% Generate grid of parameters
[parGrid, dimPar] = paramGrid(parsMdl);


%%
[DecodeRes, TheoryRes] = simCongOpppoNet(parsMdl, 1);
DecodeRes = repmat(DecodeRes, size(parGrid));
TheoryRes = repmat(TheoryRes, size(parGrid));

tStart = clock;

% parpool(2);
parfor iterPar = 1: numel(parGrid)
    fprintf('Progress: %d/ %d\n', iterPar, numel(parGrid));
    
    mdlpars = parGrid(iterPar);
    
    % Perform the simulation for many trials
    [DecodeRes(iterPar), TheoryRes(iterPar)] = simCongOpppoNet(mdlpars);
end

tEnd = clock;
%% Save

savePath = fullfile(Path_RootDir, 'Data');
mkdir(savePath);

str = datestr(now, 'yymmddHHMM');
fileName = ['testNets_', str(1:6), ...
    '_', str(7:end) '.mat'];

save(fullfile(savePath, fileName), '-v7.3')

%% Comparison between the decoded results with theoretical predictions

for iter = 1: 7
    hAxe(iter) = subplot(2,4, iter);
    hold on
end

IdxStim = 1:7;
% Compare the posterior of integration model with the one decoded from congruent neurons
plot(hAxe(1), [TheoryRes(IdxStim,:).StimEstim_Int]*180/pi, [DecodeRes(IdxStim,:).meanCongCell]*180/pi, 'o')
plot(hAxe(2), [TheoryRes(IdxStim,:).KappaEstim_Int], [DecodeRes(IdxStim,:).KappaCongCell], 'o')
plot(hAxe(3), [TheoryRes(IdxStim,:).meanBFx1]*180/pi, [DecodeRes(IdxStim,:).meanOppoCell]*180/pi, 'o')
plot(hAxe(4), [TheoryRes(IdxStim,:).KappaBFx1], [DecodeRes(IdxStim,:).KappaOppoCell], 'o')

% Compare the likelihood ratio and Bayes factor with the one decoded from opposite neurons
BF_x = [TheoryRes(IdxStim,:).BF_x];
BF_x_Ocells = [DecodeRes(IdxStim,:).BF_x_Ocells];
OFR = [TheoryRes(IdxStim,:).OFR];
plot(hAxe(5), log(BF_x(1,:)), log(BF_x_Ocells(1,:)), 'o')
plot(hAxe(6), log(BF_x(2,:)), log(BF_x_Ocells(2,:)), 'o')
plot(hAxe(7), log(BF_x(1,:))+log(BF_x(2,:))+log(OFR), ...
    log(BF_x_Ocells(1,:))+log(BF_x_Ocells(2,:))+log(OFR), 'o')

% Plot diagonal lines
for iter = 1: 7
    axisLim = axis(hAxe(iter));
    axisLim = ceil(abs(axisLim)) .* sign(axisLim);
    axisLim([1,3]) = min(axisLim([1,3]));
    axisLim([2,4]) = min(axisLim([2,4]));
    plot(hAxe(iter), axisLim(1:2), axisLim(3:4), '--k')
    axis(hAxe(iter), axisLim)
end

axis(hAxe, 'square')
xlabel(hAxe(1), 'Post. mean(theory)')
x2 = parsMdl.Stim(2,IdxStim([1,end]))*180/pi;
title(hAxe(1), sprintf('x_1=0, x_2=[%d, %2.0f]', x2(1), x2(2)))
ylabel(hAxe(1), 'Post. mean(neuron)')
xlabel(hAxe(2), 'Post. conc.(theory)')
ylabel(hAxe(2), 'Post. conc.(neuron)')
xlabel(hAxe(3), 'BF mean(num.)')
ylabel(hAxe(3), 'BF mean(neuron)')
xlabel(hAxe(4), 'BF conc.(num.)')
ylabel(hAxe(4), 'BF conc.(neuron)')

xlabel(hAxe(5), 'LR(x_1) (num.)')
ylabel(hAxe(5), 'LR(x_1) (neuron.)')
xlabel(hAxe(6), 'LR(x_2) (num.)')
ylabel(hAxe(6), 'LR(x_2) (neuron.)')
xlabel(hAxe(7), 'BF(x_1,x_2) (num.)')
ylabel(hAxe(7), 'BF(x_1,x_2) (neuron.)')


