% Decode the posterior and Bayes factor from congruent and opposite
% neurons, and compare with theoretical predictions.

% Wen-Hao Zhang, April 22, 2019
% wenhao.zhang@pitt.edu

% Load model parameters
parsProbMdl;

% Set model parameters
parsMdl.Stim = [0, 30]*pi/180; % unit: rad
parsMdl.R = [30, 30]; % Peak of tuning curves
parsMdl.tBin = 1/parsMdl.beta; % The time window of observing spikes

parsMdl.R = parsMdl.R * parsMdl.tBin; % This makes the mean spike count in the time window is R

%% Generate feedforward inputs from two sensory modalities (Poisson spikes)

TuningFF = @(s, R) R*exp(parsMdl.TunKappa * (cos(s - parsMdl.PrefStim)-1) );

RateFF = [...
    TuningFF(parsMdl.Stim(1), parsMdl.R(1)), ...
    TuningFF(parsMdl.Stim(2), parsMdl.R(2))];

%%
nTrials = 100;
% BF_x = zeros(2, nTrials);
% BF_x_Ocells = zeros(2, nTrials);

OFR = 4*pi/ (parsMdl.Ls * parsMdl.LR * sqrt(parsMdl.TunKappa) * parsMdl.beta);
OFR = repmat(OFR, 1, nTrials);

EstimNum = struct(...
    'meanBFx1', zeros(1, nTrials), ...
    'KappaBFx1', zeros(1, nTrials), ...
    'meanBFx2', zeros(1, nTrials), ...
    'KappaBFx2', zeros(1, nTrials), ...
    'BF_x', zeros(2, nTrials));
EstimTheory = repmat(mapEstimFunc([]), 1, nTrials);


DecodeRes = struct(...
    'meanCongCell', zeros(1, nTrials), ...
    'KappaCongCell', zeros(1, nTrials), ...
    'meanOppoCell', zeros(1, nTrials), ...
    'KappaOppoCell', zeros(1, nTrials), ...
    'BF_x_Ocells', zeros(2, nTrials));

x_grid = parsMdl.PrefStim;

for iterTrial = 1: nTrials
    SpkFF = poissrnd(RateFF);
    
    % Decode the position and concentration from feedforward inputs
    [parsMdl.x, parsMdl.Kappa] = popVecDecoder(SpkFF, parsMdl);
    parsMdl.Lambda = sum(SpkFF, 1)';
    parsMdl = getDependentPars(parsMdl, SpkFF);
    
    % Theoretical prediction of Bayes factor of direction
    % MAP estimate of model's parameters
    ParamMAP = mapEstimFunc(parsMdl);
    
    EstimTheory(iterTrial) = ParamMAP;
    
    %% Marginalized model evidence of moving direction (angle)
    % Numerically marginalize the dimension of Lambda
    
    % ------------- cue 1 --------------------
    EvidInt_x1 = DirGenMdl(x_grid, parsMdl.rho(1), ParamMAP.StimEstim_Int, ParamMAP.REstim_Int, parsMdl);
    EvidSeg_x1 = DirGenMdl(x_grid, parsMdl.rho(1), ParamMAP.StimEstim_Seg(1), ParamMAP.REstim_Seg(1), parsMdl);
    
    % Fit the evidence by using a von Mises distribution first
    EvdSeg = sum(EvidSeg_x1(:) .* exp(1j * parsMdl.PrefStim(:)))/ sum(EvidSeg_x1);
    EvdInt = sum(EvidInt_x1(:) .* exp(1j * parsMdl.PrefStim(:)))/ sum(EvidInt_x1);
    
    % Parameter of Bayes factor
    BF_xVec = mrl2Kappa(abs(EvdSeg)) * exp(1j*angle(EvdSeg)) - ...
        mrl2Kappa(abs(EvdInt)) * exp(1j*angle(EvdInt));
    EstimNum.meanBFx1(iterTrial) = angle(BF_xVec);
    EstimNum.KappaBFx1(iterTrial) = abs(BF_xVec);
    
    % Value of Bayes factor
    EstimNum.BF_x(1, iterTrial) = vonMisesPdf(parsMdl.x(1), angle(EvdSeg), mrl2Kappa(abs(EvdSeg))) ./ ...
        vonMisesPdf(parsMdl.x(1), angle(EvdInt), mrl2Kappa(abs(EvdInt)));
    %     BF_x(1,iterTrial) = vonMisesPdf(parsMdl.x(1), angle(BF_xVec), abs(BF_xVec));
    
    % ------------- cue 2 --------------------
    EvidInt_x2 = DirGenMdl(x_grid, parsMdl.rho(2), ParamMAP.StimEstim_Int, ParamMAP.REstim_Int, parsMdl);
    EvidSeg_x2 = DirGenMdl(x_grid, parsMdl.rho(2), ParamMAP.StimEstim_Seg(2), ParamMAP.REstim_Seg(2), parsMdl);
    
    % Fit the evidence by using a von Mises distribution first
    EvdSeg = sum(EvidSeg_x2(:) .* exp(1j * parsMdl.PrefStim(:)))/ sum(EvidSeg_x2);
    EvdInt = sum(EvidInt_x2(:) .* exp(1j * parsMdl.PrefStim(:)))/ sum(EvidInt_x2);
    
    % Parameter of Bayes factor
    BF_xVec = mrl2Kappa(abs(EvdSeg)) * exp(1j*angle(EvdSeg)) - ...
        mrl2Kappa(abs(EvdInt)) * exp(1j*angle(EvdInt));
    EstimNum.meanBFx2(iterTrial) = angle(BF_xVec);
    EstimNum.KappaBFx2(iterTrial) = abs(BF_xVec);
    
    % Value of Bayes factor
    EstimNum.BF_x(2, iterTrial) = vonMisesPdf(parsMdl.x(2), angle(EvdSeg), mrl2Kappa(abs(EvdSeg))) ./ ...
        vonMisesPdf(parsMdl.x(2), angle(EvdInt), mrl2Kappa(abs(EvdInt)));
    %     BF_x(2,iterTrial) = vonMisesPdf(parsMdl.x(2), angle(BF_xVec), abs(BF_xVec));
    
    % Occam factor ratio
    OFR(iterTrial) = OFR(iterTrial)/ (parsMdl.rho(1)*parsMdl.rho(2))^(1/4);
    
    %% Neuronal response
    SpkCongCell = sum(SpkFF,2);
    SpkOppoCell = SpkFF(:,1) + [SpkFF(end/2+1:end,2); SpkFF(1:end/2,2)];
    SpkOppoCell = (SpkOppoCell/2); % This might be a problem of dividing spike count by 2
    
    % Decode the distribution from congruent and opposite neurons
    [meanCongCell, KappaCongCell] = popVecDecoder(SpkCongCell, parsMdl);
    [meanOppoCell, KappaOppoCell] = popVecDecoder(SpkOppoCell, parsMdl);
    
    NormConst1 = 2*pi*besseli(0, KappaCongCell/2)*besseli(0, KappaOppoCell)/besseli(0, ParamMAP.KappaEstim_Seg(1));
    NormConst2 = 2*pi*besseli(0, KappaCongCell/2)*besseli(0, KappaOppoCell)/besseli(0, ParamMAP.KappaEstim_Seg(2));
    DecodeRes.BF_x_Ocells(1,iterTrial) = vonMisesPdf(parsMdl.x(1), meanOppoCell, KappaOppoCell) * NormConst1;
    DecodeRes.BF_x_Ocells(2,iterTrial) = vonMisesPdf(parsMdl.x(2), meanOppoCell+pi, KappaOppoCell) * NormConst2;
    
    DecodeRes.meanCongCell(iterTrial) = meanCongCell;
    DecodeRes.KappaCongCell(iterTrial) = KappaCongCell;
    DecodeRes.meanOppoCell(iterTrial) = meanOppoCell;
    DecodeRes.KappaOppoCell(iterTrial) = KappaOppoCell;
    
    %% Plot
    %     cla; plotBayesFactor; pause;
end
clear BF_2Vec meanCongCell KappaCongCell meanOppoCell KappaOppoCell
clear EvidInt_x1 EvidSeg_x1 EvidInt_x2 EvidSeg_x2 EvdSeg EvdInt

%% Comparison between the decoded results with theoretical predictions

for iter = 1: 7
    hAxe(iter) = subplot(2,4, iter);
    hold on
end
% for iter = 1:3
%     hAxe(iter+4) = subplot(2,3,iter+3);
%     hold on
% end
% Compare the posterior of integration model with the one decoded from congruent neurons
plot(hAxe(1), [EstimTheory.StimEstim_Int]*180/pi, [DecodeRes.meanCongCell]*180/pi, 'o')
plot(hAxe(2), [EstimTheory.KappaEstim_Int], [DecodeRes.KappaCongCell], 'o')
plot(hAxe(3), [EstimNum.meanBFx1]*180/pi, [DecodeRes.meanOppoCell]*180/pi, 'o')
plot(hAxe(4), [EstimNum.KappaBFx1], [DecodeRes.KappaOppoCell], 'o')

% Compare the likelihood ratio and Bayes factor with the one decoded from opposite neurons
plot(hAxe(5), log(EstimNum.BF_x(1,:)), log(DecodeRes.BF_x_Ocells(1,:)), 'o')
plot(hAxe(6), log(EstimNum.BF_x(2,:)), log(DecodeRes.BF_x_Ocells(2,:)), 'o')
plot(hAxe(7), log(EstimNum.BF_x(1,:))+log(EstimNum.BF_x(2,:))+log(OFR), ...
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


