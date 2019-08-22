function [DecodeRes, TheoryRes] = simCongOpppoNet(parsMdl, bInit)
% Simulate a simple network consisting of congruent and opposite neurons
% for many trials. And the theoretical predictions will also be calculated
% for each trial.

% Wen-Hao Zhang
% May 20, 2019

if ~exist('bInit', 'var')
    bInit = 0; % Default value
    % bInit = 0: run simulation
    % bInit = 1: return an initialized structs DecodeRes and TheoryRes
end
%% Generate feedforward inputs from two sensory modalities (Poisson spikes)
TuningFF = @(s, R) R*exp(parsMdl.TunKappa * (cos(s - parsMdl.PrefStim)-1) );

RateFF = [...
    TuningFF(parsMdl.Stim(1), parsMdl.R(1)), ...
    TuningFF(parsMdl.Stim(2), parsMdl.R(2))];

%% Initialze arrays
nTrials = parsMdl.nTrials;

OFR = 4*pi/ (parsMdl.Ls * parsMdl.LR * sqrt(parsMdl.TunKappa) * parsMdl.beta);
OFR = repmat(OFR, 1, nTrials);

TheoryRes = struct(...
    'meanBFx1', zeros(1, nTrials), ...
    'KappaBFx1', zeros(1, nTrials), ...
    'meanBFx2', zeros(1, nTrials), ...
    'KappaBFx2', zeros(1, nTrials), ...
    'BF_x', zeros(2, nTrials), ...
    'OFR', OFR);
EstimTheory = repmat(mapEstimFunc([]), 1, nTrials);

DecodeRes = struct(...
    'meanCongCell', zeros(1, nTrials), ...
    'KappaCongCell', zeros(1, nTrials), ...
    'meanOppoCell', zeros(1, nTrials), ...
    'KappaOppoCell', zeros(1, nTrials), ...
    'BF_x_Ocells', zeros(2, nTrials));

if bInit
    for varName = fieldnames(EstimTheory)'
        lenFields = length(EstimTheory(1).(varName{1}));
        NetStat = [EstimTheory.(varName{1})];
        NetStat = reshape(NetStat, lenFields, []);
        TheoryRes.(varName{1}) = squeeze(NetStat); % First several dims are the same as parGrids. The last few dims are the same as varName
    end

   return; % return an empty DecodeRes and TheoryRes for initialization
end
x_grid = parsMdl.PrefStim;

%%
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
    TheoryRes.meanBFx1(iterTrial) = angle(BF_xVec);
    TheoryRes.KappaBFx1(iterTrial) = abs(BF_xVec);
    
    % Value of Bayes factor
    TheoryRes.BF_x(1, iterTrial) = vonMisesPdf(parsMdl.x(1), angle(EvdSeg), mrl2Kappa(abs(EvdSeg))) ./ ...
        vonMisesPdf(parsMdl.x(1), angle(EvdInt), mrl2Kappa(abs(EvdInt)));
    
    % ------------- cue 2 --------------------
    EvidInt_x2 = DirGenMdl(x_grid, parsMdl.rho(2), ParamMAP.StimEstim_Int, ParamMAP.REstim_Int, parsMdl);
    EvidSeg_x2 = DirGenMdl(x_grid, parsMdl.rho(2), ParamMAP.StimEstim_Seg(2), ParamMAP.REstim_Seg(2), parsMdl);
    
    % Fit the evidence by using a von Mises distribution first
    EvdSeg = sum(EvidSeg_x2(:) .* exp(1j * parsMdl.PrefStim(:)))/ sum(EvidSeg_x2);
    EvdInt = sum(EvidInt_x2(:) .* exp(1j * parsMdl.PrefStim(:)))/ sum(EvidInt_x2);
    
    % Parameter of Bayes factor
    BF_xVec = mrl2Kappa(abs(EvdSeg)) * exp(1j*angle(EvdSeg)) - ...
        mrl2Kappa(abs(EvdInt)) * exp(1j*angle(EvdInt));
    TheoryRes.meanBFx2(iterTrial) = angle(BF_xVec);
    TheoryRes.KappaBFx2(iterTrial) = abs(BF_xVec);
    
    % Value of Bayes factor
    TheoryRes.BF_x(2, iterTrial) = vonMisesPdf(parsMdl.x(2), angle(EvdSeg), mrl2Kappa(abs(EvdSeg))) ./ ...
        vonMisesPdf(parsMdl.x(2), angle(EvdInt), mrl2Kappa(abs(EvdInt)));
    
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
end

%% List of output arguments

for varName = fieldnames(EstimTheory)'
    lenFields = length(EstimTheory(1).(varName{1}));
    NetStat = [EstimTheory.(varName{1})];
    NetStat = reshape(NetStat, lenFields, []);
    TheoryRes.(varName{1}) = squeeze(NetStat); % First several dims are the same as parGrids. The last few dims are the same as varName
end
