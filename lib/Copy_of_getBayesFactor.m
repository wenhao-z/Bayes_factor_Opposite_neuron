function MdlSelectRes = getBayesFactor(dat, ParamMAP, parsMdl, mode)
% Compute the evidence of integration model and segregation model,
% and Bayes factor

% Wen-Hao Zhang,
% Apr 9, 2019

% mode is an indicator:
% mode = 1: calculate the value of the Bayes factor in explaining the observations
% mode = 2: calculate the reconstructed likelihood function given best estimate of parameter

if ~exist('mode', 'var')
    mode = 1;
end
beta = parsMdl.beta;

x1 = dat(1,:);
x2 = dat(2,:);
Lambda1 = dat(3,:);
Lambda2 = dat(4,:);

% Occam factor
OccamFactorInt = pi/ (parsMdl.Ls * parsMdl.LR * sqrt(parsMdl.TunKappa*parsMdl.rho) * beta);
OccamFactorSeg = OccamFactorInt^2 * 4;

% Model evidence of segregation model
EvidenceSeg_Cue1 = CueGenMdl(x1, Lambda1, ParamMAP.StimEstim_Seg(1,:), ParamMAP.REstim_Seg(1,:), parsMdl);
EvidenceSeg_Cue2 = CueGenMdl(x2, Lambda2, ParamMAP.StimEstim_Seg(2,:), ParamMAP.REstim_Seg(2,:), parsMdl);
EvidenceSeg = EvidenceSeg_Cue1 .* EvidenceSeg_Cue2 .* OccamFactorSeg;

% Model evidence of integration model
EvidenceInt_Cue1 = CueGenMdl(x1, Lambda1, ParamMAP.StimEstim_Int, ParamMAP.REstim_Int, parsMdl);
EvidenceInt_Cue2 = CueGenMdl(x2, Lambda2, ParamMAP.StimEstim_Int, ParamMAP.REstim_Int, parsMdl);
EvidenceInt = EvidenceInt_Cue1 .* EvidenceInt_Cue2 .* OccamFactorInt;

% Bayes factor
BayesFactor_Cue1 = EvidenceSeg_Cue1./EvidenceInt_Cue1;
BayesFactor_Cue2 = EvidenceSeg_Cue2./EvidenceInt_Cue2;
BayesFactor = EvidenceSeg ./ EvidenceInt;
% BayesFactor = exp(log(EvidenceSeg) - log(EvidenceInt));

%% Fold the results into an output struct
MdlSelectRes.BayesFactor = BayesFactor;
MdlSelectRes.BayesFactor_Cue1 = BayesFactor_Cue1;
MdlSelectRes.BayesFactor_Cue2 = BayesFactor_Cue2;

MdlSelectRes.EvidenceSeg_Cue1 = EvidenceSeg_Cue1;
MdlSelectRes.EvidenceSeg_Cue2 = EvidenceSeg_Cue2;
MdlSelectRes.EvidenceInt_Cue1 = EvidenceInt_Cue1;
MdlSelectRes.EvidenceInt_Cue2 = EvidenceInt_Cue2;

MdlSelectRes.EvidenceSeg = EvidenceSeg;
MdlSelectRes.EvidenceInt = EvidenceInt;

MdlSelectRes.OccamFactorSeg = OccamFactorSeg;
MdlSelectRes.OccamFactorInt = OccamFactorInt;