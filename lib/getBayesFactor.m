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
switch length(parsMdl.rho)
    case 1
        rho1 = parsMdl.rho;
        rho2 = parsMdl.rho;
    case 2
        rho1 = parsMdl.rho(1);
        rho2 = parsMdl.rho(2);
end
OccamFactorInt = pi/ (parsMdl.Ls * parsMdl.LR * sqrt(parsMdl.TunKappa* sqrt(rho1*rho2)) * beta);
OccamFactorSeg = OccamFactorInt^2 * 4;

% Model evidence of segregation model
EvidSeg_Cue1 = CueGenMdl(x1, Lambda1, rho1, ParamMAP.StimEstim_Seg(1,:), ParamMAP.REstim_Seg(1,:), parsMdl);
EvidSeg_Cue2 = CueGenMdl(x2, Lambda2, rho2, ParamMAP.StimEstim_Seg(2,:), ParamMAP.REstim_Seg(2,:), parsMdl);
EvidSeg = EvidSeg_Cue1 .* EvidSeg_Cue2 .* OccamFactorSeg;

% Model evidence of integration model
EvidInt_Cue1 = CueGenMdl(x1, Lambda1, rho1, ParamMAP.StimEstim_Int, ParamMAP.REstim_Int, parsMdl);
EvidInt_Cue2 = CueGenMdl(x2, Lambda2, rho2, ParamMAP.StimEstim_Int, ParamMAP.REstim_Int, parsMdl);
EvidInt = EvidInt_Cue1 .* EvidInt_Cue2 .* OccamFactorInt;

% Bayes factor
BayesFactor_Cue1 = EvidSeg_Cue1./EvidInt_Cue1;
BayesFactor_Cue2 = EvidSeg_Cue2./EvidInt_Cue2;
BayesFactor = EvidSeg ./ EvidInt;
% BayesFactor = exp(log(EvidenceSeg) - log(EvidenceInt));

%% Fold the results into an output struct
MdlSelectRes.BayesFactor = BayesFactor;
MdlSelectRes.BayesFactor_Cue1 = BayesFactor_Cue1;
MdlSelectRes.BayesFactor_Cue2 = BayesFactor_Cue2;

MdlSelectRes.EvidenceSeg_Cue1 = EvidSeg_Cue1;
MdlSelectRes.EvidenceSeg_Cue2 = EvidSeg_Cue2;
MdlSelectRes.EvidenceInt_Cue1 = EvidInt_Cue1;
MdlSelectRes.EvidenceInt_Cue2 = EvidInt_Cue2;

MdlSelectRes.EvidenceSeg = EvidSeg;
MdlSelectRes.EvidenceInt = EvidInt;

MdlSelectRes.OccamFactorSeg = OccamFactorSeg;
MdlSelectRes.OccamFactorInt = OccamFactorInt;