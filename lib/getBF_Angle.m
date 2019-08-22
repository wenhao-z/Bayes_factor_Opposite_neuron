function MdlSelectRes = getBF_Angle(dat, ParamMAP, parsMdl, mode)
% Compute the evidence of integration model and segregation model,
% and Bayes factor

% Wen-Hao Zhang,
% Apr 25, 2019

% mode is an indicator:
% mode = 1: calculate the value of the Bayes factor in explaining the observations
% mode = 2: calculate the reconstructed likelihood function given best estimate of parameter

if ~exist('mode', 'var')
    mode = 1;
end
beta = parsMdl.beta;

x1 = dat(1,:);
x2 = dat(2,:);

% Occam factor
switch length(parsMdl.rho)
    case 1
        rho1 = parsMdl.rho;
        rho2 = parsMdl.rho;
    case 2
        rho1 = parsMdl.rho(1);
        rho2 = parsMdl.rho(2);
end

% Occam factor
OFInt = pi/ (parsMdl.Ls * parsMdl.LR * sqrt(parsMdl.TunKappa* sqrt(rho1*rho2)) * beta);
OFSeg = OFInt^2 * 4;

% Model evidence of segregation model
EvidSeg_x1 = DirGenMdl(x1, rho1, ParamMAP.StimEstim_Seg(1,:), ParamMAP.REstim_Seg(1,:), parsMdl);
EvidSeg_x2 = DirGenMdl(x2, rho2, ParamMAP.StimEstim_Seg(2,:), ParamMAP.REstim_Seg(2,:), parsMdl);
EvidSeg_x = EvidSeg_x1 .* EvidSeg_x2 * OFSeg;

% Model evidence of integration model
EvidInt_x1 = DirGenMdl(x1, rho1, ParamMAP.StimEstim_Int, ParamMAP.REstim_Int, parsMdl);
EvidInt_x2 = DirGenMdl(x2, rho2, ParamMAP.StimEstim_Int, ParamMAP.REstim_Int, parsMdl);
EvidInt_x = EvidInt_x1 .* EvidInt_x2 * OFInt;

% Bayes factor
LR_x1 = EvidSeg_x1./EvidInt_x1;
LR_x2 = EvidSeg_x2./EvidInt_x2;
BayesFactor = EvidSeg_x./ EvidInt_x;

%% Fold the results into an output struct
MdlSelectRes.BayesFactor = BayesFactor;
MdlSelectRes.LR_x1 = LR_x1;
MdlSelectRes.LR_x2 = LR_x2;

MdlSelectRes.EvidSeg_x1 = EvidSeg_x1;
MdlSelectRes.EvidSeg_x2 = EvidSeg_x2;
MdlSelectRes.EvidInt_x1 = EvidInt_x1;
MdlSelectRes.EvidInt_x2 = EvidInt_x2;

MdlSelectRes.EvidSeg_x = EvidSeg_x;
MdlSelectRes.EvidInt_x = EvidInt_x;

MdlSelectRes.OFSeg = OFSeg;
MdlSelectRes.OFInt = OFInt;
