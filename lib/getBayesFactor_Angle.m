function MdlSelect_Angle = getBayesFactor_Angle(x1, MdlSelectRes, ParamMAP)
% Calculate the marginalized model evidence and Bayes factor over moving
% direction (angle)

% Wen-Hao Zhang
% April 17, 2019


% Sum the 3rd and 4th dim of the variables in MdlSelectRes in order to
% marginalize the dimension corresponding to Lamba1 and Lambda2
MdlSelect_Angle = MdlSelectRes;
for varName = fieldnames(MdlSelect_Angle)'
    if numel(MdlSelect_Angle.(varName{1})) <= 1
       MdlSelect_Angle = rmfield(MdlSelect_Angle, varName{1});
       continue
    end
    if ~strcmp(varName{1}, 'BayesFactor')
%         MdlSelect_Angle.(varName{1}) = sum(sum(MdlSelect_Angle.(varName{1}), 3), 4);
        MdlSelect_Angle.(varName{1}) = sum(MdlSelect_Angle.(varName{1}), 2);
    end
end

% Calculate the approximated marginalized model evidence and Bayes factor
MdlSelect_Angle.EvidenceInt_Cue_Theory = vonMisesPdf(x1, ParamMAP.StimEstim_Int, ParamMAP.KappaEstim_Int/2);
MdlSelect_Angle.EvidenceSeg_Cue1_Theory = vonMisesPdf(x1, ParamMAP.StimEstim_Seg(1), ParamMAP.KappaEstim_Seg(1));
MdlSelect_Angle.EvidenceSeg_Cue2_Theory = vonMisesPdf(x1, ParamMAP.StimEstim_Seg(2), ParamMAP.KappaEstim_Seg(2));

%% Bayes factor of cue 1
% CAUTION: Naive way in calculating Bayes factor. Large numerical error!!
BFx_Naive = MdlSelect_Angle.EvidenceSeg_Cue1 ./ MdlSelect_Angle.EvidenceInt_Cue1;
MdlSelect_Angle.BFx_Naive = BFx_Naive/sum(BFx_Naive);

% Fit the model evidence by a von Mises distribution first, and divide two
% von Mises distributions
EvdSeg = MdlSelect_Angle.EvidenceSeg_Cue1 / sum(MdlSelect_Angle.EvidenceSeg_Cue1);
EvdInt = MdlSelect_Angle.EvidenceInt_Cue1 / sum(MdlSelect_Angle.EvidenceInt_Cue1);

EvdSeg = sum(EvdSeg(:) .* exp(1j * x1(:)));
EvdInt = sum(EvdInt(:) .* exp(1j * x1(:)));

BFx_VMFit = vonMisesPdf(x1, angle(EvdSeg), mrl2Kappa(abs(EvdSeg))) ./ ...
    vonMisesPdf(x1, angle(EvdInt), mrl2Kappa(abs(EvdInt)));
MdlSelect_Angle.BFx_VMFit = BFx_VMFit/sum(BFx_VMFit);

% Theoretical approximation of marginalized Bayes factor
% diffVec = parsMdl.Kappa .* exp(1j * parsMdl.x)/2;
% diffVec = -diff(diffVec);

diffVec = ParamMAP.KappaEstim_Seg(1) .* exp(1j*ParamMAP.StimEstim_Seg(1)) ...
    -ParamMAP.KappaEstim_Int(1) .* exp(1j*ParamMAP.StimEstim_Int(1))/2;
BFx_Theory = vonMisesPdf(x1, angle(diffVec), abs(diffVec));
MdlSelect_Angle.BFx_Theory = BFx_Theory/sum(BFx_Theory);

clear diffVec
