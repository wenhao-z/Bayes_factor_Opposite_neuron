% Demo the whole likelihood function using the best fit parameters of two
% models

% Wen-Hao Zhang, Apr 10, 2019

% Set model parameters

% Load model parameters
% parsProbMdl;

parsMdl.Lambda = 30*[1; 1];
parsMdl.x = [0; 20] * pi/180; % Unit: rad

parsMdl.LR = 100; % The width of peak firing rate space

parsMdl.cueCond = 0; % 0: combined cue; 
%                      1: cue 1; 2: cue 2

switch parsMdl.cueCond
    case 1
        parsMdl.Lambda(2) = 0; % cue 1 condition, switch off cue 2
    case 2
        parsMdl.Lambda(1) = 0; % cue 2 condition, switch off cue 1
end

% Calculate the dependent parameters
parsMdl = getDependentPars(parsMdl);

%% Model evidence

% MAP estimate of model's parameters
ParamMAP = mapEstimFunc(parsMdl);

% The values of the reconstructed cues
x_recon = linspace(-pi, pi, parsMdl.N+1);
Lambda_recon = 0: parsMdl.LR + 3*sqrt(parsMdl.LR);
[x_recon, Lambda_recon] = ndgrid(x_recon, Lambda_recon);
datGrid = [x_recon(:), x_recon(:), Lambda_recon(:), Lambda_recon(:)]';

szDat = size(x_recon);

% Compute model evidence and Bayes factor
MdlSelectRes = getBayesFactor(datGrid, ParamMAP, parsMdl, 2);

for varName = fieldnames(MdlSelectRes)'
    if numel(MdlSelectRes.(varName{1})) >1 
        MdlSelectRes.(varName{1}) = reshape(MdlSelectRes.(varName{1}), szDat);
    end
end
datGrid = reshape(datGrid, [size(datGrid,1), szDat]);
clear varName

% Calculate the marginalized evidence and Bayes factor for moving direction
MdlSelect_Angle = getBayesFactor_Angle(x_recon(:,1), MdlSelectRes, ParamMAP);

%% Plot

figure
for iter = 1:6
   hAxe(iter) = subplot(2,3,iter);
   hold on
end

cLim = log([MdlSelectRes.EvidenceSeg(:); ...
    MdlSelectRes.EvidenceInt(:)]);
cLim = [min(cLim), max(cLim)];

nLines_Contour = 15;

% Evidence of the segregation model
if parsMdl.cueCond == 2
    imagesc(hAxe(1), x_recon(:,1)*180/pi, Lambda_recon(1,:), log(squeeze(MdlSelectRes.EvidenceSeg_Cue1))')
else
    contourf(hAxe(1), squeeze(x_recon*180/pi), squeeze(Lambda_recon), ...
        log(squeeze(MdlSelectRes.EvidenceSeg_Cue1)), nLines_Contour)
end
plot(hAxe(1), ParamMAP.StimEstim_Seg(1), ParamMAP.REstim_Seg(1)*parsMdl.beta, 'xk')
plot(hAxe(1), parsMdl.x*180/pi, parsMdl.Lambda, 'ok')
axes(hAxe(1)); axis square; axis xy
xlabel('x_1')
ylabel('\Lambda_1')
title('Model evidence (seg.)')
caxis(hAxe(1), cLim)
colorbar('north')

% Evidence of the integation model
contourf(hAxe(2), squeeze(x_recon*180/pi), squeeze(Lambda_recon), ...
    log(squeeze(MdlSelectRes.EvidenceInt_Cue1)), nLines_Contour)
plot(hAxe(2), ParamMAP.StimEstim_Int(1)*180/pi, ParamMAP.REstim_Int(1)*parsMdl.beta, 'xk')
plot(hAxe(2), parsMdl.x*180/pi, parsMdl.Lambda, 'ok')
% imagesc(hAxe(2), x1*180/pi, Lambda1, ...
%     log(squeeze(MdlSelectRes.EvidenceInt))')
axes(hAxe(2)); axis square; axis xy
xlabel('x_1')
title('Model evidence (int.)') 
caxis(hAxe(2), cLim)

% Bayes factor
if parsMdl.cueCond == 2
    imagesc(hAxe(3), x_recon(:,1)*180/pi, Lambda_recon(1,:), ...
        log(squeeze(MdlSelectRes.BayesFactor_Cue1))')
else
    contourf(hAxe(3), squeeze(x_recon*180/pi), squeeze(Lambda_recon), ...
        log(squeeze(MdlSelectRes.BayesFactor_Cue1)), nLines_Contour)
end
plot(hAxe(3), parsMdl.x*180/pi, parsMdl.Lambda, 'ok')
axes(hAxe(3)); axis square; axis xy
title('Bayes factor')
colorbar

% ------------------------------------------------------------------------
% Plot the marginalized evidence of segregation and integration model
plot(hAxe(4), x_recon(:,1)*180/pi, ...
    MdlSelect_Angle.EvidenceSeg_Cue1/sum(MdlSelect_Angle.EvidenceSeg_Cue1))
plot(hAxe(4), x_recon(:,1)*180/pi, ...
    MdlSelect_Angle.EvidenceSeg_Cue1_Theory/sum(MdlSelect_Angle.EvidenceSeg_Cue1_Theory), '--' );
axes(hAxe(4)); axis square; axis xy; 
legend('Numerical', 'Approx')

plot(hAxe(5), x_recon(:,1)*180/pi, ...
    MdlSelect_Angle.EvidenceInt_Cue1/sum(MdlSelect_Angle.EvidenceInt_Cue1))
plot(hAxe(5), x_recon(:,1)*180/pi, ...
    MdlSelect_Angle.EvidenceInt_Cue_Theory/sum(MdlSelect_Angle.EvidenceInt_Cue_Theory), '--');
axes(hAxe(5)); axis square; axis xy;
set(hAxe(4:5), 'ylim', [0, 0.15], 'ytick', 0:0.05:0.15, 'xlim', [-180, 180])

% ------------------------------------------------------------------------
% Plot the posterior distribution
postStim_int = vonMisesPdf(x_recon(:,1), ParamMAP.StimEstim_Int, ParamMAP.KappaEstim_Int);
postStim_int = postStim_int/sum(postStim_int);
plot(hAxe(6), x_recon(:,1)*180/pi, postStim_int);

% Plot the marginalized Bayes factor
plot(hAxe(6), x_recon(:,1)*180/pi, MdlSelect_Angle.BFx_VMFit);
plot(hAxe(6), x_recon(:,1)*180/pi, MdlSelect_Angle.BFx_Theory);
plot(hAxe(6), x_recon(:,1)*180/pi, MdlSelect_Angle.BFx_Naive, '--');

legend(hAxe(6), 'Post(s|D)', 'Bayes factor (VM fit)', 'Bayes factor (Theory)', ...
    'Bayes factor (naive)', 'location', 'best')
axes(hAxe(6)); axis square; axis xy; 
xlim([-180, 180])

% Set the other properties of figure
set(hAxe(1:3), 'ytick', 0:50:Lambda_recon(1,end), 'ylim', [0, 130]);
% 'ylim', [0, Lambda_recon(1,end)])
set(hAxe(1:6), 'xtick', -180:180:180, 'xlim', 180*[-1,1])

%% Plot the parallelogram of posterior and Bayes factor

vecPostSeg = ParamMAP.KappaEstim_Seg.* exp(1j*ParamMAP.StimEstim_Seg);
vecPostInt = ParamMAP.KappaEstim_Int.* exp(1j*ParamMAP.StimEstim_Int);

%Calculate the mean and concentration of BFx

BFx = MdlSelect_Angle.BFx_VMFit/ sum(MdlSelect_Angle.BFx_VMFit);
vecBFx = sum(BFx(:) .* exp(1j*x_recon(:,1))); % Mean resultant length
rhoBFx = abs(vecBFx);
vecBFx = mrl2Kappa(rhoBFx) .* exp(1j* angle(vecBFx));

figure; hold on
plot([0,real(vecPostSeg(1))], [0, imag(vecPostSeg(1))], 'g');
plot([0,real(vecPostSeg(2))], [0, imag(vecPostSeg(2))], 'g');
plot([0,real(vecPostInt)], [0, imag(vecPostInt)], 'b');
plot([0,real(vecBFx)], [0, imag(vecBFx)], 'r');

% Auxiliary lines
plot([0,real(vecPostSeg(1))]+real(vecPostSeg(2)), ...
    [0, imag(vecPostSeg(1))]+imag(vecPostSeg(2)), '--k');
plot([0,real(vecPostSeg(2))]+real(vecPostSeg(1)), ...
    [0, imag(vecPostSeg(2))]+imag(vecPostSeg(1)), '--k');
plot([0,real(vecBFx)]+real(vecPostInt)/2, ...
    [0, imag(vecBFx)]+imag(vecPostInt)/2, 'r');
axis equal
