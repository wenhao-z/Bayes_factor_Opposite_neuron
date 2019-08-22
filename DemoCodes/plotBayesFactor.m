% An intermediate code comparing the normalized Bayes factor computed
% through a varieties of ways


%% Call the packed functions which are used in the calculation of simple models
% The values of the reconstructed cues
x_recon = parsMdl.PrefStim;
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

%% The core code in the packed functions, which should be yield the same results
x_grid = parsMdl.PrefStim;
Lambda_grid = 0: parsMdl.LR + 3*sqrt(parsMdl.LR);
[x_grid, Lambda_grid] = ndgrid(x_grid, Lambda_grid);

% Compute the model evidence
EvidSeg_x1 = CueGenMdl(x_grid, Lambda_grid, parsMdl.rho(1), ParamMAP.StimEstim_Seg(1), ParamMAP.REstim_Seg(1), parsMdl);
EvidSeg_x1 = sum(EvidSeg_x1, 2);
EvidInt_x1 = CueGenMdl(x_grid, Lambda_grid, parsMdl.rho(1), ParamMAP.StimEstim_Int, ParamMAP.REstim_Int, parsMdl);
EvidInt_x1 = sum(EvidInt_x1, 2);

EvdSeg = sum(EvidSeg_x1(:) .* exp(1j * parsMdl.PrefStim(:)))/ sum(EvidSeg_x1);
EvdInt = sum(EvidInt_x1(:) .* exp(1j * parsMdl.PrefStim(:)))/ sum(EvidInt_x1);

BF_x1 = ...
    vonMisesPdf(x_grid, angle(EvdSeg), mrl2Kappa(abs(EvdSeg))) ./ ...
    vonMisesPdf(x_grid, angle(EvdInt), mrl2Kappa(abs(EvdInt)));

BF_x11 = mrl2Kappa(abs(EvdSeg)) * exp(1j*angle(EvdSeg)) - ...
    mrl2Kappa(abs(EvdInt)) * exp(1j*angle(EvdInt));
BF_x11 = vonMisesPdf(x_grid, angle(BF_x11), abs(BF_x11));

%% Decoded distributions from simple neuronal responses
diffVec = ParamMAP.KappaEstim_Seg(1) .* exp(1j*ParamMAP.StimEstim_Seg(1)) ...
    -ParamMAP.KappaEstim_Int(1) .* exp(1j*ParamMAP.StimEstim_Int(1))/2;
% diffVec1 = -diff(parsMdl.Kappa.*exp(1j*parsMdl.x))/2;
% BF_Ocells = vonMisesPdf(x_grid, angle(diffVec), abs(diffVec));

BF_Ocells = vonMisesPdf(x_grid, meanOppoCell, KappaOppoCell);

%% Plot
figure(1);
cla;
hold on
cSpec = lines(2);
plot(parsMdl.PrefStim*180/pi, BF_x1/sum(BF_x1), 'color', cSpec(1,:))
plot(parsMdl.PrefStim*180/pi, BF_x11/sum(BF_x11), '--o', 'color', cSpec(1,:))
plot(parsMdl.PrefStim*180/pi, BF_Ocells/sum(BF_Ocells), 'color', cSpec(2,:))
plot(x_recon(:,1)*180/pi, MdlSelect_Angle.BFx_Theory, '--o', 'color', cSpec(2,:))
