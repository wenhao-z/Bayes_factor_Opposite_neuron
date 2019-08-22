% Calculate the Bayes factor

% Wen-Hao Zhang,
% Apr 8 ,2019

% Set model parameters

% Load model parameters
% parsProbMdl;

parsMdl.x = [0; 0]*pi/180; % Unit: deg

Lambda1 = 0:100;
Lambda2 = 0:100;
[Lambda1_grid, Lambda2_grid] = ndgrid(Lambda1, Lambda2);
parsMdl.Lambda = [Lambda1_grid(:), Lambda2_grid(:)]';

%% Stimulus estimate (MAP) of two models

% MAP estimate of model's parameters
ParamMAP = mapEstimFunc(parsMdl);

x1_grid = parsMdl.x(1,:)*ones(size(Lambda1_grid));
x2_grid = parsMdl.x(2,:)*ones(size(Lambda1_grid));
Lambda1_grid = parsMdl.Lambda(1,:);
Lambda2_grid = parsMdl.Lambda(2,:);

dat = [x1_grid(:), x2_grid(:), Lambda1_grid(:), Lambda2_grid(:)]';

MdlSelectRes = getBayesFactor(dat, ParamMAP, parsMdl);

for varName = fieldnames(MdlSelectRes)'
    if numel(MdlSelectRes.(varName{1})) > 1
        MdlSelectRes.(varName{1}) = reshape(MdlSelectRes.(varName{1}), length(Lambda1), length(Lambda2));
    end
end
clear varName

%% Plot
% The model evidence and Bayes factor with cue inputs

% Set the colormap
logBF = log(MdlSelectRes.BayesFactor);
cmap = getColorMapPosNegDat([min(logBF(:)), max(logBF(:))]);

figure
colormap(cmap)
subplot(1,2,1)
imagesc(Lambda1_grid, Lambda2_grid, log(MdlSelectRes.BayesFactor))
% contourf(reshape(Lambda1_grid, length(Lambda1), []), ...
%     reshape(Lambda2_grid, length(Lambda2), []), ...
%     log(MdlSelectRes.BayesFactor));
axis square; axis xy
set(gca, 'xtick', Lambda1([1, (end+1)/2, end]), 'ytick', Lambda2([1, (end+1)/2, end]))
title('Log. of Bayes factor')
xlabel('\Lambda_1')
ylabel('\Lambda_2')
colorbar

subplot(1,2,2)
dLambda = diag(rot90(reshape(Lambda2_grid, length(Lambda1), []))) ...
    - diag(rot90(reshape(Lambda1_grid, length(Lambda1), [])));
plot(dLambda, log(diag(rot90(MdlSelectRes.EvidenceSeg))))
hold on
plot(dLambda, log(diag(rot90(MdlSelectRes.EvidenceInt))))
plot(dLambda, log(diag(rot90(MdlSelectRes.BayesFactor))))
plot(100*[-1, 1], zeros(1,2), '--k')
axis square
% set(gca, 'xtick', -90:90:90, 'xlim', 90*[-1,1])
legend('Seg. mdl.', 'Int. mdl.', 'Bayes factor', 'location', 'best')
ylabel('Model Evidence and Bayes Factor (log)')
xlabel('|\Lambda_1-\Lambda_2|')
