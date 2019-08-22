% Calculate the Bayes factor

% Wen-Hao Zhang,
% Apr 8 ,2019

% Set model parameters

% Load model parameters
% parsProbMdl;

x1 = linspace(-pi/2, pi/2, parsMdl.N/2+1);
x2 = linspace(-pi/2, pi/2, parsMdl.N/2+1);
[x1_grid, x2_grid] = ndgrid(x1, x2);
parsMdl.x = [x1_grid(:), x2_grid(:)]';

parsMdl.Lambda = 30*[1; 1]; % Spike count of feedforward inputs. Unit: Hz

%% Stimulus estimate (MAP) of two models

% MAP estimate of model's parameters
ParamMAP = mapEstimFunc(parsMdl);

x1_grid = parsMdl.x(1,:);
x2_grid = parsMdl.x(2,:);
Lambda1_grid = parsMdl.Lambda(1)*ones(size(x1_grid));
Lambda2_grid = parsMdl.Lambda(2)*ones(size(x1_grid));

dat = [x1_grid(:), x2_grid(:), Lambda1_grid(:), Lambda2_grid(:)]';

MdlSelectRes = getBayesFactor(dat, ParamMAP, parsMdl);

for varName = fieldnames(MdlSelectRes)'
    if numel(MdlSelectRes.(varName{1})) > 1
        MdlSelectRes.(varName{1}) = reshape(MdlSelectRes.(varName{1}), length(x1), length(x2));
    end
end
clear varName

%% Plot
% The model evidence and Bayes factor with cue inputs
logBF = log(MdlSelectRes.BayesFactor);
cmap = getColorMapPosNegDat([min(logBF(:)), max(logBF(:))]);

figure
colormap(cmap)

subplot(1,2,1)
imagesc(x1*180/pi, x2*180/pi, log(MdlSelectRes.BayesFactor))
% contourf(reshape(x1_grid, length(x1), []) * 180/pi, ...
%     reshape(x2_grid, length(x2), []) * 180/pi, ...
%     log(MdlSelectRes.BayesFactor));
axis square; axis xy
set(gca, 'xtick', x1([1, (end+1)/2, end])*180/pi, 'ytick', x2([1, (end+1)/2, end])*180/pi)
title('Log. of Bayes factor')
xlabel('x_1')
ylabel('x_2')
colorbar 

subplot(1,2,2)
dAngle = diag(rot90(reshape(x2_grid, length(x1), []))) ...
    - diag(rot90(reshape(x1_grid, length(x1), [])));
dAngle = dAngle * 180/pi;
plot(dAngle, log(diag(rot90(MdlSelectRes.EvidenceSeg))))
hold on
plot(dAngle, log(diag(rot90(MdlSelectRes.EvidenceInt))))
plot(dAngle, log(diag(rot90(MdlSelectRes.BayesFactor))))
plot(180*[-1, 1], zeros(1,2), '--k')
axis square
set(gca, 'xtick', -90:90:90, 'xlim', 90*[-1,1])
legend('Seg. mdl.', 'Int. mdl.', 'Bayes factor', 'location', 'best')
ylabel('Model Evidence and Bayes Factor (log)')
xlabel('|x_1-x_2|')
