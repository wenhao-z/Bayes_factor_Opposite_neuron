% Find the decision boundary between integration and segregation with input
% spike count

% Wen-Hao Zhang,
% Apr 25,2019

% Set model parameters

% Load model parameters
parsProbMdl;

x1 = 0;
x2 = linspace(0, pi, 2*parsMdl.N+1);
[x1_grid, x2_grid] = ndgrid(x1, x2);
parsMdl.x = [x1_grid(:), x2_grid(:)]';

LambdaArray = 1:1:150;
logBFArray = zeros(length(x2), length(LambdaArray));

%%
for iter = 1: length(LambdaArray)
    parsMdl.Lambda = LambdaArray(iter)*[1; 1]; % Spike count of feedforward inputs. Unit: Hz
    
    % Stimulus estimate (MAP) of two models
    % MAP estimate of model's parameters
    ParamMAP = mapEstimFunc(parsMdl);
    
    x1_grid = parsMdl.x(1,:);
    x2_grid = parsMdl.x(2,:);
    dat = [x1_grid(:), x2_grid(:)]';
    
    MdlSelectRes = getBF_Angle(dat, ParamMAP, parsMdl);
    
    for varName = fieldnames(MdlSelectRes)'
        if numel(MdlSelectRes.(varName{1})) > 1
            MdlSelectRes.(varName{1}) = reshape(MdlSelectRes.(varName{1}), length(x1), length(x2));
        end
    end
    clear varName

    logBFArray(:,iter) = log(MdlSelectRes.BayesFactor);
end
% Find the boundary between integration and segregation
[~, Idx] = min(abs(logBFArray),[], 1);
DecisionBoundary = x2(Idx);
DecisionBoundary(DecisionBoundary<0) = DecisionBoundary(DecisionBoundary<0) + pi/2;

%% Plot

% The model evidence and Bayes factor with cue inputs
% subplot(1,2,1)
cmap = getColorMapPosNegDat([min(logBFArray(:)), max(logBFArray(:))]);
imagesc(LambdaArray, x2*180/pi, logBFArray)
hold on
plot(LambdaArray, DecisionBoundary*180/pi)
xlabel('Input spike count \Lambda')
ylabel('Direction disparity |x_1-x_2|')
set(gca, 'ylim', [0, 180], 'ytick', 0:45:180, 'xlim', [0, 150], 'xtick', 0:50:150)
axis xy; axis square
colorbar
colormap(cmap)



