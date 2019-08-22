% An interface doing all demonstrations of Bayes factor
% Wen-Hao Zhang, 
% April 17, 2019

setWorkPath;

% Load default model parameters
parsProbMdl;


%%
flagDemo = 2;
% 1. Plot the Bayes factor as varing moving DIRECTIONS of two cues
% 2. Plot the Bayes factor as varing STRENGTH of two cues
% 3. The decision boundary between integration and segregation models with 
%    input spike counts.      
% 4. Plot the reconstructed likelihood
% 5. Verify the approximation of marginal model evidence is good
% 6. Decode the posterior and Bayes factor from congruent and opposite
%    neurons, and compare with theoretical predictions.
% 7. Plot the tuning of an example congruent and opposite neurons

switch flagDemo
    case 1
        %         demoBayesFactor_Angle;
        demoBFx_Angle;
    case 2
        demoBayesFactor_Lambda;
    case 3
        DecisionBoundary_Lambda;
    case 4
        demoBestFitLikelihood;
    case 5
        cmprMdlEvidenceMargin;
    case 6
        decodeBayesFactor;
    case 7 
        plotCONeuronTuning;
end


