% Parameter file of the simulation on neural coding level to calculate
% Bayes factor in multisensory cue processing

% Wen-Hao Zhang, Apr 9, 2019
parsMdl.Ls = 2*pi; % The width of moving direction space
parsMdl.LR = 100; % The width of peak firing rate space

% parsMdl.Stim = [0, 60]; % unit: deg
% parsMdl.R = [30, 30]; % Peak of tuning curves

parsMdl.x = [-60; 0]*pi/180; % Unit: deg
parsMdl.Lambda = [30; 30]; % Spike count of feedforward inputs. Unit: Hz

parsMdl.TunKappa = 3; % The concentration of tuning
parsMdl.N = 180; % The number of neurons

PrefStim = linspace(-pi,pi, parsMdl.N+1)'; 
PrefStim(1) = [];
parsMdl.PrefStim = PrefStim;

clear PrefStim

% --------------------------------------------
% Calculate dependent parameters of the model 
parsMdl = getDependentPars(parsMdl);
