% Plot the tuning of example congruent and opposite neurons
% Wen-Hao Zhang, April 26, 2019


% Load model parameters
% parsProbMdl;

% Set model parameters
parsMdl.Stim = -90*[1,1]*pi/180; % unit: rad
parsMdl.R = [30, 25]; % Peak of tuning curves

parsMdl.R = parsMdl.R / exp(parsMdl.TunKappa);

% Generate feedforward inputs from two sensory modalities (Poisson spikes)

TuningFF = @(s, R) R*exp(parsMdl.TunKappa * cos(s - parsMdl.PrefStim));

RateCongCell = [...
    TuningFF(parsMdl.Stim(1), parsMdl.R(1)), ...
    TuningFF(parsMdl.Stim(2), parsMdl.R(2))];

RateOppoCell = RateCongCell;
RateOppoCell(:,2) = [RateOppoCell(end/2+1:end,2); RateOppoCell(1:end/2,2)];
RateOppoCell = RateOppoCell/2; % This might be a problem of dividing spike count by 2


%% 
for iter = 1:2
   hAxe(iter) = subplot(1,2,iter);
end
plot(hAxe(1), parsMdl.PrefStim*180/pi, RateCongCell)
plot(hAxe(2), parsMdl.PrefStim*180/pi, RateOppoCell)

axis(hAxe, 'square')
set(hAxe, 'xlim', [-180, 180], 'xtick', -180:90:180)
set(hAxe(1), 'ylim', [0, 35])
set(hAxe(2), 'ylim', [0, 20])

xlabel(hAxe(1), 'Cue direction (\circ)')
ylabel(hAxe(1), 'Firing rate (Hz)')
