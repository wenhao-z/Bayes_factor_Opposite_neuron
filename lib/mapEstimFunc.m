function ParamMAP = mapEstimFunc(parsMdl)
% MAP estimate of parameter of integration and segregation model.
% Wen-Hao Zhang, Apr 9, 2019

% Ensure that parsMdl.x has unit of rad.

%% Initialize an empty struct and return
if isempty(parsMdl)
    ParamMAP = struct(...
        'StimEstim_Seg', [], ...
        'KappaEstim_Seg', [], ...
        'REstim_Seg', [], ...
        'StimEstim_Int', [], ...
        'KappaEstim_Int', [], ...
        'REstim_Int', [], ...
        'StimDiff', [], ...
        'KappaDiff', []);
    return;
end

%%
% Stimulus estimate (MAP) of segregation model
ParamMAP.StimEstim_Seg  = parsMdl.x; % unit: rad
ParamMAP.KappaEstim_Seg = parsMdl.Kappa;
ParamMAP.REstim_Seg     = parsMdl.Lambda / parsMdl.beta;

% Stimulus estimate (MAP) of integration model
postStim_int = sum(parsMdl.Kappa .* exp(1j* parsMdl.x));
ParamMAP.StimEstim_Int  = angle(postStim_int);
ParamMAP.KappaEstim_Int = abs(postStim_int);
ParamMAP.REstim_Int     = sum(parsMdl.Lambda)/ (2*parsMdl.beta);

% The estimate of Bayes factor of moving direction, which is the ratio
% between likelihood of generating the direction of cue 1 over the
% likelihood of cue 2
BF_Dir = - diff(parsMdl.Kappa .* exp(1j* parsMdl.x));
ParamMAP.StimDiff = angle(BF_Dir);
ParamMAP.KappaDiff = abs(BF_Dir);
