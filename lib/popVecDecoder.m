function [x, Kappa] = popVecDecoder(nSpk, parsMdl)
% Population vector decoding of population activities
% Wen-Hao Zhang, April 22, 2019

% nSpk: spike count of population responses
% PrefStim: the preferred stimulus of population of neurons

Spk_CompRep = sum(nSpk .* exp(1j * parsMdl.PrefStim), 1);

x = angle(Spk_CompRep); % decoded direction
Kappa = abs(Spk_CompRep) * parsMdl.TunKappa; % The concentration of the decoded distribution

x = x(:);
Kappa = Kappa(:);