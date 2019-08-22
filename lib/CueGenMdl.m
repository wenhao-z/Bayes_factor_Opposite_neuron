function probGen = CueGenMdl(x, Lambda, rho, s, R, parsMdl)
% The generative model of each cue(x, Lambda) given stimulus parameters (s,R)
% Wen-Hao Zhang, April 22, 2019


Kappa = parsMdl.TunKappa * rho * Lambda;
probGen = vonMisesPdf(x, s, Kappa) .* poisspdf(Lambda, R * parsMdl.beta);