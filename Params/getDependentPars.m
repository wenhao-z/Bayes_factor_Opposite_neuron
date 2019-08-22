function parsMdl = getDependentPars(parsMdl, nSpk)
% Calculate dependent parameters of the model
% Wen-Hao Zhang, Apr 10, 2019


Rate = exp(parsMdl.TunKappa* (cos(0*pi/180 - parsMdl.PrefStim)-1));
parsMdl.beta = sum(Rate,1); % sum of all neurons' tuning

% Mean resultant length of normalized tunings, indicating tuning width 
if exist('nSpk', 'var')
    rateComplex = nSpk .* exp(1j * parsMdl.PrefStim);
    parsMdl.rho = abs(sum(rateComplex,1))./ sum(nSpk, 1);
    parsMdl.rho = parsMdl.rho(:);
else
    rateComplex = Rate .* exp(1j * parsMdl.PrefStim);
    parsMdl.rho = abs(sum(rateComplex))/ sum(Rate);
end

parsMdl.Kappa = parsMdl.TunKappa .* parsMdl.rho .* parsMdl.Lambda; % Concentration of likelihood of direction