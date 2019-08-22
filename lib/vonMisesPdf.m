function pdf = vonMisesPdf(x, mu, kappa)
% The probability density function of a von Mises distribution
% Wen-Hao Zhang, Apr 9, 2019
% Unievrsity of Pittsburgh

% When kappa > 700, exp(kappa) will be Inf due to numerical error.
% For large kappa, I approximate the pdf by using Gaussian distribution
Idx = (kappa < 700);

% Note that x and mu should be of unit rad!
pdf = exp(kappa(Idx) .* cos(x-mu));
pdf = pdf./ (2*pi*besseli(0, kappa(Idx)));

% When kappa is too large, the pdf of von Mises distribution will be NaN
% due to numerical error;
if sum(~Idx) > 0
    dx = angle(exp(1i * (x-mu)));
    pdf2 = exp(-dx.^2 .* kappa(~Idx)/2);
    pdf2 = pdf2 .* sqrt(kappa(~Idx)) / sqrt(2*pi);
    pdf = cat(ndims(kappa), pdf, pdf2);
end