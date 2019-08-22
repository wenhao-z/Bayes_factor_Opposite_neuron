function probGen = DirGenMdl(x, rho, s, R, parsMdl)
% The generative model of the heading direction x of each cue given stimulus parameters (s,R)
% Wen-Hao Zhang, April 22, 2019

% Lambda = 0: parsMdl.LR + 4*sqrt(parsMdl.LR);

LambdaAvg = parsMdl.beta*parsMdl.LR;
Lambda = 0: LambdaAvg + 3*sqrt(LambdaAvg);
Kappa = parsMdl.TunKappa * rho * Lambda;

if size(s) ~= size(R)
    error('The size of s and R should be the same!');
end

if size(x) == size(s)
    % The function only calculates the value of Bayes factor on input (a value)
    szX = size(x);
    
    x = x(:);
    s = s(:);
    R = R(:);
    
    Lambda = shiftdim(Lambda(:), -1); % [1, nLambda] Array
    Kappa = shiftdim(Kappa(:), -1);
    dimAvg = 2;
    szProbGen = szX;
else
    % Calculate the whole Bayes factor on new inputs (a function)
    Lambda = shiftdim(Lambda(:), -2);
    Kappa = shiftdim(Kappa(:), -1);
    
    x = x(:);
    s = shiftdim(s(:), -1); % [1, n_s] array
    R = shiftdim(R(:), -1); % [1, n_R] array
    
    Lambda = shiftdim(Lambda(:), -2); % [1, 1, nLambda] Array
    Kappa = shiftdim(Kappa(:), -2);
    dimAvg = 3;
    szProbGen = [numel(x), numel(s)];
end

% Calculate the likelihood
LH_Lambda = poisspdf(Lambda, R * parsMdl.beta);
if sum(isinf(LH_Lambda))
    error('Numerical error in Poisson pdf')
end

LH_x = vonMisesPdf(x, s, Kappa);
probGen = sum(LH_x .* LH_Lambda, dimAvg);

probGen = reshape(probGen, szProbGen);

