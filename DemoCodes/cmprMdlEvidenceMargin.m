% Compare the marginalized Model evidence over moving direction and
% theoretical approximation

% Wen-Hao Zhang,
% Apr 17,2019

% Load model parameters
% parsProbMdl;

% Set the range of reconstructed cues
RArray = 0:100; % Note that this R is the firing rate of the spike count of all neurons, 
%                 not the tuning height as specified in the papers.
R_demo = 30;
x_recon = linspace(-pi/2, pi/2, parsMdl.N+1);

[x_recon,Lambda_recon] = ndgrid(x_recon, 0: RArray(end) + 4*sqrt(RArray(end)));
Kappa_recon = parsMdl.TunKappa * parsMdl.rho * Lambda_recon;

Evid_Array = zeros(1, length(RArray));
Evid2d_demo = zeros(size(x_recon));

for iterPar = 1: length(RArray)
%     fprintf('Progress:%d\n', iterPar)
    Evid = vonMisesPdf(x_recon, 0, Kappa_recon) ...
        .* poisspdf(Lambda_recon, RArray(iterPar));
    if RArray(iterPar) == R_demo
       Evid2d_demo = Evid; 
    end
    Evid = sum(Evid,2);
    
    Evid = Evid/sum(Evid);
    Evid_Array(iterPar) = sum(Evid(:) .* exp(1j * x_recon(:,1)));
end

%%
% Estimate the position and concentration of model evidence under each
% model parameter

hAxe(1) = subplot(2,2,1); hold on;
hAxe(2) = subplot(2,2,3); hold on;
hAxe(3) = subplot(2,2,[2,4]); hold on;
% for iter = 1: 3
%    hAxe(iter) = subplot(1,3,iter);
%    hold on
% end
nLines_Contour = 10;

contourf(hAxe(1), x_recon*180/pi, Lambda_recon, Evid2d_demo, nLines_Contour)
set(hAxe(1), 'ylim', [10, 50], 'xlim', 30*[-1,1])
xlabel(hAxe(1), 'Cue direction (\circ)')
ylabel(hAxe(1), 'Spike count')
colorbar(hAxe(1))

EvidNum = sum(Evid2d_demo, 2);
EvidNum = EvidNum/sum(EvidNum);
EvidTheory = vonMisesPdf(x_recon, 0, parsMdl.TunKappa * parsMdl.rho * R_demo);
EvidTheory = EvidTheory/sum(EvidTheory);

plot(hAxe(2), x_recon(:,1)*180/pi, EvidNum, '-o');
plot(hAxe(2), x_recon(:,1)*180/pi, EvidTheory);
set(hAxe(2), 'xlim', 30*[-1, 1], 'ylim', [0, 0.07])
set(hAxe(1:2), 'xtick', -30:15:30)
xlabel(hAxe(2), 'Cue direction (\circ)')
ylabel(hAxe(2), 'Probability distribution')
legend(hAxe(2), 'Num.', 'Theory')

plot(hAxe(3), mrl2Kappa(abs(Evid_Array)), RArray*parsMdl.TunKappa * parsMdl.rho, 'o')
KappaLim = RArray*parsMdl.TunKappa * parsMdl.rho;
plot(hAxe(3), [0, KappaLim], [0, KappaLim], '--k')
xlabel(hAxe(3), {'Concentration (Num.)'})
ylabel(hAxe(3), {'Concentration (Theo. approx.)'})
axis(hAxe(3), 'square')
% axis(hAxe, 'square')

