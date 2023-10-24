% last edited by CV, October 2023

%%%%
%  compute resolvent for l63, includes lines to compute NLSA model if necessary, recommended to follow with pltRosslerTiles for plots of results
%
%%%%% 
% relevant resulting arrays:
% zeta : approximate eigenfunctions
% frequencies : corresponding frequencies 
% x  =  getSrcData( model ); : model data
% plottimes = (1:size(toplot,1))*delt; : time array corresponding to data
% 
%%%%

experiment = '64k_dt0.01_nEL0'; % experiment name
computeNLSA = false; % recompute NLSA eigenfunctions?
computeResolvent = true; % recompute resolvent?

% get NLSA model, compute if necessary
if computeNLSA
    demoKoopmanRKHS
end
[ model, In, Out ] = demoKoopmanForecastRKHS_nlsaModel( experiment ); 

resolventmeth = "circshift" % Koopman approximation, using circular shift
savestart = experiment + "_" + resolventmeth % naming conventions
qend = 50; % how long to numerically integrate for (T_\ell in the paper), 50 time units is more than enough
delt = In.dt; % should be the same sampling interval used in the model
z0 = 1; 
tt = 2e-6; % value of tau for G, original was 5
batchnum = 20; % batching parameter (use if diffeigs is big and dt is small)
diffeigs = 2000; % number of NLSA eigenfunctions to use
% batchnum = 20;
M = floor(diffeigs / 6); % rank of final operator
pltTiles = [2, 4, 30] % recommended eigenfunctions to plot to replicate results

% compute and construct resolvent
getPolarOpSqrt

% make folder and name for figures
namecon =  string(savestart) +  ...
 "_z=" + string(z0) + "_q=" + string(qend) + ...
 "_tau=" + string(tt) + "_eigs=" + string(diffeigs) + ... 
 "_M=" + string(M);

savestart = savestart
if ~exist('figs/' + savestart, 'dir')
    mkdir('figs/' + savestart)
end

% sort eigenfunctions based on \eps_{T_c} criteria with T_c = nL * delt
nL = 200;
[eee_s, sIdx] = eigOrder(zeta, frequencies, nL, diffeigs, delt, true, true);
gamma = gamma(sIdx);
frequencies = frequencies(sIdx);
zeta = zeta(:,sIdx);
