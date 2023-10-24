% last edited by CV, October 2023

%%%%
%  compute resolvent for Rossler system, includes lines to compute NLSA model if necessary, recommended to follow with pltRosslerTiles for plots of results
%
%%%%% 
% relevant resulting arrays:
% zeta : approximate eigenfunctions
% frequencies : corresponding frequencies 
% x  =  getSrcData( model ); : model data
% plottimes = (1:size(toplot,1))*delt; : time array corresponding to data
% 
%%%%

% experiment = "64k_dt0.04_nEL0"; 
experiment = "6.4k_dt0.04_nEL0"; 
computeNLSA = false; % recompute NLSA eigenfunctions?
computeResolvent = true; % recompute resolvent?

% get NLSA model, compute if necessary
if computeNLSA
    demoKoopmanRKHS
end
[model, In] = demoKoopmanRKHS_nlsaModel('rossler_6.4k');

% set resolvent parameters
switch experiment
    case '6.4k_dt0.04_nEL0'
        diffeigs = 201; % number of NLSA eigenfunctions to use
        M = floor(diffeigs / 3); % rank of final operator
        batchnum = 1; % batching (use if diffeigs is very large)
        pltTiles = [2, 10, 30] % recommended eigenfunctions to plot to replicate results
    case '64k_dt0.04_nEL0'
        diffeigs = 2001; % number of NLSA eigenfunctions to use
        M = floor(diffeigs / 3); % rank of final operator
        batchnum = 1; % batching (use if diffeigs is very large)
        pltTiles = [4, 2, 18] % recommended eigenfunctions to plot to replicate results
    
    otherwise
        error('Invalid experiment')
end

% the following parameters are shared between experiments
resolventmeth = "circshift" % Koopman approximation, using circular shift
savestart = experiment + "_" + resolventmeth
qend = 50; % how long to numerically integrate for (T_\ell in the paper), 50 time units is more than enough
delt     = In.dt % dt of data
z0 = 1; % z to compute resolvent at
tt = 1e-4; % value of \tau, this one seems to work well

% compute and construct resolvent
getPolarOpSqrt

% make folder and name for figures
namecon = string(savestart) +  ...
"_z=" + string(z0) + "_q=" + string(qend) + ...
"_tau=" + string(tt) + "_eigs=" + string(diffeigs) + ... 
"_M=" + string(M);

% create figure folder if needed
if ~exist('figs/' + savestart, 'dir')
    mkdir('figs/' + savestart)
end

% sort eigenfunctions based on \eps_{T_c} criteria with T_c = nL * delt
nL = 1000;
[eee_s, sIdx] = eigOrder(zeta, frequencies, nL, diffeigs, delt, true, true);

gamma = gamma(sIdx);
frequencies = frequencies(sIdx);
zeta = zeta(:,sIdx);