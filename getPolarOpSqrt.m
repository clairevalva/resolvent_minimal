% last edited by CV, October 2023

%%%%
% computes resolvent with parameters from the ``driver'' file
%%%%

sn = string(savestart) + "_z=" + string(z0) + "_q=" + string(qend) + ".mat";
sn = "resdat/" + resolventmeth + "_" + sn

if ~exist('resdat', 'dir')
    mkdir('resdat')
end

if ~computeResolvent
    loaded = load(sn);
    R_z = loaded.R_z;
    G = Gtau(model, tt, false);
    
else
    disp( 'Computing positive resolvent...' ); t = tic; 
    [R_z, G] = posApproxResolvent(model, z0, delt, qend, batchnum, tt, resolventmeth, diffeigs);
    toc( t );
    save(sn,'R_z');
end

diffOp = getDiffusionOperator( model);
[phi_orig, mu] = getEigenfunctions( diffOp );


R_z = R_z(1:diffeigs,1:diffeigs);
G = G(1:diffeigs, 1:diffeigs);
center = 1/(2*z0);

[U, P] = polardecomp(R_z);
[P_sqrt, residual] = sqrtm(P);
Pt = P_sqrt * G * P_sqrt;

[c,gamma_polar] = eig(Pt);
gamma_polar = real(gamma_polar);

[zeta, gamma_polar, c] = sortEig(c, gamma_polar, diffOp);
gamma = ainv(gamma_polar, z0); 
frequencies = imag(1./gamma);


phi_orig = phi_orig(:, 1:diffeigs);
s = abs(gamma_polar);
[s, idxs] = sort(s, 'descend');
psi = zeta(:, idxs);
[N, L] = size(psi);

Psi = phi_orig' * psi / N;
PsiBar = conj(Psi);

Splus = Psi(:, 2 : M) * diag(s(2 : M)) * Psi(:, 2 : M)' ; % N \times N = number of diffusion eigenfunctions
Sminus = - conj(Psi(:, 2 : M)) * diag(s(2 : M)) * conj(Psi(:, 2 : M))' ;
S = Splus + Sminus;

[c, s_new] = eig(S, 'vector');
psi_new = phi_orig * c;

ispos = s_new >= 0;
rho_new = zeros(size(s_new));
rho_new(ispos) = ainv(z0, s_new(ispos));
rho_new(~ispos) = conj(ainv(z0, -s_new(~ispos)));
omega_new = rinv(z0, rho_new);

% final wanted frequencies and eigenfunctions (called zeta)
frequencies = real(omega_new);
zeta = psi_new;


function l = ainv(z, x)
    l = x .* exp(i * (pi / 2 + asin(x * z)));
end

function omega = rinv(z, rho)
    omega = (rho .^ -1 + z) / i;
end



