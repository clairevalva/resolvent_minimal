function G = Gtau(model, tau, alteta)
    %   G:  G^(1/2) is the eigenvalues representing Gtau in the phi basis. 
    %   phi_i = getDiffusionEigenfunctions(model);

    lam_i = getDiffusionEigenvalues(model);
    if alteta
        % \eta_j = \log(\Lambda_j)/ \log(\Lambda_1)
        eta_i = log(lam_i) ./ log(lam_i(2));
        eta_i(1) = 0;
    else
        eta_i = lam_i.^(-1) - 1;
        eta_i = eta_i ./ eta_i(2);
    end
    lam_it = exp(-1*eta_i*tau);
    G = diag((lam_it).^(1/2));
end