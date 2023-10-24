function [R_z, G] = posApproxResolvent(model, z0, delt, qend, batchnum, tau, method, L)
    q = [0:qend/delt];
    G = Gtau(model, tau, false);

    diffOp = getDiffusionOperator( model);
    [ phi_circ, lambda_circ ] = getEigenfunctions( diffOp ); 
    
    phi_c = filterfreq(phi_circ);
    lam_c = lambda_circ;

    if method == "power"
        % uses power method, left shift, no svd limiting
        R_z_power = resolventOperatorPower(z0, delt, q, phi_c, lam_c, batchnum);
        R_z = R_z_power;
    elseif method == "pos"
        % uses sequences of left shift 
        R_z_plus = resolventOperator(z0, delt, q, phi_c, lam_c, batchnum); 
        R_z = R_z_plus;
    elseif method == "neg"
        % uses sequences of right shift 
        R_z_minus = resolventOperator(z0, delt, -1*q, phi_c, lam_c, batchnum);
        R_z = R_z_minus;
    elseif method == "svdcon"
        % uses power method, left shift, limits svd to 1
        R_z= resolventOperatorPowerConstrain(z0, delt, q, phi_c, lam_c, batchnum);
    elseif method == "circshift"
        % uses power method, circular shift
        R_z = resolventOperatorCirc(z0, delt, q, phi_c, lam_c, batchnum);
    else
        disp(" no valid method chosen, try 'pos'")
    end
end