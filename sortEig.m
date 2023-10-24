function [zeta, gamma, c] = sortEig(c, gamma, diffOp)
    gamma         = diag( gamma ).';
    Lambda = getEigenvalues( diffOp );
    M = size(c,1);
    
    Lambda = Lambda(1:M,:);
    [ phi2, mu ] = getEigenfunctions( diffOp );
    
    % Compute Dirichlet energies of eigenvectors
    idxPhi = 1:(size(Lambda, 1));
    Lambda = Lambda( idxPhi );
    Lambda = Lambda( : );
    c = c(1:M, 1:M);
    E = sum( abs( c) .^ 2 ./ Lambda );
    idxE = 1:(size(Lambda, 1)); % should be the same as idxPhi
    nEig = numel( idxE );

    % Sort results in order of increasing Dirichlet energy
    [ E, idxE ] = sort( E, 'ascend' );
    E = E( 1 : nEig );
    idxE = idxE( 1 : nEig );
    gamma = gamma( idxE );
    c = c( :, idxE );

    phi = phi2( :, idxPhi ); 
    zeta = phi * c;

end

% function [zeta, gamma, c] = sortEig(c, gamma, diffOp)
%     gamma         = diag( gamma ).';
%     Lambda = getEigenvalues( diffOp );
%     M = size(c,1);
    
%     Lambda = Lambda(1:M,:);
%     [ phi2, mu ] = getEigenfunctions( diffOp );
    
%     % Compute Dirichlet energies of eigenvectors
%     idxPhi = 2:(size(Lambda, 1));
%     Lambda = Lambda( idxPhi );
%     Lambda = Lambda( : );
%     c = c(2:M, 2:M);
%     E = sum( abs( c) .^ 2 ./ Lambda );
%     idxE = 2:(size(Lambda, 1)); % should be the same as idxPhi
%     nEig = numel( idxE );

%     % Sort results in order of increasing Dirichlet energy
%     [ E, idxE ] = sort( E, 'ascend' );
%     E = E( 1 : nEig );
%     idxE = idxE( 1 : nEig );
%     gamma = gamma( idxE );
%     c = c( :, idxE );

%     phi = phi2( :, idxPhi ); 
%     zeta = phi * c;

% end