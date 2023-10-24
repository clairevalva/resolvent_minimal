function U = koopmanOperatorCircular( q, phi, mu, nPar )
    % KOOPMANOPERATOR Compute an approximation of the Koopman operator in a basis 
    % of observables.
    %
    % Input arguments:
    %
    % q:    Scalar or vector of integers containing the shifts to 
    %       apply. (Now can accept non-negative integers)
    % phi:  Array of size [ nS nL ] where nS is the number of samples and nL the 
    %       number of basis functions.
    %
    % mu:   Inner product weight array of size [ nS, 1 ]. The phi's 
    %       are assumed to be orthonormal with respect to the mu inner product. 
    %
    % nPar: Number of workers for parallel for loop calculation. Calculation 
    %       reverts to serial for loop if nPar is set to 0 or is unspecified. 
    %
    % Output arguments:
    %
    % U:    If q is scalar, U is an [ nL, nL ] sized matrix representing the 
    %       Koopman operator in the phi basis. If q is a vector, U is an array of 
    %       size [ nL nL nQ ], where nQ is the length of q, such that U( :, :, i ) 
    %       contains the Koopman operator matrix for shift q( i ).
    %
    % Modified 2022/02/28.
    %
    %
    % test the negative directions 
    % Upos = koopmanOperator(1, phi_c, lam_c, 0 );
    % Uneg = koopmanOperator(-1, phi_c, lam_c, 0 );
    % U0 = koopmanOperator(0, phi_c, lam_c, 0 );
    % do get that Uneg*Upos \appox U0 = I
    
    % Return U as a matrix for a single shift q
    nS = length(mu);
    if isscalar( q )
        if q < 0
            % backward/negative shift operator
            % [mytest, testStart]
            b = abs(q);
            
            ind1 = circshift(1:nS, -1*b);
            U = phi( ind1, : )' * ( phi( : , : ) .* mu( ind1 ) );
        else
            ind1 = circshift(1:nS, -1*q);
            % forward time, original operation
            U = phi( : , : )' * ( phi(ind1, : ) .* mu( : ) );
        end
        return
    end
    
    % Recursive call for a vector of shifts
    if nargin <= 3 
        nPar = 0
    end
    
    nL = size( phi, 2 );
    nQ = numel( q );
    U  = zeros( nL, nL, nQ );
    
    parfor( iQ = 1 : nQ, nPar )
        U( :, :, iQ ) = koopmanOperator( q( iQ ), phi, mu );
    end
    