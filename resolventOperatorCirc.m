function R = resolventOperatorCirc(z, delt, q, phi2, mu,batchnum)
    U0 = koopmanOperatorCircular(1, phi2, mu, 0 );
    
    numadd = size(q, 2);
    ts = q*delt;
    tend = ts(end);

    % Uts = zeros(numadd, size(U0));
    tosum = zeros(size(U0));
    lens = ceil(numadd / batchnum);
    for j = 1:batchnum
        if j == batchnum
            arrstart = (j-1)*lens + 1;
            arrend = numadd;
        else
            arrstart = (j-1)*lens + 1;
            arrend = j*lens;
        end  

        bnds = arrstart:arrend;
        ttt = ts(bnds);
        qqq = q(bnds);

        lenbatch = size(qqq,2)
        Uts = squeeze(zeros(size(U0,2),size(U0,2), lenbatch));
        for k = 1:lenbatch
            kpow = q(k);
            Uts(:,:,k) = U0^kpow;
        end 

        eee = exp(-1*z*ttt);
        
        V = bsxfun(@times,Uts,reshape(eee,1,1,[]));
        V = permute(V, [3,1,2]);
        sumhold = squeeze(simps(ts(arrstart:arrend), V));
        size(V)

        
        tosum = tosum + sumhold;
    end   
    R = tosum;     
    
end