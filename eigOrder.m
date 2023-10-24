function [eee_worst_s, sortIdx] = eigOrder(zeta, frequencies, L, diffeigs, delt, useavg, v2)

    ACs = zeros(diffeigs, L + 1);
    goalACs = zeros(diffeigs, L + 1);
    lags = 0:L;
    compts = delt*lags;
    
    for j = 1:diffeigs
        cf = xcorr(zeta(:,j), L, "normalized");
        ACs(j,:) = cf(L + 1:end);
        goalACs(j,:) = exp(1j*frequencies(j)*compts);
    end

    eee_ts = real(1 - abs(goalACs .* ACs));

    if useavg
        usemm = ceil(0.3 / delt);
        if v2
            testf = double(abs(frequencies) >= 200);
            size(testf)
            eee_worst = mean(eee_ts, 2) + testf;
            size(eee_worst)
        else
            touse = movmean(eee_ts, usemm, 2);
            testf = double(abs(frequencies) >= 200);
            % size(testf)
            % size(touse)
            eee_worst = max(touse, [], 2) + testf;
            % size(eee_worst)
        end
        
    else % eee_worst = max(eee_ts, [], 2) + double(abs(diag(frequencies)) >= 150);
        touse = eee_ts ;
        eee_worst = max(eee_ts , [], 2);
        % size(eee_worst)
    end

    % clear eee_ts, goalACs, ACs;

    [eee_worst_s, sortIdx] = sort(eee_worst,'ascend');
    
end