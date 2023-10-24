function Xpos = filterfreq(X)
    % filters out the positive frequencies from each X(:,i) row
    N = size(X,1);
    rows = size(X,2);
    f = fftfreq(N,1);
    zmat = (f < 0);

    Xpos = 1j .* zeros(size(X));

    for r = 1:rows
        ffted = fft(X(:,r));
        ffted(zmat) = 0;
        Xpos(:,r) = ifft(ffted);
    end
end