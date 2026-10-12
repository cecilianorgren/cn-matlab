function e = gmm__l2err(w, mu, Sig)
    % w: 1xK weights, mu: Kxd means, Sig: dxdxK covariances.
    % Normalized L2 error between the mixture and its moment-matched single Gaussian.
    K  = numel(w);
    wM = sum(w);
    muM = (w*mu)/wM;                                   % 1xd
    SigM = zeros(size(Sig,1));
    for k = 1:K
        dm = mu(k,:) - muM;
        SigM = SigM + w(k)*(Sig(:,:,k) + dm'*dm)/wM;   % law of total variance
    end
    ff = 0;  fg = 0;
    for a = 1:K
        for b = 1:K
            ff = ff + w(a)*w(b)*mvnpdf(mu(a,:), mu(b,:), Sig(:,:,a)+Sig(:,:,b));
        end
        fg = fg + w(a)*wM*mvnpdf(mu(a,:), muM, Sig(:,:,a)+SigM);
    end
    gg = wM^2*mvnpdf(muM, muM, 2*SigM);
    e  = sqrt(max(ff - 2*fg + gg, 0) / ff);
end