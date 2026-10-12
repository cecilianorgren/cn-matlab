function e = l2err_partition(w, mu, Sig, labels)
    % w: 1xK weights, mu: Kxd means, Sig: dxdxK covariances,
    % labels: 1xK group index of each component, e.g. [1 1 2 3 3]
    ids = unique(labels);  J = numel(ids);  d = size(mu,2);
    W = zeros(1,J);  muG = zeros(J,d);  SigG = zeros(d,d,J);
    for j = 1:J
        idx = find(labels == ids(j));
        W(j)     = sum(w(idx));
        muG(j,:) = (w(idx)*mu(idx,:)) / W(j);
        S = zeros(d);
        for k = idx
            dm = mu(k,:) - muG(j,:);
            S  = S + w(k)*(Sig(:,:,k) + dm'*dm) / W(j);   % law of total variance
        end
        SigG(:,:,j) = S;
    end
    ff = overlap(w, mu, Sig, w, mu, Sig);
    fg = overlap(w, mu, Sig, W, muG, SigG);
    gg = overlap(W, muG, SigG, W, muG, SigG);
    e  = sqrt(max(ff - 2*fg + gg, 0) / ff);
end

function s = overlap(w1, mu1, S1, w2, mu2, S2)
    % sum_{a,b} w1(a) w2(b) N(mu1(a); mu2(b), S1(:,:,a) + S2(:,:,b))
    s = 0;
    for a = 1:numel(w1)
        for b = 1:numel(w2)
            s = s + w1(a)*w2(b)*mvnpdf(mu1(a,:), mu2(b,:), S1(:,:,a) + S2(:,:,b));
        end
    end
end