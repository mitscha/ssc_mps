function labels = SpectClust(Z,L)
    N = size(Z,1);

    Dsum = diag(1./sqrt(sum(Z)+eps));
    LN = eye(N) - Dsum*Z*Dsum;
    [~,S,V] = svd(LN);
    
    if nargin < 2
        singvals = diag(S);
        [~,minidx] = min(diff(singvals(1:(end-1))));
        L = N - minidx;
    end
    
    VL = V(:,N-L+1:N);
    VL = normr(VL);
    [labels,~] = kmeans(VL,L,'maxiter',1000,'replicates',50,'EmptyAction','singleton');
end