% Implementation of the SSC-MP algorithm
% 
% X: m x N matrix of N data points
% smax: max. number of MP iterations
% pmax: sparsity level
% dosc: perform spectral clustering if true
% L: number of clusters
% labels: labels of the data points
% C: coefficient matrix

function [labels, C] = MPSSC(X,tau,smax,pmax,dosc,L)
    N = size(X,2);
    C = zeros(N);

    if nargin < 4
        pmax = N;
    end
    if nargin < 5
        dosc = true;
    end
    if nargin < 6
        L = 0;
    end
    
    for n = 1:N
       r = X(:,n);
       Xc = X;
       Xc(:,n) = nan;
       c = zeros(N,1);
       
       % MP iterations
       nit = 0;
       while (sqrt(r'*r) > tau) && (nit < smax) && (sum(abs(c) > 0) < pmax)
           [~,idx] = nanmax(abs(Xc'*r));
           cidx = Xc(:,idx)'*r/(Xc(:,idx)'*Xc(:,idx));
           c(idx) = c(idx) + cidx;
           r = r - cidx*Xc(:,idx);
           nit = nit + 1;
       end
       C(:,n) = c;
    end
    
    if dosc

        labels = SpectClust(abs(C)+abs(C)',L);
    else
        labels = NaN;    
    end
end