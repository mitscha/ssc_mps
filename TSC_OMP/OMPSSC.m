% Implementation of SSC-OMP with fixed number of iterations
% 
% X: m x N matrix of N data points
% smax: number of OMP iterations
% dosc: perform spectral clustering if true
% L: number of clusters
% labels: labels of the data points
% C: coefficient matrix

function [labels, C] = OMPSSC(X,tau,smax,dosc,L)
    N = size(X,2);
    m = size(X,1);
    C = zeros(N);

    if nargin < 4
        dosc = true;
    end
    if nargin < 5
        L = 0;
    end
    
    for n = 1:N
       r = X(:,n);
       Xc = X;
       Xc(:,n) = nan;
       Xs = zeros(m,smax);
       lambda = zeros(smax,1);
       
       % OMP iterations
       for s = 1:smax
           [~,idx] = nanmax(abs(Xc'*r));
           if s == 1;
               xsp = Xc(:,idx)/sqrt(Xc(:,idx)'*Xc(:,idx));
           else
               xsp = Xc(:,idx) - Xs(:,1:(s-1))*(Xs(:,1:(s-1))'*Xc(:,idx));
               xsp = xsp/sqrt(xsp'*xsp);
           end
           r = r - xsp*(xsp'*r);
           Xs(:,s) = xsp;
           Xc(:,idx) = nan;
           lambda(s) = idx;
           if sqrt(r'*r) <= tau
                break
           end
       end

       C(lambda(lambda > 0.0),n) = X(:,lambda(lambda > 0.0))\X(:,n); %pinv(X(:,lambda(lambda > 0.0)))*X(:,n);
    end
    
    if dosc
        labels = SpectClust(abs(C) + abs(C)',L);
    else
        labels = NaN;    
    end
        
end