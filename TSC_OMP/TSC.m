% 
% This function implements the TSC algorithm from the paper 
% ``Robust subspace clustering via thresholding'' by Reinhard Heckel and Helmut B??lcskei
% Reinhard Heckel, 2013
%
% Modified version by Michael Tschannen 2015: External spectral clustering step
%
% X: m x N matrix of N data points
% q: input parameter of TSC
% L: number of clusters, optional. If not provided, L is estimated via the eigengap heuristic%
% labels: labels of the data points
% Z: adjacency matrix 
% nL: estimated number of clusters

function Z = TSC(X,q)

% normalize the data points
X = normc(X);

[m,N] = size(X);

Z = zeros(N,N);

for i=1:N
    corvec = abs(X'*X(:,i));
    corvec(i) = 0; % so TSC will not select it
    [el,order] = sort(corvec, 'descend');
    Z(i, order(1:q) ) = exp(-2*acos(el(1:q))); % better than squared arcsin
end

Z = Z + Z';

end


