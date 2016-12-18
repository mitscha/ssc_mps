% Script to reproduce Tables 1 and 2 in "Noisy subspace clustering via 
% matching pursuits" by Michael Tschannen and Helmut Boelcskei

clear;

addpath ./TSC_OMP
addpath ./greedysc_codes
addpath ./SSC_ADMM_v1.1

datapath = './SSC_ADMM_v1.1/YaleBCrop025.mat';

s = RandStream('mcg16807','Seed',100);
RandStream.setGlobalStream(s);

Ls = [2 3 5 8 10];
nsubj = 38;
nimg = 64;
nrep = 100;
m = 2016;

% SSC-OMP/MP parameters
smax = 5;
tau = 0;

% SSC parameters
affine = false;
alpha = 20;

% TSC parameters
q = 3;

% NSN parameters
K = 20;
kmax = 20;
epsilon = 0.0001;

ces = zeros(4,length(Ls),nrep);
tcomps = ces;

load(datapath,'Y')

% Benchmark algorithms
for l = 1:length(Ls)
    for n = 1:nrep
        try
        idx = randsample(nsubj,Ls(l));
        
        X = [];
        gt = [];
        for i = 0:(Ls(l)-1)
            X = [X squeeze(Y(:,:,idx(i+1)))];
            gt = [gt (i+1)*ones(1,size(Y,2))]; 
        end
        X = normc(X);
        
        % SSC-OMP
        tic
        [labelsomp,~] = OMPSSC(X,tau,smax,true,Ls(l));
        tcomps(1,l,n) = toc;
        
        % SSC-MP
        tic
        [labelsmp,~] = MPSSC(X,tau,smax,smax,true,Ls(l));
        tcomps(2,l,n) = toc;
        
        % SSC
        tic
        Cssc = admmOutlier_mat_func(X,affine,alpha);
        Zssc = abs(Cssc(1:size(Cssc,2),:));
        labelsssc = SpectClust(Zssc + Zssc',Ls(l));
        tcomps(3,l,n) = toc;
        
        % TSC
        tic
        Atsc = TSC(X,q);
        labelstsc = SpectClust(Atsc,Ls(l));
        tcomps(4,l,n) = toc;
        
        % NSN
        tic
        Znsn = NSN(X,K,kmax,epsilon);
        labelsnsn = SpectClust(Znsn+Znsn',Ls(l));
        tcomps(5,l,n) = toc;
        
        ces(1,l,n) = computece(labelsomp,gt);
        ces(2,l,n) = computece(labelsmp,gt);
        ces(3,l,n) = computece(labelsssc,gt);
        ces(4,l,n) = computece(labelstsc,gt);
        ces(5,l,n) = computece(labelsnsn,gt);
        catch e
            fprintf('%s',e.message);
            fprintf('l=%i n=%i idx=\n',l,n);
            disp(idx);
            exit
        end
        
    end
end

cesavg = squeeze(mean(ces,3));
tcompsavg = squeeze(mean(tcomps,3));
cesstd = squeeze(std(ces,0,3));
tcompsstd = squeeze(std(tcomps,0,3));

fprintf('Average clustering errors (SSC-OMP | SSC-MP | SSC | TSC | NSN) x (L=2 | L=3 | L=5 | L=8 | L=10):\n')
disp(cesavg)
fprintf('Clustering error standard deviations:\n')
disp(cesstd)

fprintf('Average running times:\n')
disp(tcompsavg)
fprintf('Running time standard deviations:\n')
disp(tcompsstd)

save('faceclustering.mat','cesavg','tcompsavg','cesstd','tcompsstd')

