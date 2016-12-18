% Function to reproduce Figure 2 in "Noisy subspace clustering via 
% matching pursuits" by Michael Tschannen and Helmut Boelcskei

function phasediag()
    addpath ./TSC_OMP
    s = RandStream('mcg16807','Seed',1);
    RandStream.setGlobalStream(s);

    % Clustering problem parameters
    d = 20;
    m = 200;
    L = 3;
    
    ts = 1:1:20;
    rhos = 0.5:0.2:4;
    sigmans = 0:0.1:2;
    
    tfix = 8;
    rhofix = 3;
    sigmafix = 1;

    nrep = 20;

    % Algorithm parameters
    tau = 0;
    smax = 10;

    [ce_tvsrho_omp, ce_tvsrho_mp] = evaluateParams(d,m,L,ts,sigmafix,rhos,tau,smax,nrep);
    [ce_tvssig_omp, ce_tvssig_mp] = evaluateParams(d,m,L,ts,sigmans,rhofix,tau,smax,nrep);
    [ce_rhovssig_omp, ce_rhovssig_mp] = evaluateParams(d,m,L,tfix,sigmans,rhos,tau,smax,nrep);

    saveheatmap(ce_tvsrho_omp,sqrt(ts/d),rhos,'CE-tvsrho_omp.dat')
    saveheatmap(ce_tvsrho_mp,sqrt(ts/d),rhos,'CE-tvsrho_mp.dat')
    saveheatmap(ce_tvssig_omp,sqrt(ts/d),sigmans,'CE-tvssig_omp.dat')
    saveheatmap(ce_tvssig_mp,sqrt(ts/d),sigmans,'CE-tvssig_mp.dat')
    saveheatmap(ce_rhovssig_omp,sigmans,rhos,'CE-rhovssig_omp.dat');
    saveheatmap(ce_rhovssig_mp,sigmans,rhos,'CE-rhovssig_mp.dat');

    save('phasediag.mat','ce_tvsrho_omp','ce_tvsrho_mp','ce_tvssig_omp','ce_tvssig_mp','ce_rhovssig_omp','ce_rhovssig_mp');

end


function [cesomp, cesmp] = evaluateParams(d,m,L,ts,sigmans,rhos,tau,smax,nrep)
    cesomp = zeros(length(ts),length(sigmans),length(rhos));
    cesmp = cesomp;
    for t = 1:length(ts)
        for s = 1:length(sigmans)
            for r = 1:length(rhos)
                for n = 1:nrep
                    X = zeros(m,round(L*d*rhos(r)));
                    gt = zeros(1,round(L*d*rhos(r)));
                    U = orth(randn(m,L*(d-ts(t))+ts(t)));
                    
                    % Generate data set and ground truth
                    for l = 0:(L-1)
                        X(:,int32(l*rhos(r)*d+1):1:int32((l+1)*rhos(r)*d)) = ...
                            [U(:,1:ts(t)) U(:,(ts(t)+l*(d-ts(t))+1):(ts(t)+(l+1)*(d-ts(t))))]*normc(randn(d,int32(rhos(r)*d)));
                        gt(int32(l*d*rhos(r)+1):1:int32((l+1)*d*rhos(r))) = l+1;
                    end

                    if sigmans(s) > 0.0
                        X = X + sigmans(s)/sqrt(m)*randn(size(X));
                    end

                    % Evaluate algorithm performance
                    [labelsomp,~] = OMPSSC(X,tau,smax,true,L);
                    [labelsmp,~] = MPSSC(X,tau,smax,smax,true,L);
                    ceomp = computece(labelsomp,gt);
                    cemp = computece(labelsmp,gt);
                    
                    cesomp(t,s,r) = cesomp(t,s,r) + ceomp;
                    cesmp(t,s,r) = cesmp(t,s,r) + cemp;
                end
            end
        end
    end
    
    % Average clustering error
    cesomp = squeeze(cesomp/nrep);
    cesmp = squeeze(cesmp/nrep);
end

