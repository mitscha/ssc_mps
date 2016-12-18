% Function to reproduce Figure 4 in "Noisy subspace clustering via 
% matching pursuits" by Michael Tschannen and Helmut Boelcskei

function itersensitivity(mode,sigman)
    datapath = './SSC_ADMM_v1.1/YaleBCrop025.mat';
    addpath ./TSC_OMP

    s = RandStream('mcg16807','Seed',1);
    RandStream.setGlobalStream(s);

    switch(mode)
        % Face clustering params
        case 1
            L = 3;
            nsubj = 38;
            filename = strcat('itersens_faces');
            smaxs = 5:5:100;
            nrep = 100;
            load(datapath)
            
        % Synthetic data params
        case 2
            m = 80; 
            dint = 3;
            dl = 15*ones(1,3);
            rho = 4;
            nl = rho*dl;
            L = 3;
            if nargin < 2
                sigman = 0.1;
            end
            smaxs = 4:4:60;
            nrep = 50;
            filename = strcat('itersens_synth_',num2str(sigman,'%0.2f'));
    end

    tau = 0;

    ces = zeros(3,length(smaxs),nrep);
    tcomps = ces;
    tpr = ces;
    fpr = ces;
    tnorm = ces;
    fnorm = ces;

    for n = 1:nrep
        for s = 1:length(smaxs)

            X = [];
            gt = [];
            switch(mode)
                % Generate synthetic data set
                case 1
                    idx = randsample(nsubj,L);
                    for l = 1:L
                        X = [X squeeze(Y(:,:,idx(l)))];
                        nl(l) = size(Y(:,:,idx(l)),2);
                        gt = [gt l*ones(1,nl(l))];
                    end
                % Generate face clustering data set
                case 2
                    Uint = orth(randn(m,dint));
                    for l = 1:L
                        U = [Uint orth(randn(m,dl(l)-dint))];
                        X = [X U*normc(randn(dl(l),nl(l)))];
                        gt = [gt l*ones(1,nl(l))];
                    end

                    if sigman > 0.0
                        X = X + sigman/sqrt(m)*randn(size(X));
                    end

            end

            dataidx = [0 cumsum(nl)];

            X = normc(X);
            
            tic
            [labelsomp,Comp] = OMPSSC(X,tau,smaxs(s),true,L);
            tcomps(1,s,n) = toc;
            tic
            [labelsmp,Cmp] = MPSSC(X,tau,smaxs(s),smaxs(s),true,L);
            tcomps(2,s,n) = toc;
            tic
            [labelsmppmax,Cmppmax] = MPSSC(X,tau,inf,smaxs(s),true,L);
            tcomps(3,s,n) = toc;

            ces(1,s,n) = computece(labelsomp,gt);
            ces(2,s,n) = computece(labelsmp,gt);
            ces(3,s,n) = computece(labelsmppmax,gt);

            for l = 1:L
                for j = 1:3
                    switch(j)
                        case 1 
                            currblock = abs(Comp(:,(dataidx(l)+1):dataidx(l+1)));
                        case 2
                            currblock = abs(Cmp(:,(dataidx(l)+1):dataidx(l+1)));
                        case 3
                            currblock = abs(Cmppmax(:,(dataidx(l)+1):dataidx(l+1)));
                    end
                    % Compute performance metrics
                    tp = sum(sum(currblock((dataidx(l)+1):dataidx(l+1),:) > 0.0));
                    fp = sum(currblock(:) > 0.0) - tp;
                    tpr(j,s,n) = tpr(j,s,n) + tp/(nl(l)*(nl(l)-1)*L);
                    fpr(j,s,n) = fpr(j,s,n) + fp/((size(currblock,1)-nl(l))*nl(l)*L);

                    tn = sum(sum(currblock((dataidx(l)+1):dataidx(l+1),:)));
                    fn = sum(currblock(:)) - tn;
                    tnorm(j,s,n) = tnorm(j,s,n) + tn/(nl(l)*L);
                    fnorm(j,s,n) = fnorm(j,s,n) + fn/(nl(l)*L);

                end
            end
        end
    end

    % Average performance metrics
    cesavg = squeeze(mean(ces,3));
    tcompsavg = squeeze(mean(tcomps,3));
    tpravg = squeeze(mean(tpr,3));
    fpravg = squeeze(mean(fpr,3));
    tnormavg = squeeze(mean(tnorm,3));
    fnormavg = squeeze(mean(fnorm,3));

    save(strcat(filename,'.mat'),'smaxs','ces','tcomps','tpr','fpr','tnorm','fnorm');
    
    dlmwrite(strcat('CEs-',filename,'.dat'),...
        [smaxs' cesavg'],'delimiter',' ');
    dlmwrite(strcat('CEstd-',filename,'.dat'),...
        [smaxs' squeeze(std(ces,0,3))'],'delimiter',' ');
    dlmwrite(strcat('TPFPs-',filename,'.dat'),...
        [smaxs' tpravg' fpravg'],'delimiter',' ');
    dlmwrite(strcat('TF1norm-',filename,'.dat'),...
        [smaxs' tnormavg' fnormavg'],'delimiter',' ');
    dlmwrite(strcat('Tcomps-',filename,'.dat'),...
        [smaxs' squeeze(tcompsavg)'],'delimiter',' ');
end



