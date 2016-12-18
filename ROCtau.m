% Script to reproduce Figure 3 in "Noisy subspace clustering via matching 
% pursuits" by Michael Tschannen and Helmut Boelcskei

clear;

addpath ./TSC_OMP

s = RandStream('mcg16807','Seed',1);
RandStream.setGlobalStream(s);

% Clustering problem parameters
m = 300;
dint = 4;
dl = [20 40 60 80];
rho = 4;
nl = rho*dl;
L = length(dl);
dataidx = [0 cumsum(nl)];
sigmans = [0.2 0.5];

% Algorithm parameters
tausomp = 10.^(-linspace(0.1,2,15));
tausmp = 10.^(-linspace(0.1,1,15));

nrep = 20;

truediscomp = zeros(length(tausomp), length(dl), length(sigmans));
truediscmp = truediscomp;
tpromp = truediscomp;
tprmp = truediscomp;
fpromp = truediscomp;
fprmp = truediscomp;


for t = 1:length(tausomp)
    for s = 1:length(sigmans)
        for n = 1:nrep
            
            % Generate data set
            X = zeros(m,dataidx(end));
            Uint = orth(randn(m,dint));
            for l = 1:L
                U = [Uint orth(randn(m,dl(l)-dint))];
                X(:,(dataidx(l)+1):dataidx(l+1)) = U*normc(randn(dl(l),nl(l)));
            end
            
            % Add noise
            if sigmans(s) > 0.0
                X = X + sigmans(s)/sqrt(m)*randn(size(X));
            end
            
            [~,Comp] = OMPSSC(X,tausomp(t),dataidx(end)-1,false,0);
            [~,Cmp] = MPSSC(X,tausmp(t),inf,dataidx(end)-1,false,0);
            
            % Compute performance metrics
            for l = 1:L
                blockomp = abs(Comp(:,(dataidx(l)+1):dataidx(l+1))) > 0;
                tdomp = sum(blockomp((dataidx(l)+1):dataidx(l+1),:));
                % #TP OMP
                truediscomp(t,l,s) = truediscomp(t,l,s) +...
                    sum(tdomp)/nl(l);
                % TPR OMP
                tpromp(t,l,s) = tpromp(t,l,s) + sum(tdomp)/(nl(l)*dl(l));
                % FPR OMP
                fpromp(t,l,s) = fpromp(t,l,s) + sum(sum(blockomp) - tdomp)/(nl(l)*(m-dl(l)));
                
                blockmp = abs(Cmp(:,(dataidx(l)+1):dataidx(l+1))) > 0;
                tdmp = sum(blockmp((dataidx(l)+1):dataidx(l+1),:));
                % #TP MP
                truediscmp(t,l,s) = truediscmp(t,l,s) +...
                    sum(tdmp)/nl(l);
                % TPR MP
                tprmp(t,l,s) = tprmp(t,l,s) + sum(tdmp)/(nl(l)*dl(l));
                % FPR MP
                fprmp(t,l,s) = fprmp(t,l,s) + sum(sum(blockmp) - tdmp)/(nl(l)*(m-dl(l)));
            end
        end
    end
end

% Average performance metrics
truediscomp = truediscomp/nrep;
truediscmp = truediscmp/nrep;
tpromp = tpromp/nrep;
fpromp = fpromp/nrep;
tprmp = tprmp/nrep;
fprmp = fprmp/nrep;

% Save data
for s = 1:length(sigmans)
    
    dlmwrite(strcat('TDOMP-',num2str(sigmans(s),'%0.2f'),'.dat'),...
        [tausomp' squeeze(truediscomp(:,:,s))],'delimiter',' ');
    dlmwrite(strcat('TPROMP-',num2str(sigmans(s),'%0.2f'),'.dat'),...
        [tausomp' squeeze(tpromp(:,:,s))],'delimiter',' ');
    dlmwrite(strcat('FPROMP-',num2str(sigmans(s),'%0.2f'),'.dat'),...
        [tausomp' squeeze(fpromp(:,:,s))],'delimiter',' ');
    
    dlmwrite(strcat('TDMP-',num2str(sigmans(s),'%0.2f'),'.dat'),...
        [tausmp' squeeze(truediscmp(:,:,s))],'delimiter',' ');
    dlmwrite(strcat('TPRMP-',num2str(sigmans(s),'%0.2f'),'.dat'),...
        [tausmp' squeeze(tprmp(:,:,s))],'delimiter',' ');
    dlmwrite(strcat('FPRMP-',num2str(sigmans(s),'%0.2f'),'.dat'),...
        [tausmp' squeeze(fprmp(:,:,s))],'delimiter',' ');
    
    save('ROCtau.mat','sigmans','tausomp','truediscomp','tpromp','fpromp',...
        'tausmp','truediscmp','tprmp','fprmp');
end

