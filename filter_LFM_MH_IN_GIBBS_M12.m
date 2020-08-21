%% Test LFM Nonlinear

% Tim Rogers 2020 - tim.rogers@sheffield.ac.uk

clear all
close all
clc

rng(1)

%% Simulate Nonlinear System

fs = 2048;
dt = 1/fs;
secs = 0.5;
t = 0:dt:(secs-dt);
T = length(t);

m = 1;
k = 1e4;
c = 20;
k3 = 1e9;

F = sin(sqrt(k/m)*2*pi*t);
dF = sqrt(k/m)*2*pi*cos(sqrt(k/m)*2*pi*t);


sf2 = 20;
ll = 0.1;
kern = Matern32(struct('sf2',sf2,'ll',ll,'sn2',1e-12),0);
K = kern.calc_K(t',t');

F = chol(K,'lower')*randn(length(F),1); F = F';

lambda = 1./ll;
q = 2*sf2./ll;

Agp = exp(-1/ll*dt);
Qgp = q-Agp*q*Agp;

sp2 = 0;
LQ = [sqrt(dt)*sqrt(sp2)*eye(2),[0;0];[0 0],chol(Qgp,'lower')];
Q = LQ*LQ';

duffing_deriv = @(xx,tt) [xx(:,3)/m-k/m*xx(:,1)-c/m*xx(:,2)-k3/m*xx(:,1).^3];
fxx = @(xx,tt) [xx(:,2),...
    xx(:,3)/m-k/m*xx(:,1)-c/m*xx(:,2)-k3/m*xx(:,1).^3,...
    -lambda*xx(:,3)];
X = NaN(1,3,T+1);
X(1,:,1) = zeros(1,3);


for tt = 1:T
    X(:,:,tt+1) = X(:,:,tt) + dt*fxx(X(:,:,tt),tt) + (LQ*randn(3,1))';
end

x = squeeze(X(1,1,2:end))';
xd = squeeze(X(1,2,2:end))';
F = squeeze(X(1,3,2:end))';
xdd = duffing_deriv(squeeze(X(1,:,2:end))',1:T)';

y = [x;xd;xdd;F];
y0 = zeros(3,1);

figure
for pp = 1:3
    
    subplot(4,1,pp)
    plot(squeeze(X(1,pp,2:end)));
       
end

yt = squeeze(X(:,:,2:end))';


%% Set Up Filter Using Euler Now With Added LFM

% States 1 & 2 are dynamics 3 & 4 are forcing
K = 1e5;
Np = 50;

nav = 1;

hyps = [sf2,ll];
lik = NaN(K,length(Np),size(hyps,1));
lik_curr = NaN(K,nav);
lik_prop = NaN(K,nav);

% Fix random seed
rng(2)

xsamps = NaN(T,3,K);
hyp_est = NaN(K,2);

% Initial guess of hyperparameters perturbed
hyp_est(1,:) = keep_pos([sf2,ll]+[10 -0.09]);


props.xstar = yt;
xstar = zeros(size(yt));
xsamps(:,:,1) = xstar;

mcmc_prop = @(xx) xx + [1,1e-1].*randn(1,2);

hyp_prior = @(xx)1;% logmvnpdf(xx,keep_pos([sf2,ll]),0.01*eye(2));

hyp_prop=hyp_est(1,:);
lik(1) = prop_lik12(hyp_prop,m,k,c,k3,dt,xstar);
lik(1) = -1e8;
stepp = 5;

%% Running metropolis in Gibbs for PF
naccept = 0;
for kk = 1:K
    
    %% Gibbs Step
    
    fprintf('Iteration: %i / %i\t',kk,K);
    
    if mod(kk,stepp) == 0 || kk == 1
        
        if mod(kk,100)==0
            figure(201)
            subplot(411)
            plot(lik(2:end))
            hold on
            plot(lik_prop(2:end))
            hold off
            subplot(4,1,2:4)
            semilogy(keep_pos(hyp_est(1:kk,:),1)./hyps)
            drawnow
        end
    
    pf = lfmlikelihood12(hyp_est(kk,:),m,k,c,k3,y,dt,F,t,xstar,Np);
    ii = resamp(pf.w(:,end),'multi',1);
    xsamps(:,:,kk+1) = squeeze(pf.xpaths(ii,:,:))';

    xstar = squeeze(xsamps(:,:,kk+1));
        %xstar = squeeze(sum(bsxfun(@times,pf.w(:,end),pf.xpaths),1))';

    else
        
        xsamps(:,:,kk+1) = xsamps(:,:,kk);

    end
    
    
    %% M-H Step
    
 
        xlik = xstar(2:end,:);

        hyp_prop =  mcmc_prop(hyp_est(kk,:));
        lik_curr(kk,:) =  prop_lik12(hyp_est(kk,:),m,k,c,k3,dt,xlik);
        [lik_prop(kk,:)] = prop_lik12(hyp_prop,m,k,c,k3,dt,xlik);
        
        alpha = min(exp(lik_prop(kk,:)-lik_curr(kk,:)),1);
        
        if  rand() < alpha
            
            naccept = naccept+1;
            
            fprintf('ACCEPT\t %i/%i',naccept,kk)
            
            hyp_est(kk+1,:) = hyp_prop;
            lik(kk+1) = lik_prop(kk,:);
            
        else
            
            fprintf('REJECT\t %i/%i',naccept,kk)
            
            hyp_est(kk+1,:) = hyp_est(kk,:);
            lik(kk+1) = lik_curr(kk);
            
            
        end
        
    
    fprintf('\n')
    
    
end


%%

figure(101)
KK = kk;
spoint = kk/2;

for pp = 1:3
    subplot(3,1,pp)
    hold off
    plot(squeeze(xsamps(:,pp,spoint:stepp:KK)),'Color',[0.1 0.1 1 0.2])
    hold on
    plot(yt(:,pp),'r')
    plot(nanmean(squeeze(xsamps(:,pp,spoint:KK)),2),'Color',[0 1 0])
%     ylim([1.5*min(yt(:,pp)) 1.5*max(yt(:,pp))])
end
subplot(311)
title('Sampled Paths')
linkaxes(get(gcf,'children'),'x')


%%


figure
plot(lik)
hold on
plot(lik_prop)

figure
for nn = 1:length(Np)
    subplot(length(Np),1,nn)
    for ss = 1:size(hyps,1)
        histogram(lik(:,nn,ss),100,'EdgeColor','none');
        hold on
    end
    
end


figure
for nn = 1:length(Np)
    subplot(length(Np),1,nn)
    for ss = 1:size(hyps,2)
        histogram(keep_pos(hyp_est(:,ss),1),100,'EdgeColor','none');
        hold on
    end
    
end

%%

figure
plot(keep_pos(hyp_est,1)./hyps)

