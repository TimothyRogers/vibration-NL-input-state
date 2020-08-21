function [pf] = lfmlikelihood12(hyps,m,k,c,k3,y,dt,F,t,xstar,Np)

% Tim Rogers 2020 - tim.rogers@sheffield.ac.uk


hyps_real = keep_pos(hyps,1);
sf2 = hyps_real(1);
ll = hyps_real(2);

lambda = 1./ll;
q = 2*sf2./ll;
% Agp = expm([0 1;-lambda^2 -2*lambda]*dt);
% 
% emldt = exp(-dt*lambda);
% Agp = emldt*[(1+dt*lambda), dt; -dt*lambda^2, (1-dt*lambda)];
% Qgp = [sf2*(1-(dt^2*lambda^2 + (1+dt*lambda)^2)*exp(-2*dt*lambda)), 2*sf2*dt^2*lambda^3*exp(-2*dt*lambda);
%      2*sf2*dt^2*lambda^3*exp(-2*dt*lambda), sf2*(lambda^2-(dt^2*lambda^4+lambda^2*(1-dt*lambda)^2)*exp(-2*dt*lambda))];

Agp = exp(-1/ll*dt);
Qgp = q-Agp*q*Agp;

sp2 = 0;
LQ = [sqrt(dt)*sqrt(sp2)*eye(2),[0;0];[0 0],chol(Qgp,'lower')];
Q = LQ*LQ';

duffing_deriv = @(xx,tt) [xx(:,3)/m-k/m*xx(:,1)-c/m*xx(:,2)-k3/m*xx(:,1).^3];
fx = @(xx,tt) xx + dt.* [xx(:,2),...
    xx(:,3)/m-k/m*xx(:,1)-c/m*xx(:,2)-k3/m*xx(:,1).^3,...
    -lambda*xx(:,3)];


% Q = [5e-8*eye(2); LQ = chol(Q,'lower');
proposal = @(xx,tt) fx(xx,tt) + (LQ*randn(size(xx')))';
proposalpdf = @(xt,xx,tt) logmvnpdf(xt,fx(xx,tt),Q+1e-16*eye(3));

% Likelihood is a Gaussian observe displacment and velocity
R =0.01*cov(y');
% R = 0.01*eye(3);
% R(3,3) = 1;
obs_mod = 'accel';
switch obs_mod
    case 'disp'
        likelihood = @(yy,xx,tt) logmvnpdf(yy,xx(:,1),R(1,1));
    case 'vel'
        likelihood = @(yy,xx,tt) logmvnpdf(yy,xx(:,2),R(2,2));
    case 'accel'
        likelihood = @(yy,xx,tt) logmvnpdf(yy,duffing_deriv(xx,tt),R(3,3));
    case 'force'
        likelihood = @(yy,xx,tt) logmvnpdf(yy,xx(:,3),0.01^2);
    case 'dispandvel'
        likelihood = @(yy,xx,tt) logmvnpdf(yy,xx(:,1:2),R(1:2,1:2));
end

P0gp = sf2;
x0 = [0;0;0];%sqrt(k/m)*2*pi];
P0y = 1e-6*cov(y(1:2,:)');
P0 = [P0y,zeros(2,1);zeros(1,2),P0gp];

dF = NaN(size(F));

% psamps = x0 + chol(P0,'lower')*randn(size(x0,1),1000);
% figure
% for pp = 1:4
%     subplot(4,1,pp)
%     histogram(psamps(pp,:))
%     if ~isnan(yt(pp,1))
%         xline(yt(pp,1),'r','LineWidth',2);
%     end
% end

%%
% Assign Proposal
dynamics.proposal = proposal;
dynamics.proposalpdf = proposalpdf;
% Assign Likelihood
dynamics.likelihood = likelihood;
% Select quantity to observe
switch obs_mod
    case 'disp'
        dynamics.y = y(1,:) + sqrt(R(1,1))*randn(size(y(1,:)));
    case 'vel'
        dynamics.y = y(2,:) + sqrt(R(2,2))*randn(size(y(2,:)));
    case 'accel'
        dynamics.y = y(3,:) + sqrt(R(3,3))*randn(size(y(3,:)));  
    case 'force'
        dynamics.y = y(4,:) + 0.01*randn(size(y(4,:)));
    case 'dispandvel'
        dynamics.y = y(1:2,:);
end

% Number of particles
props.N_particles = Np;
% Time Vector
props.t = t;
% Systematic Resampling
props.sampler = 'syst';
% Prior Sampling
props.prior = @(nn) repmat(x0',nn,1) + (chol(P0,'lower')*randn(3,nn))';

% Reference trajectory
% props.xstar = yt';
% props.xstar(:,4) = 150*randn(size(F));
props.xstar = xstar;

pf = bootstrapCPF_AS(dynamics,props);
% priorpdf = @(xx) logmvnpdf(xx,x0',P0);
lik = pf.logLik;


end