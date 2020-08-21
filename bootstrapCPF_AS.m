function [pf] = bootstrapCPF_AS(dynamics,props)

% Conditional Bootstap Particle Filter 
%
%  Tim Rogers 2020 - tim.rogers@sheffield.ac.uk


prop = dynamics.proposal; % Proposal distribution sampler from f(x,t)
proppdf = dynamics.proposalpdf; % Proposal distribution pdf f(x,t)
lik = dynamics.likelihood; % Log Likelihood model p(y | g(x,t))
y = dynamics.y; % Measurements

xstar = props.xstar; % Trajectory to condition on
N = props.N_particles; % Number of particles

t = props.t; % Time vector
T = length(t);
sampler = props.sampler; % Resampling strategy
x0 = feval(props.prior,N-1); % Samples from prior
w0 = lik(y(1),x0,t(1)); % Weights of prior
% w0 = w0./sum(w0); % Normalise

[~,D] = size(x0);

% Preallocate
x = NaN(N,D,T);
w = NaN(N,T);
a = NaN(N,T);

% Initial Conditions
x(1:N-1,:,1) = x0; % Samples from prior
x(N,:,1:T) = xstar'; % Start of reference trajectory
xpaths = x;
w(1:N-1,1) = w0;
w(end,1) = lik(y(1),xstar(1,:),t(1));
w(:,1) = normLogWeight(w(:,1));
a(:,1) = 1:N;
logLik = 0;


for k = 2:T
    
    %% RESAMPLE-MOVE
    
    % Resample good particles
    a(1:N-1,k) = resamp(w(:,k-1),sampler,N-1);
    
    
    % Propgation except reference
    x(1:N-1,:,k) = prop(squeeze(x(a(1:N-1,k),:,k-1)),t(k-1));
    
    % Last particle is reference
    x(N,:,k) = xstar(k,:);
    
       
    %% Ancestor Sampling
    
    % Propagate Reference
    xp = [squeeze(x(1:N-1,:,k));prop(xstar(k-1,:),t(k-1))];
    wanc = proppdf(squeeze(x(N,:,k)),xp,t(k));
    a(N,k) = resamp(normLogWeight(wanc),'multi',1);
    
    xpaths(:,:,1:k-1) = xpaths(a(:,k),:,1:k-1);
    xpaths(:,:,k) = x(:,:,k);
      
    %% Weighting
    w(1:N,k) = lik(y(k),squeeze(x(1:N,:,k)),t(k));
    logLik = logLik + logSumExp(w(:,k)) - log(N); % Log Likelihood of Filter
    w(:,k) = normLogWeight(w(:,k));
      
end

xpaths(:,:,end) = x(:,:,end);
for k = T-1:-1:1    
    xpaths(:,:,k) = x(a(:,k+1),:,k);
end
    


% Outputs
pf.a = a;
pf.x = x;
pf.w = w;
pf.xpaths = xpaths;
pf.logLik = logLik;


end