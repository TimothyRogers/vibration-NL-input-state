function [lik] = prop_lik12(hyps,m,k,c,k3,dt,xstar)

% Tim Rogers 2020 - tim.rogers@sheffield.ac.uk


% sf2 = exp(hyps(1));
% ll = exp(hyps(2));

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
fx = @(xx,tt) xx + dt*[xx(:,2),...
    xx(:,3)/m-k/m*xx(:,1)-c/m*xx(:,2)-k3/m*xx(:,1).^3,...
    -lambda*xx(:,3)];

% Q = [5e-8*eye(2); LQ = chol(Q,'lower');
proposal = @(xx,tt) fx(xx,tt) + (LQ*randn(size(xx')))';
proposalpdf = @(xt,xx,tt) logmvnpdf(xt,fx(xx,tt),Q+1e-16*eye(3));

px = (Agp*(xstar(2:end-1,3))')';

% lik = sum(logmvnpdf(xstar(3:end,3),xstar(2:end-1,3)+dt*xstar(2:end-1,4),Qgp(1,1)));

% lik = sum(logmvnpdf(xstar(3:end,3),(Agp*(xstar(2:end-1,3))')',Qgp));
lik = sum(proposalpdf(xstar(3:end,:),(xstar(2:end-1,:)),Q));

end