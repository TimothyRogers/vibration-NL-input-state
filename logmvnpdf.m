function [loglik] = logmvnpdf(x,mu,Sig,LS)

%  Tim Rogers 2020 - tim.rogers@sheffield.ac.uk


% x - N x D matrix of observations
% mu - 1(N) x D matrix of mean(s)
% Sig - D x D Covariance matrix


D = size(x,2);

% if all(size(Sig)~=D)
%     error('Require square covariance matrix of size D')
% end
if nargin < 4
    LS = chol(Sig,'lower');
end
res = bsxfun(@minus,x,mu);
R = res/LS';

loglik = - D/2*log(2*pi) ...
    - sum(log(diag(LS))) ...
    - 0.5*sum(R.^2, 2);

% res = x-mu;
% v = res/(LS');
% loglik = -D/2*log(2*pi) ...
%          -sum(log(diag(LS))) ...
%          -0.5*diag(v*v');
