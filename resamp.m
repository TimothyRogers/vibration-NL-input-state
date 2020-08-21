function [ ind ] = resamp( x, strategy, ns )
%RESAMP resampling for SIR algorithm
%   Inputs:
%       x - vector of probabilities N x 1 must sum to 1
%       np - number of particles
%       strategy - string selecting resampling strategy
%           ~ 'multi' - multinomial resampling
%           ~ 'strat' - stratified resampling
%           ~ 'syst' - systematic resampling
%
%   Outputs:
%       ind - index of sampled particles
%
%  Tim Rogers 2020 - tim.rogers@sheffield.ac.uk


np = size(x,1);
if nargin < 3 || isempty(ns)
    ns = np;
end
cdf = cumsum(x); % Cumulative probability

switch lower(strategy)
    case 'multi'
        rn = rand(ns,1);% random draws for each particle
    case 'strat'
        rn = rand(ns,1);
        bnds = (linspace(0,np-1,np)./np)';
        rn = bnds+rn./np;
    case 'syst'
        rn = repmat(rand(1),ns,1);
        bnds = (linspace(0,np-1,ns)./np)';
        rn = bnds+rn./np;
end

ind = np+1-sum(bsxfun(@le,rn,cdf'),2);