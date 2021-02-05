function [ self, r ] = calc_dist( self, xp, xq )
% Calculate pariwise distances for Kernels

switch self.dmetric
    case 'euclidean'
            %r = pdist2(xp,xq);
            r = pdist2(bsxfun(@rdivide,xp,sqrt(self.hyps.ll)),bsxfun(@rdivide,xq,sqrt(self.hyps.ll)));
    case 'manhattan'
        
    otherwise
        if isdouble(self.dmetric)
            
        else
            error('Unknown Distance Metric')
        end
end

if nargout<2
    self.r = r;
end

end

