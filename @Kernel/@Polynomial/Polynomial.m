function [ K ] = Polynomial( self, xp, xq, hyps, idiff )

% POLYNOMIAL kernel function

if nargin< 5
    idiff = 0;
end

if ~isempty(self.p)
    p = self.hyps.p;
else
    p = 3;
end

if size(xp,2) == size(xq,2)
    D = size(xp,2);
else
    error('Input dimensionality not consistent');
end

if nargin < 4
    % Get hyps from self
    sf2 = self.hyps.sf2;
    c = self.hyps.c;
    sn2 = self.hyps.sn2;
else
    % Get from hyps
    sf2 = hyps.sf2;
    c = hyps.c;
    sn2 = hyps.sn2;    
end
    

if ~isempty(self.active_dims)
    xp = xp(:,self.active_dims);
    xq = xq(:,self.active_dims);
end
    
npts_p = size(xp,1);
npts_q = size(xq,1 );

if self.ARD == false
    
    ip = xp*xq';
    
    switch idiff
        case 0
            K = sf2*(ip+c).^p;
            K = K + sn2*eye(npts_p,npts_q);
        case 1
            K = 2*sf2*(ip + c).^p;
        case 2
            K = c*p*sf2*(ip + c).^(p-1);            
    end
    
else
    
    error('Not yet implemented');
    
end
   



end

