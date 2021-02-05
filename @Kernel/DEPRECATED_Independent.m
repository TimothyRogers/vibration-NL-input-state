function [ K ] = Independent( self, xp, xq, hyps, idiff )

if nargin< 5
    idiff = 0;
end

if size(xp,2) == size(xq,2)
    D = size(xp,2);
else
    error('Input dimensionality not consistent');
end

if nargin < 4
    sf2 = self.hyps.sf2;
    ll = self.hyps.ll;
    sn2 = self.hyps.sn2;
else
    sf2 = hyps.sf2;
    ll = hyps.ll;
    sn2 = hyps.sn2;
end

if ~isempty(self.active_dims)
    xp = xp(:,self.active_dims);
    xq = xq(:,self.active_dims);
end
    

npts_p = size(xp,1);
npts_q = size(xq,1 );


type = self.subtype;

K = zeros(npts_p,npts_q);
for i = 1:D
    hyps_temp.sf2 = sf2(i);
    hyps_temp.ll = ll(i);
    hyps_temp.sn2 = 0;
    K = K+self.(type)(xp(:,i),xq(:,i),hyps_temp,idiff);
end

K = K+eye(size(K))*sn2;


end

