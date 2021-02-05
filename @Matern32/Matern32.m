classdef Matern32 < Kernel
    
    % MATERN32 - Class defining Matern32 covariance function, inherits from
    % the Kernel superclass
    
    properties
        hyps
    end
    
    methods
        
        % Constructor
        function self = Matern32(hyps,ARD,active_dims)
            
            % Call superconstructor
            if nargin<3
                active_dims = [];
            end
            self = self@Kernel(hyps,ARD,active_dims);
            
            % Hyps error handling - needs refining
            if nargin > 0
                self.hyps = hyps;
            else
                error('Must Specify Hyperparameters');
            end
            
            
        end
        
        % Get full K matrix
        function [ K ] = calc_K( self, varargin )
            
            % Enforce data to be considered
            if nargin < 2
                error('Need data input to kernel')
            else
                xp = varargin{1};
                xq = varargin{2};
            end
            
            % Default derivative to none
            if nargin< 5
                idiff = 0;
            else
                idiff = varargin{4};
            end
            
            % Check input dimentisonality
            if size(xp,2) == size(xq,2)
                D = length(self.active_dims);
            else
                error('Input dimensionality not consistent');
            end
            
            % Use kernel hyps unless others are specified
            if nargin < 4
                sf2 = self.hyps.sf2;
                ll = self.hyps.ll;
                sn2 = self.hyps.sn2;
            else
                
                sf2 = varargin{3}.sf2;
                ll = varargin{3}.ll;
                sn2 = varargin{3}.sn2;
            end
            
            % If no noise then set to 0
            if isempty(sn2)
                sn2 = 0;
            end
            
            % If active dims not defined used all dims
            if ~isempty(self.active_dims)
                xp = xp(:,self.active_dims);
                xq = xq(:,self.active_dims);
            end
            
            % Get points in each data set
            npts_p = size(xp,1);
            npts_q = size(xq,1 );
            
            % If not ARD use ISO kernel
            if self.ARD == false
                
                % Pairwise distances
                r = pdist2(xp,xq);
                
                % Calc K
                switch idiff
                    case 0
                        % K
                        K = sqrt(3).*r./ll;
                        K = sf2*((1+K).*exp(-K))+sn2*eye(npts_p,npts_q);
                        
                    case 1
                        % dK/dsf2
                        K = sqrt(3)*r./ll;
                        K = sf2*(1+K).*exp(-K);
                    case 2
                        % dK/dll
                        K = sqrt(3*r/ll);
                        K = (3*sf2.*r.^2.*exp(-K))./ll^2;
                    case 3
                        % dK/dsn2
                        K = sn2*eye(npts_p,npts_q);
                end
                
                
                
            else
                
                
                %% Matern Kernel nu = 3/2 with ARD
                
                % Pairwise with ARD lengthscales
                D = length(ll);
                r = pdist2(bsxfun(@rdivide,xp,sqrt(ll)),bsxfun(@rdivide,xq,sqrt(ll)));
                
                % Calc K
                switch idiff
                    case 0
                        % K
                        K = sf2*((1+r*sqrt(3)).*exp(-sqrt(3)*r))+sn2*eye(npts_p,npts_q);
                        
                    case 1
                        % dK/d(ln(sf2))
                        K = sf2*((1+r*sqrt(3)).*exp(-sqrt(3)*r))+sn2*eye(npts_p,npts_q);
                        
                    case D+2
                        % dK/d(ln(sn2))
                        K = sn2*eye(npts_p,npts_q);
                        
                    otherwise
                        % dK/d(ln(ll))
                        
                        xl2 = pdist2(xp(:,idiff-1)./ll(idiff-1),xp(:,idiff-1)./ll(idiff-1)).^2;
                        expK = exp(-sqrt(3)*r);
                        K = -sqrt(3).*xl2.*expK + sqrt(3).*xl2.*(1+sqrt(3)*r)*expK;
                        K = K*ll(idiff-1);
                end
                
                
            end
            
        end
        
        % Get diagonal of covariance matrix - r = 0 so diagonal is signal
        % variance plus noise variance
        function Kd = diag_K(self,varargin)
            % Handle no noise case
            if isempty(self.hyps.sn2)
                self.hyps.sn2 = 0;
            end
            Kd = self.hyps.sf2+self.hyps.sn2;
        end
        
    end
    
end

