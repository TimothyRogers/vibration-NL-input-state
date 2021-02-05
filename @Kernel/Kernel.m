classdef Kernel < handle
    %KERNEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nhyp
        type
        ARD = false;
        active_dims = [];
    end
    
    methods
        function self = Kernel(hyps, ARD, active_dims)
            
            if nargin ~= 0
                
                if nargin > 2 && ~isempty(ARD)
                    self.ARD = ARD;
                else
                    self.ARD = false;
                end
                
                hypnames = fieldnames(hyps);
                for i = 1:length(hypnames)
                    hyps_cell{i} = hyps.(hypnames{i});
                end
                
                size_matrix(1,:) = cellfun(@(x) size(x,1),hyps_cell,'uni',false);
                size_matrix(2,:) = cellfun(@(x) size(x,2),hyps_cell,'uni',false);
                size_matrix = cell2mat(size_matrix);
                
                self.nhyp = sum(size_matrix(2,:));
                
                
                if nargin > 2 && ~isempty(active_dims)
                    self.active_dims = active_dims;
                end
                
            end
            
        end
        
        function self = set_hyps(self,hyps)
            self.hyps = hyps;
        end
        
        function self = plus(self,K2)
            self = SumKern(self,K2);
        end
        
        function self = times(self,K2)
            self = MultKern(self,K2);
        end
        
        function self = mtimes(self,K2)
            self = MultKern(self,K2);
        end
        
        %        K = RBF(self,xp,xq,hyps,idiff);
        %        K = Matern12(self,xp,xq,hyps,idiff);
        %        K = Matern32(self,xp,xq,hyps,idiff);
        %        K = Matern52(self,xp,xq,hyps,idiff);
        %        K = Linear(self,xp,xq,hyps,idiff);
        %        K = Polynomial(self,xp,xq,hyps,idiff,degree);
        %        K = Independent(self,xp,xq,hyps,idiff);
        
    end
    
    methods (Static)
        
    end
    
end

