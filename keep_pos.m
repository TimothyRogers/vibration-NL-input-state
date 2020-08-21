function [xp] = keep_pos(x,reverse)

% TR 2020

if nargin < 2
    reverse = false;
end

if all(x < 1000)  
    if ~reverse
        xp = log(exp(x)-1);
    else
        xp = log(1+exp(x));
    end
else
    xp = x;
end

% if reverse
%     xp = exp(x);
% else
%     xp = log(x);
% end

end