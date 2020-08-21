function [ wn ] = normLogWeight( w )

%  Tim Rogers 2020 - tim.rogers@sheffield.ac.uk


wstar = max(w);
wn = exp(w-wstar);
wn = wn./sum(wn);
end

