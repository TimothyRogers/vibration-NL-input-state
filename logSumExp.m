function [lse] = logSumExp(x)

%  Tim Rogers 2020 - tim.rogers@sheffield.ac.uk


mx = max(x);
lse = log(sum(exp(x-mx)))+mx;


end