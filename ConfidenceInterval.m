function [data] = ConfidenceInterval(S0)
meanS0 = mean(S0);
stdS0 = std(S0);
pos = meanS0 + 2 * stdS0;
neg = meanS0 - 2 * stdS0;
ciTwoSigma = [pos neg];
SEM = stdS0 / sqrt(length(S0));
ts = tinv([0.025  0.975],length(S0)-1);
ci95 = meanS0 + ts * SEM;
data = zeros(2, 2);
data(1,:) = ciTwoSigma;
data(2,:) = ci95;
end

