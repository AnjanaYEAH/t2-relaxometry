function [a] = AIC(RESNORM, N, K)
%k is number of measurements
%N is the number of parameters, we add 1 since we are estimating variance
%S is the signal
%A is the actual measurements
%RESNORM is sum of sum((A - S)^2)
a = (2 * (N+1)) + (K * log(RESNORM/K));
end