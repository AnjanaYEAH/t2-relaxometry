function fSumRes = TwoCompSSD(x, vY, vEchoTime)

% extract parameters
S0 = x(1)^2;
v1 = 1/(x(2)^2+1);
T2_1 = x(3)^2;
T2_2 = x(4)^2;

v2 = 1-v1;

% Synthesize the signals according to the model
vSignal = S0*(v1*exp(-vEchoTime/T2_1) + v2*exp(-vEchoTime/T2_2));

% Compute the sum of square differences
fSumRes = sum((vY - vSignal).^2);
