function fSumRes = ThreeCompSSD(x, vY, vEchoTime)

% extract parameters
S0 = x(1)^2;
v1 = (sin(x(2))^2);
v2 = (sin(x(3))^2);
v3 = sin(x(4))^2;
T2_1 = x(5)^2;
T2_2 = x(6)^2;
T2_3 = x(7)^2;

a = [v1, v2, v3];
a = sort(a);
total = v1 + v2 + v3;
v3 = a(1)/total;
v2 = a(2)/total;
v1 = a(3)/total;
% Synthesize the signals according to the model
vSignal = S0*(v1*exp(-vEchoTime/T2_1) + v2*exp(-vEchoTime/T2_2) + v3*exp(-vEchoTime/T2_3));

% Compute the sum of square differences
fSumRes = sum((vY - vSignal).^2);
