clear all
clc

% created by Z.K.Yang, 11/29
% this file is for veifying theoretical interference

%% parameter setting
sam = 1e6;
P = 1;    % constant power case
H = exprnd(1,1,sam);
D = 15 + 10 * rand(1,sam);
% D = 15 + 5 * exprnd(1,1,sam);

alpha = 4;
G = exprnd(1,1,sam);

v1 = 2.5;
v2 = 5;

v1_opt = -3:0.1:3;
v2_opt = 0:0.1:5;

lambda = 4e-4:1e-4:15e-4;
SIR_opt = zeros(1,length(lambda));
%% main code - mean SIR, theoretical result
% S1 and I1 are presented as numerator and denominator
S1 = P * mean(H) * mean(D.^(-alpha)) * gamma(1+alpha/2);
I1 = ( pi * lambda * gamma(1-2/alpha) * mean(P^(2/alpha)) * mean(G.^(2/alpha)) ).^(alpha/2);
SIR1 = S1./I1;

%% main code - mean SIR, power control policy
S2 = mean(H.^(v1+1)) * mean(D.^(v2-alpha)) * gamma(1+alpha/2);
I2 = ( pi * lambda * gamma(1-2/alpha) * mean( ( (H.^(v1)) .* (D.^(v2)) ).^(2/alpha) ) * mean(G.^(2/alpha)) ).^(alpha/2);
SIR2 = S2./I2;

% theoretical result for power control
I2_the = gamma(2+v1)*mean(D.^(v2))*2*pi*lambda/(alpha-2);
%% plot section
% plot(lambda, I1,'b--');
% hold on
% plot(lambda, I2, 'r--');
% hold on
plot(lambda, I2_the, 'ro');

I2_the./I2