clear all
clc

% created by Z.K.Yang, 11/28
% this file is for finding the worse/optimal parameter for PC, no optimal
% trying

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

v3 = -1.5;
v4 = 0;

lambda = 4e-4:1e-4:15e-4;
SIR_opt = zeros(1,length(lambda));
%% main code - mean SIR, theoretical result
S1 = P * mean(H) * mean(D.^(-alpha)) * gamma(1+alpha/2);
I1 = ( pi * lambda * gamma(1-2/alpha) * mean(P^(2/alpha)) * mean(G.^(2/alpha)) ).^(alpha/2);
SIR1 = S1./I1;

%% main code - mean SIR, power control policy
S2 = mean(H.^(v1+1)) * mean(D.^(v2-alpha)) * gamma(1+alpha/2);
I2 = ( pi * lambda * gamma(1-2/alpha) * mean( ( (H.^(v1)) .* (D.^(v2)) ).^(2/alpha) ) * mean(G.^(2/alpha)) ).^(alpha/2);
SIR2 = S2./I2;

%% main code - mean SIR, power control policy, with rho <0
S3 = mean(H.^(v3+1)) * mean(D.^(v4-alpha)) * gamma(1+alpha/2);
I3 = ( pi * lambda * gamma(1-2/alpha) * mean( ( (H.^(v3)) .* (D.^(v4)) ).^(2/alpha) ) * mean(G.^(2/alpha)) ).^(alpha/2);
SIR3 = S3./I3;

%% plot section
plot (lambda,SIR1,'b-o');
hold on 
plot (lambda,SIR2,'r-o');
hold on
plot (lambda, SIR3, 'm-o');