clear all
clc
%% 11/29, verifying approximate distribution of interference, gamma distribution

% parameter setting
sam=1e6;
lambda = 1e-4;
alpha =4;
v1 = 2.5;
v2 = 5;
% D = 15 + 10*rand(1,sam);
D = 20;
beta = 1;
k = 4*gamma(2+v1)^2*(alpha-1)*mean(D.^(v2))*pi*lambda/gamma(3+2*v1)/(alpha-2)^2
theta = gamma(3+2*v1)*(alpha-2)/(2*gamma(2+v1)*(alpha-1))
a = beta^(1/(1+v1))*D^((v2-alpha)/(v1+1))

P_suc = 0;
for n =1:200
    P_suc = P_suc + (-a)^n/factorial(n) *theta^(k+( n/(1+v1) ))*gamma(k+( n/(1+v1) ));
end

P_suc