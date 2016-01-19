clear all
clc
% this code help us know what is best for distance-aware power control

tic
sam = 1e6;
d = 15+10*rand(1,sam);
d2 = 15+10*rand(1,sam);
% d = 20;
% d2 = 20;
alpha = 4;
v = -2:1:100;
for i=1:length(v)
    tau_gain(i) = mean(d.^(alpha)) / ( mean(d.^(alpha-v(i))) * mean(d.^(v(i))) );
    tau_gain(i) = tau_gain(i)^(2/alpha);
end

%% plot the success probability here
lambda = 1e-4:1e-4:15e-4;
for i = 1:length(lambda)
    P_suc(i) = mean( exp(-lambda(i)*pi*pi/2*d.^(2)) );
    P_suc2(i) = mean( exp(-lambda(i)*pi*pi/2*sqrt( d.^(4-v(17)).*d2.^(v(17)) ) ) );
    P_suc3(i) = mean( exp(-lambda(i)*pi*pi/2*sqrt( d.^(4-v(30)).*d2.^(v(30)) ) ) );
    P_suc4(i) = mean( exp(-lambda(i)*pi*pi/2*sqrt( d.^(4-v(40)).*d2.^(v(40)) ) ) );
end

plot(lambda,P_suc,'r--','LineWidth',2);
hold on
plot(lambda,P_suc2,'b--','LineWidth',2);
hold on
plot(lambda,P_suc3,'m--','LineWidth',2);
hold on 
plot(lambda,P_suc4,'g--','LineWidth',2);
grid on
xlabel('Intensity(\lambda, Transmitters/m^2)');
ylabel('Success Probability');
legend('Wtihout Power Control','Power Control with \nu = 17','Power Control with \nu = 30','Power Control with \nu = 40');
toc