clear all
clc

% created by Z.K.Yang, 11/04
% this file is for verifying fitness of the
% simulation and the theoretical result

%% parameter setting
lambda = 4e-4:1e-4:15e-4;
SimRad=500;
SimArea=SimRad^2*pi;

sam = 1e5;
NumSam=50;
NumDrops=50;
% D0=20;
D0 = 15 + 10 * rand(1,NumSam);

alpha=4;
P=1;

D = 15 + 10 * rand(1,sam);
G = exprnd(1,1,sam);
H = exprnd(1,1,sam);

v1 = 2.5;
v2 = 5;
%% main simulation part
for m=1:length(lambda)
    MeanTX=round(lambda(m)*SimArea);
    NumTX=poissrnd(MeanTX,1,NumDrops); 
    temp=0; % this is for original 
    temp2 = 0 ;% this is for pc
	for k=1:NumDrops
		Rad=sqrt(rand(1,NumTX(k)))*SimRad;
		Ang=rand(1,NumTX(k))*2*pi;
		Xt=Rad.*cos(Ang)+j*Rad.*sin(Ang);

		H0=exprnd(1,1,NumSam);
		S0=H0.*D0.^(-alpha); % signal power
		Ht=exprnd(1,NumTX(k),NumSam);
		I0=abs(Xt).^(-alpha)*Ht; %interference power
		SIR0=S0./I0;  

		P0 = H0.^(v1).* D0.^(v2);
		Dt = 15 + 10 * rand(NumTX(k),NumSam);
		S1_sim = P0 .*S0;
		I1_sim = abs(Xt).^(-alpha)*(Ht.^(1+v1).*Dt.^(v2));
		SIR1_sim = S1_sim./I1_sim;

		temp=temp+ mean(SIR0);
		temp2 = temp2 +mean(SIR1_sim);

		% temp=temp+ mean(S0) / mean(I0);
		% temp2 = temp2 +mean(S1_sim) / mean(I1_sim);

	end
	q(m)=temp/NumDrops;
	q2(m)=temp2/NumDrops;
end

%% main code - mean SIR, theoretical result
S1 = P * mean(H) * mean(D.^(-alpha)) * gamma(1+alpha/2);
I1 = ( pi * lambda * gamma(1-2/alpha) * mean(P^(2/alpha)) * mean(G.^(2/alpha)) ).^(alpha/2);
SIR1 = S1./I1;

%% main code - mean SIR, power control policy
S2 = mean(H.^(v1+1)) * mean(D.^(v2-alpha)) * gamma(1+alpha/2);
I2 = ( pi * lambda * gamma(1-2/alpha) * mean( ( (H.^(v1)) .* (D.^(v2)) ).^(2/alpha) ) * mean(G.^(2/alpha)) ).^(alpha/2);
SIR2 = S2./I2;

%% plot section
plot(lambda,q,'ob','LineWidth',5); % q simulation,original
hold on
plot(lambda,SIR1,'b--'); % theoretical, origin
hold on
plot(lambda,q2,'or','LineWidth',1.5); % simulation, pc
hold on
% plot(lambda,SIR2,'r--'); % theoretical, pc

grid on
legend('Simulation Result','theoretical result');
xlabel('Node Intensity (\lambda)');
ylabel('SIR');
grid;