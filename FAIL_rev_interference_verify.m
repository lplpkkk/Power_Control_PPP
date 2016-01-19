clear all
clc

% created by Z.K.Yang, 11/07
% this file is for verifying the problem in the letter paper

%% parameter setting
lambda = 4e-4:1e-4:15e-4;
SimRad=500;
SimArea=SimRad^2*pi;

sam = 1e5;
NumSam=500;
NumDrops=500;
% D0=20;
D0 = 15 + 10 * rand(1,NumSam);

alpha=4;
P=1;

D = 15 + 10 * rand(1,sam);
G = exprnd(1,1,sam);
H = exprnd(1,1,sam);

v1 = 2.5;
v2 = 5;
tau = pi * gamma(1-2/alpha) * gamma (1+2/alpha);
%% main simulation part
for m=1:length(lambda)
    MeanTX=round(lambda(m)*SimArea);
    NumTX=poissrnd(MeanTX,1,NumDrops); 
    temp=0; % this is for original 
    temp2 = 0 ;% this is for pc
    temp3 = 0; % this is for 

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
		t1_c = exp(-I0);
        t1 = exp(-I1_sim);
        
		temp=temp+ mean(1./I0);
		temp2 = temp2 +mean(1./I1_sim);

	end
	q(m)=temp/NumDrops;
	q2(m)=temp2/NumDrops;
end

t2_h = mean(exprnd(1,1,NumSam).^(2*v1/alpha));
t2_d = mean((15+10*rand(1,NumSam)).^(2*v2/alpha));
r2_c = exp(-lambda*tau);
r2 = exp(-lambda*tau*t2_h*t2_d);

result_c = gamma(1+alpha/2) ./ (pi*lambda*gamma(1-2/alpha)*1*gamma(1+2/alpha)).^(alpha/2);
result_pc = gamma(1+alpha/2) ./ (pi*lambda*gamma(1-2/alpha)*t2_h*t2_d*gamma(1+2/alpha)).^(alpha/2);
%% plot section
% plot(lambda,q,'ob','LineWidth',5); % const sim
% hold on
% plot(lambda,result_c,'b--'); % theoretical, origin
% hold on
plot(lambda,q2,'or','LineWidth',1.5); % pc sim
hold on
plot(lambda,result_pc,'r--'); % theoretical, pc

grid on
legend('Simulation Result','theoretical result');
xlabel('Node Intensity (\lambda)');
ylabel('mean(1/I)');
grid;