%% this program is the old one for power control.
% and modified by Z.K.Yang, 11/17 
clear all
clc
%% parameter setting
lambda = 3e-4:1e-4:3e-4;
SimRad=500;
SimArea=SimRad^2*pi;
NumSam=5000;
NumDrops=50;
D=20;
alpha=4;
P=1;
G=1;
tau = pi * P * G * gamma(1-2/alpha);
beta=1;
S = P * G * D^(-alpha);

q=zeros(1,length(lambda)); 
rate = zeros(1,length(lambda));
qth=q;
c=2*pi*gamma(2/alpha)*gamma(1-2/alpha)/alpha;
qth=1-exp(-lambda.*D^2.*c.*beta^(2/alpha)); %Theoretical result of the outage probability

for m=1:length(lambda)
    MeanTX=round(lambda(m)*SimArea);
    NumTX=poissrnd(MeanTX,1,NumDrops); 
    temp=0;
    temp2=0;
	for k=1:NumDrops
		Rad=sqrt(rand(1,NumTX(k)))*SimRad;
		Ang=rand(1,NumTX(k))*2*pi;
		Xt=Rad.*cos(Ang)+j*Rad.*sin(Ang);
		Ht=exprnd(1,NumTX(k),NumSam);

		%% 10/22 0:08 
		D_PC = 15+10*rand(NumTX(k),NumSam);
		v_1 = 2.5; % for fading
		v_2 = 5; % for distance 
		abs_power = gamma(1+v_1)*3.8792*10^(6);
		P_PC = Ht.^(v_1) .* D_PC.^(v_2)/abs_power;

		D_t = 15+10*rand(1,NumSam);
		h_t = exprnd(1,1,NumSam);
		P_t = h_t.^(v_1) .* D_t.^(v_2)/abs_power;

		S0=P_t.*h_t.*D_t.^(-alpha); % signal power
		I0=abs(Xt).^(-alpha)*(Ht.*P_PC); %interference power
		%%

		S_o=h_t.*D_t.^(-alpha);
		I_o=abs(Xt).^(-alpha)*Ht; 

		SIR0=S0./I0;
		SIR_o=S_o./I_o;    
		temp=temp+mean(I0);
		temp2=temp2+mean(SIR_o);
	end
	q(m)=temp/NumDrops
	q0(m)=temp2/NumDrops;
end

%% plot section;
plot(lambda,q,'b-o');
hold on
plot(lambda,q0,'r-o');
% for i =1:length(lambda)
%     fun = @(y) (1-exp(-y))./y .*exp(-lambda(i)*tau.*S^(-2/alpha).*y.^(2/alpha)); 
%     rate(i) = integral(fun,0,inf);
% end

legend('Simulation Result','Integration result');
xlabel('Node Intensity (\lambda)');
ylabel('Ergodic rate');
grid;