%% this program is for calculating the SINR gain for Power control policy 
% and modified by Z.K.Yang, 11/10
clear all
clc
%% parameter setting
lambda = 4e-4:1e-4:15e-4; % TX intensity
L=550;      
SimArea=(2*L)^2; % Area of the network
Nsamp=6500; % Generate Nsamp samples for each intensity
alpha=4;
MeanSINR=zeros(1,length(lambda));
m=2; % Nakagami-m Fading 
Dmax=25;
Dmin=15;
rho=1.5; u=5;
sigma20=1e-10;
meanPt=1/mean((gamrnd(m,1/m,1e5,1).^rho).*((Dmin+(Dmax-Dmin)*rand(1e5,1)).^u));
%% main code
for i=1:length(lambda)
	NumTX=poissrnd(lambda(i)*SimArea,[Nsamp,1]);
	SumSINR=0;
    SumSINR_npc=0;
    SumI=0;
    SumS=0;
    q_pc=0;
    q_npc=0;
	for k=1:Nsamp
		% Generate the positions
		Xt=unifrnd(-L,L,NumTX(k),2);
		DisXt=sqrt(Xt(:,1).^2+Xt(:,2).^2); 
		Indx=find(DisXt<1);
		G=gamrnd(m,1/m,[NumTX(k),1]); % Generate fading samples of interference channels
		G(Indx)=0;
        % Generate the interference with pc
		Di=Dmin+(Dmax-Dmin)*rand(NumTX(k),1);
		Hi=gamrnd(m,1/m,[NumTX(k),1]);
		Pi=(Hi.^rho).*(Di.^u)*meanPt;
        % Generate the signal with pc
		H=gamrnd(m,1/m,1); 
		D=Dmin+(Dmax-Dmin)*rand(1);
		Pt=(H^rho)*(D^u)*meanPt;
        % calculate the SIR
        SumI=SumI+sum(Pi.*G.*(DisXt.^(-alpha)))+sigma20;
        SumS=SumS+Pt*H*D^(-alpha);
		SumSINR=SumSINR+sum((Pt*H*D^(-alpha))./(sum(Pi.*G.*(DisXt.^(-alpha)))+sigma20));
        SumSINR_npc=SumSINR_npc+sum((H*D^(-alpha))./(sum(G.*(DisXt.^(-alpha)))+sigma20));
 
        SINR_pc=sum((Pt*H*D^(-alpha))./(sum(Pi.*G.*(DisXt.^(-alpha)))+sigma20));
        if (SINR_pc>1)
            q_pc=q_pc+1;
        end
        % this part is for non-pc
        SINR_npc=sum((H*D^(-alpha))./(sum(G.*(DisXt.^(-alpha)))+sigma20));
       if(SINR_npc>1)
           q_npc=q_npc+1;
       end
	end
	% Calculate Mean SINR for each intensity
	MeanSINR(i)=SumSINR/Nsamp;
    MeanSINR_npc(i)=SumSINR_npc/Nsamp;
end

%% theoretical section

sam = 1e5;
NumSam=500;
D = 15 + 10 * rand(1,sam);
v1 = 1.5;
v2 = 5;
% _c : constant power control scheme
rev_interference_c = gamma(1+alpha/2) ./ (pi*lambda*gamma(1-2/alpha)*1*gamma(1+2/alpha)).^(alpha/2);
signal_c = mean(D.^(-alpha));
SIR_c = signal_c.*rev_interference_c;
% _pc : channel-aware power control scheme
t2_h = mean(exprnd(1,1,NumSam).^(2*v1/alpha));
t2_d = mean((15+10*rand(1,NumSam)).^(2*v2/alpha));
rev_interference_pc = gamma(1+alpha/2) ./ (pi*lambda*gamma(1-2/alpha)*t2_h*t2_d*gamma(1+2/alpha)).^(alpha/2);
signal_pc = mean(exprnd(1,1,NumSam).^(1+v1))*mean(D.^(v2-alpha));
SIR_pc = signal_pc.*rev_interference_pc;

%% Plot section
plot(lambda(1:i),MeanSINR(1:i),'bo','LineWidth',2,'MarkerSize',10);
hold on
plot(lambda(1:i),SIR_pc(1:i),'b--','LineWidth',2,'MarkerSize',10);
hold on
plot(lambda(1:i),MeanSINR_npc(1:i),'ro','LineWidth',2,'MarkerSize',10);
hold on
plot(lambda(1:i),SIR_c(1:i),'r--','LineWidth',2,'MarkerSize',10);

xlabel('Intensity($\lambda$, Transmitters/m$^2$)');
ylabel('Mean SINR');
grid on;