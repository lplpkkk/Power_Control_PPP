%% this program is for calculating the SINR gain for Power control policy 
% and modified by Z.K.Yang, 11/10
clear all
clc
%% parameter setting
lambda = 3e-4:3e-4:30e-4; % TX intensity
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
    outage_pc(i)=q_pc/Nsamp
    outage_npc(i)=q_npc/Nsamp
end

%% Plot section
plot(lambda(1:i),outage_pc(1:i),'b-s','LineWidth',2,'MarkerSize',10);
hold on
plot(lambda(1:i),outage_npc(1:i),'r-s','LineWidth',2,'MarkerSize',10);
xlabel('Intensity($\lambda$, Transmitters/m$^2$)');
ylabel('Mean SINR');
grid on;