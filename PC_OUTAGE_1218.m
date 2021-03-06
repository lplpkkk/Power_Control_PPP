%% this program is for calculating the SINR gain for Power control policy 
% and modified by Z.K.Yang, 11/10
clear all
clc
%% parameter setting
lambda = 10e-4:1e-4:15e-4; % TX intensity
L=550;      
SimArea=(2*L)^2; % Area of the network
Nsamp=6500; % Generate Nsamp samples for each intensity
alpha=4;
MeanSINR=zeros(1,length(lambda));
m=1; % Nakagami-m Fading 
Dmax=25;
Dmin=15;
Dave=(Dmax+Dmin)/2;
rho=0; u=20;
P_comp=Dave^u;
sigma20=1e-10;
meanPt=1/mean((gamrnd(m,1/m,1e5,1).^rho).*((Dmin+(Dmax-Dmin)*rand(1e5,1)).^u));
%% main code
for i=1:length(lambda)
	NumTX=poissrnd(lambda(i)*SimArea,[Nsamp,1]);
    Signal1=0;
    Signal2=0;
	SumSINR=0;
    SumSINR_npc=0;
    SumI=0;
    SumS=0;
    q_pc=0;
    q_npc=0;
    T5=0;
    T6=0;
    signal_gain=0;
    interference_gain=0;
    sinr_gain=0;
    sinr_gain2=0;
	for k=1:Nsamp
        
		% Generate the positions
		Xt=unifrnd(-L,L,NumTX(k),2);
		DisXt=sqrt(Xt(:,1).^2+Xt(:,2).^2); 
		Indx=find(DisXt<1);
		G=gamrnd(m,1/m,[NumTX(k),1]); % Generate fading samples of interference channels
		G(Indx)=0;
        
        % Generate the interference with pc
		Di=Dmin+(Dmax-Dmin)*rand(NumTX(k),1);
% 		Hi=gamrnd(m,1/m,[NumTX(k),1]);
		Pi=(Di.^u);
        
        % Generate the signal with pc
		H=gamrnd(m,1/m,1); 
		D=Dmin+(Dmax-Dmin)*rand(1);
        Pt=D.^u;
        
%         Signal1 = Signal1+mean( (D.^(-alpha+u))/mean(Di.^(u)) );
%         Signal2 = Signal2+D.^(-alpha);      
%         MagGain=mean(Signal1)/mean(Signal2)
        
        test_th1 = mean( (D.^(-alpha+u))/ P_comp);
        test_th2 = mean(D.^(-alpha));
        MagGain_th = test_th1/test_th2;

        % calculate the SIR
%         SumI=SumI+sum(Pi.*G.*(DisXt.^(-alpha)))+sigma20;
        SumI=SumI+sum(G.*(DisXt.^(-alpha)))+sigma20;
        SumS=SumS+Pt*H*D^(-alpha);
		SumSINR=SumSINR+sum((Pt*H*D^(-alpha))./(sum(Pi.*G.*(DisXt.^(-alpha)))+sigma20));
        SINR_pc=(Pt*H*D^(-alpha))./(sum(Pi.*G.*(DisXt.^(-alpha)))+sigma20);
        T1=(sum(G.*(Pi.*DisXt.^(-alpha)))+sigma20);
        T3=sum(G.*Pi.*DisXt.^(-alpha));
        T5=T5+Pt*H*D^(-alpha);
        if (SINR_pc>1)
            q_pc=q_pc+1;
        end
        
        % this part is for non-pc
        SINR_npc=(H*D^(-alpha))./(sum(G.*(DisXt.^(-alpha)))+sigma20);
        
        SumSINR_npc=SumSINR_npc+sum((H*D^(-alpha))./(sum(G.*(DisXt.^(-alpha)))+sigma20));
        T2=(sum(G.*(DisXt.^(-alpha)))+sigma20);
        T4=sum(G.*DisXt.^(-alpha));
        T6=T6+H*D^(-alpha);
        signal_gain=signal_gain+(T5/T6);
        interference_gain=interference_gain+(T1);
        a1=(SINR_pc/SINR_npc);
        sinr_gain = sinr_gain+(SINR_pc/SINR_npc);
        sinr_gain2 = sinr_gain2 + (signal_gain/interference_gain);
         (signal_gain/interference_gain);
       if(SINR_npc>1)
           q_npc=q_npc+1;
       end
	end
	% Calculate Mean SINR for each intensity
    IG = interference_gain/Nsamp;
    MeanI=SumI/Nsamp
    MeanS1 = T5/Nsamp;
    MeanS2 = T6/Nsamp;
    SG=MeanS1/MeanS2;
    SINRG2 = sinr_gain2/Nsamp;
    SINRG = sinr_gain/Nsamp;
	MeanSINR(i)=SumSINR/Nsamp;
    MeanSINR_npc(i)=SumSINR_npc/Nsamp;
    outage_pc(i)=q_pc/Nsamp;
    outage_npc(i)=q_npc/Nsamp;
end

%% Plot section
plot(lambda(1:i),outage_pc(1:i),'b-s','LineWidth',2,'MarkerSize',10);
hold on
plot(lambda(1:i),outage_npc(1:i),'r-s','LineWidth',2,'MarkerSize',10);
xlabel('Intensity($\lambda$, Transmitters/m$^2$)');
ylabel('Outage Prob');
grid on;