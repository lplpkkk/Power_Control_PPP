clear all
clc

tic
%% Parameter Setting Section
sam = 1e6;
d = 15+10*rand(1,sam);
d2 = 15+10*rand(1,sam);
alpha = 4;
m=1;
Dmax=25;
Dmin=15;
v = -2:1:100;
lambda = 1e-4:1e-4:15e-4;
L=550;      
SimArea=(2*L)^2; % Area of the network
Nsamp=5000; % Generate Nsamp samples for each intensity
s_pc=zeros(1,length(lambda));s_npc=s_pc;
u=20;
sigma20=1e-10;

%% Section for numerical try tau
for i=1:length(v)
    tau_gain(i) = mean(d.^(alpha)) / ( mean(d.^(alpha-v(i))) * mean(d.^(v(i))) );
    tau_gain(i) = tau_gain(i)^(2/alpha);
end

%% Theoretical Section
for i = 1:length(lambda)
    P_suc(i) = mean( exp(-lambda(i)*pi*pi/2*d.^(2)) );
    P_suc2(i) = mean( exp(-lambda(i)*pi*pi/2*sqrt( d.^(4-v(23)).*d2.^(v(23)) ) ) );
    P_suc3(i) = mean( exp(-lambda(i)*pi*pi/2*sqrt( d.^(4-v(30)).*d2.^(v(30)) ) ) );
    P_suc4(i) = mean( exp(-lambda(i)*pi*pi/2*sqrt( d.^(4-v(40)).*d2.^(v(40)) ) ) );
end

%% Simulation Section
for i=1:length(lambda)
	NumTX=poissrnd(lambda(i)*SimArea,[Nsamp,1]);
    SumI=0;SumS=0;SumSINR=0;SumSINR_npc=0;
	for k=1:Nsamp
        
		% Generate the positions
		Xt=unifrnd(-L,L,NumTX(k),2);
		DisXt=sqrt(Xt(:,1).^2+Xt(:,2).^2); 
		Indx=find(DisXt<1);
		G=gamrnd(m,1/m,[NumTX(k),1]); % Generate fading samples of interference channels
		G(Indx)=0;
        
        % Generate the interference with pc
		Di=Dmin+(Dmax-Dmin)*rand(NumTX(k),1);
		Pi=(Di.^u);
        
        % Generate the signal with pc
		H=gamrnd(m,1/m,1); 
		D=Dmin+(Dmax-Dmin)*rand(1);
        Pt=D.^u;
        
        % Test..checking distance aware power control scheme gain
%         test_th1 = mean( (D.^(-alpha+u))/ P_comp);
%         test_th2 = mean(D.^(-alpha));
%         MagGain_th = test_th1/test_th2;

        % Calculate the SIR
		SumSINR=SumSINR+sum((Pt*H*D^(-alpha))./(sum(Pi.*G.*(DisXt.^(-alpha)))+sigma20));
%         SINR_pc=(Pt*H*D^(-alpha))./(sum(Pi.*G.*(DisXt.^(-alpha)))+sigma20);
         SINR_pc=(Pt*H*D^(-alpha))./1e-5;
         if (SINR_pc>1)
            s_pc(i)=s_pc(i)+1;
         end
         
         % Calculate the SIR_npc
        SumSINR_npc=SumSINR_npc+sum((H*D^(-alpha))./(sum(G.*(DisXt.^(-alpha)))+sigma20));
%         SINR_npc=(H*D^(-alpha))./(sum(G.*(DisXt.^(-alpha)))+sigma20);  
        SINR_npc=(H*D^(-alpha))./1e-5;
        if (SINR_npc>1)
           s_npc(i)=s_npc(i)+1;
        end
        
         % Test..checking magnitude of mean interference and signal power
        SumI=SumI+sum(Pi.*G.*(DisXt.^(-alpha)))+sigma20;
        SumS=SumS+Pt*H*D^(-alpha);
%         T1=(sum(G.*(Pi.*DisXt.^(-alpha)))+sigma20);
%         T3=sum(G.*Pi.*DisXt.^(-alpha));
%         T5=T5+Pt*H*D^(-alpha);        
%         T2=(sum(G.*(DisXt.^(-alpha)))+sigma20);
%         T4=sum(G.*DisXt.^(-alpha));
%         T6=T6+H*D^(-alpha);
%         signal_gain=signal_gain+(T5/T6);
%         a1=(SINR_pc/SINR_npc);
%         sinr_gain = sinr_gain+(SINR_pc/SINR_npc);
%         sinr_gain2 = sinr_gain2 + (signal_gain/interference_gain);
%          (signal_gain/interference_gain);
       
    end
    
	% Calculate Mean of everything
    suc_pc(i)=s_pc(i)/Nsamp;
    suc_npc(i)=s_npc(i)/Nsamp;
end

%% plot the success probability here, thoretical first

% Theoretical plotting
plot(lambda,P_suc,'r--','LineWidth',2);
hold on
plot(lambda,P_suc2,'b--','LineWidth',2);
hold on
plot(lambda,P_suc3,'m--','LineWidth',2);
hold on 
plot(lambda,P_suc4,'g--','LineWidth',2);
grid on

% Simulation plotting
plot(lambda, suc_npc,'ro','MarkerSize',8);
hold on
plot(lambda, suc_pc, 'bo','MarkerSize',8);
hold on

% Labelling and Legending
xlabel('Intensity(\lambda, Transmitters/m^2)');
ylabel('Success Probability');
legend('Wtihout Power Control','Power Control with \nu = 17','Power Control with \nu = 30','Power Control with \nu = 40','Simulation-non power control','Simulation-power control with 20');
toc