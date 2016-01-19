%% this program is for verifying some bugs..
% and modified by Z.K.Yang, 12/21
clear all
clc
%% parameter setting
a = -
%% main code
for i=1:length(lambda)
	for k=1:Nsamp
        d1=15+10*rand(1);
        d2=15+10*rand(1);
		a=a+d1*d2;
        b=
    end
   
end

%% Plot section
plot(lambda(1:i),outage_pc(1:i),'b-s','LineWidth',2,'MarkerSize',10);
hold on
plot(lambda(1:i),outage_npc(1:i),'r-s','LineWidth',2,'MarkerSize',10);
xlabel('Intensity($\lambda$, Transmitters/m$^2$)');
ylabel('Outage Prob');
grid on;