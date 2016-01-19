clear all
clc

% created by Z.K.Yang, 11/04
% this file is for verifying how much gain will provide in power control policy

%% parameter setting
sam = 1e6;
P = 1;    % constant power case
H = exprnd(1,1,sam);
D = 15 + 10 * rand(1,sam);
% D = 15 + 5 * exprnd(1,1,sam);

alpha = 4;
G = exprnd(1,1,sam);

v1 = 2.5;
v2 = 5;

v1_opt = -3:0.1:3;
v2_opt = 0:0.1:5;

lambda = 4e-4:1e-4:15e-4;
SIR_opt = zeros(1,length(lambda));
%% main code - mean SIR, theoretical result
S1 = P * mean(H) * mean(D.^(-alpha)) * gamma(1+alpha/2);
I1 = ( pi * lambda * gamma(1-2/alpha) * mean(P^(2/alpha)) * mean(G.^(2/alpha)) ).^(alpha/2);
SIR1 = S1./I1;

%% main code - mean SIR, power control policy
S2 = mean(H.^(v1+1)) * mean(D.^(v2-alpha)) * gamma(1+alpha/2);
I2 = ( pi * lambda * gamma(1-2/alpha) * mean( ( (H.^(v1)) .* (D.^(v2)) ).^(2/alpha) ) * mean(G.^(2/alpha)) ).^(alpha/2);
SIR2 = S2./I2;

%% main code - this is for finding optimal power control policy
S_temp = mean(H.^(v1_opt(1)+1)) * mean(D.^(v2_opt(1)-alpha)) * gamma(1+alpha/2);
I_temp = ( pi * lambda * gamma(1-2/alpha) * mean( ( (H.^(v1_opt(1))) .* (D.^(v2_opt(1))) ).^(2/alpha) ) * mean(G.^(2/alpha)) ).^(alpha/2);
SIR_temp = S_temp./I_temp;


parfor k = 1:length(lambda)

	for i = 1:length(v1_opt)
		for j = 1:length(v2_opt)

				S_opt = mean(H.^(v1_opt(i)+1)) * mean(D.^(v2_opt(j)-alpha)) * gamma(1+alpha/2);
				I_opt = ( pi * lambda(k) * gamma(1-2/alpha) * mean( ( (H.^(v1_opt(i))) .* (D.^(v2_opt(j))) ).^(2/alpha) ) * mean(G.^(2/alpha)) ).^(alpha/2);
				SIR_opt(k) = S_opt/I_opt;

				if (SIR_opt(k) > SIR_temp(k))
						SIR_temp(k) = SIR_opt(k);
						i_opt = i;
						j_opt = j;
				end
		end
	end

	disp(['when lambda is:' num2str(lambda(k)) ',the optimal policy is:' num2str(v1_opt(i_opt)) ',' num2str(v2_opt(j_opt))]);
end
%% plot section
plot (lambda,SIR1,'b-o');
hold on 
plot (lambda,SIR2,'r-o');
hold on
plot(lambda,SIR_temp,'g-o','LineWidth',3);
grid on;

ylabel('mean SIR');
xlabel('\lambda');
legend('original result');