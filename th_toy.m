%% Parameter setting
ro=20;
ave=20^ro;
alpha=4;
sam=1e6;
counter=0;
counter_npc=0;
counter_pc=0;
%% main code
for i=1:sam
%     D=unifrnd(15,25,1,1);
    D=15+10*rand(1);
    PL_npc=mean(D.^(-alpha));
    PL_pc=mean((D.^(-alpha+ro))/ave);
    % PC_gain = mean(PL_pc)/mean(PL_npc)
    PC_gain = PL_pc/PL_npc;
    
    if(PC_gain>1)
        counter=counter+1;
       
        counter_pc=counter_pc+1;
    end
    
    if(PL_npc/1e-5>1)
        counter_npc=counter_npc+1;
    end
    
    if(PL_pc/1e-5>1)
        counter_pc=counter_pc+1;
    end
    
end

Psuc_gain=counter/sam
Psuc_pc=counter_pc/sam
Psuc_npc=counter_npc/sam