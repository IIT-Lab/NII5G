function [Log_Utility,User,Log_Data] = Compute_Utility_2(Num_Cell,K,User,Cells,iter,posNumber,Log_Data)

%%% ------------------------------------------------------------
%%% This function computes log utility according to current User rate
%%% ------------------------------------------------------------
alpha = 0.9;
Log_Utility = 0;
sumrate = 0;
wgt_sumrate = 0;
for iCell = 1:Num_Cell
    
    for k = 1:K
        
        User(iCell,k).long_avg_rate = alpha*User(iCell,k).long_avg_rate + (1-alpha)*User(iCell,k).inst_rate;
        User(iCell,k).avg_rate = (User(iCell,k).avg_rate*(iter - 1) + User(iCell,k).inst_rate)/iter;
        Log_Data(posNumber,iter).User(iCell,k).rate = User(iCell,k).inst_rate;
        Log_Data(posNumber,iter).User(iCell,k).Cluster = User(iCell,k).ServingCluster;
        Log_Data(posNumber,iter).User(iCell,k).beam_tx = User(iCell,k).beam_tx;
        %Log_Data(posNumber,iter).User(iCell,k).beam_rx = User(iCell,k).beam_rx;
        %backhaul = backhaul + User(iCell,k).inst_rate;
        sumrate = sumrate + User(iCell,k).inst_rate;
        Log_Utility = Log_Utility + log10(User(iCell,k).long_avg_rate);
    end
end
    
Log_Data(posNumber,iter).utility = Log_Utility;
%Log_Data(posNumber,iter).sum_backhaul = backhaul;
Log_Data(posNumber,iter).sumrate = sumrate;

