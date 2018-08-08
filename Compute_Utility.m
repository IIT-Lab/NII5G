function [Log_Utility,User,Log_Data] = Compute_Utility(Num_Cell,K,User,Cells,iter,posNumber,Log_Data)

%%% ------------------------------------------------------------
%%% This function computes log utility according to current User rate
%%% ------------------------------------------------------------
alpha = 0.9;
Log_Utility = 0;
backhaul = 0;
sumrate = 0;
wgt_sumrate = 0;
for iCell = 1:Num_Cell
    cnt_scheduled_user = 0;
    
    Log_Data(posNumber,iter).Cells(iCell).Scheduled_User = Cells(iCell).Scheduled_User;
    
    for k = 1:K
        
        if sum(ismember(k,Cells(iCell).Scheduled_User)) > 0 %if user k is scheduled
            User(iCell,k).long_avg_rate = alpha*User(iCell,k).long_avg_rate + (1-alpha)*User(iCell,k).inst_rate;
            User(iCell,k).avg_rate = (User(iCell,k).avg_rate*(iter - 1) + User(iCell,k).inst_rate)/iter;
            cnt_scheduled_user = cnt_scheduled_user + 1;
            Log_Data(posNumber,iter).User(iCell,cnt_scheduled_user).rate = User(iCell,k).inst_rate;
            Log_Data(posNumber,iter).User(iCell,cnt_scheduled_user).Cluster = User(iCell,k).ServingCluster;
            Log_Data(posNumber,iter).User(iCell,cnt_scheduled_user).beam_tx = User(iCell,k).beam_tx;
            Log_Data(posNumber,iter).User(iCell,cnt_scheduled_user).beam_rx = User(iCell,k).beam_rx;
            backhaul = backhaul + User(iCell,k).inst_rate*length(User(iCell,k).ServingCluster);
            sumrate = sumrate + User(iCell,k).inst_rate;
            wgt_sumrate = wgt_sumrate + User(iCell,k).inst_rate*User(iCell,k).wgt_wsr;
        else
            
            User(iCell,k).long_avg_rate = alpha*User(iCell,k).long_avg_rate;
            User(iCell,k).avg_rate = User(iCell,k).avg_rate*(iter - 1)/iter;
            
        end
        
        Log_Utility = Log_Utility + log10(User(iCell,k).long_avg_rate);

    end
end
    
Log_Data(posNumber,iter).utility = Log_Utility;
Log_Data(posNumber,iter).sum_backhaul = backhaul;
Log_Data(posNumber,iter).sumrate = sumrate;
Log_Data(posNumber,iter).wgt_sumrate = wgt_sumrate;