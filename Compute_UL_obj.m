function [UL_obj,User] = Compute_UL_obj(L,User,Cells,lambda,Chn,L_Macro,L_Pico,Num_TxAnt_Pico,per_BS_power_constraint_matrix)

%%
UL_obj = sum(sum(per_BS_power_constraint_matrix.*lambda));
Total_BS = L_Macro+L_Pico;
for l= 1:L
    for ik = 1:length(Cells(l).Scheduled_User)
        k = Cells(l).Scheduled_User(ik);
        
        H = (User(l,k).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)';
        diag_lambda = zeros(User(l,k).ServingAnt,User(l,k).ServingAnt);
        for iBS = 1:length(User(l,k).ServingCluster)
                curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
                curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
                diag_lambda = diag_lambda + lambda(curr_cell,curr_BS)*User(l,k).ServeBS(iBS).CoopMatrix;
        end
            
        User(l,k).Kz_UL = User(l,k).Kz_UL_intf + diag_lambda;
        
        UL_obj = UL_obj + real((User(l,k).wgt_mse)^2*H*pinv(User(l,k).Kz_UL)*H');
        
    end
end
       
