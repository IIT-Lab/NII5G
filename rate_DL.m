function [User, sumrate, wgt_sumrate, Per_BS_Bkhaul,curr_user_num,Cells,max_rate_diff] = ...
    rate_DL(L,Q,User,Cells,Chn,noise_dB,intf_2nd,L_Macro,L_Pico,Num_TxAnt_Pico)

%%% -----------------------------------------------------------------------
%%% This function computes rate for given transmit beamformers
%%% -----------------------------------------------------------------------
ant_gain = 1;
gap = 1;
sumrate = 0;
wgt_sumrate = 0;
Total_BS = L_Macro+L_Pico;
Per_BS_Bkhaul = zeros(L,Total_BS);
curr_user_num = 0;
max_rate_diff = 0;
for l = 1:L
    user_index_to_delete = [];
    for ik = 1:length(Cells(l).Scheduled_User)
        k = Cells(l).Scheduled_User(ik);
        Kz_DL = eye(Q)* 10^(noise_dB/10) + intf_2nd;
        for m = 1:L
            for ikk = 1:length(Cells(m).Scheduled_User)
                kk = Cells(m).Scheduled_User(ikk);
                try
                    beam = User(m,kk).beam_tx;
                    H_k = Get_Chn(Chn,User(m,kk).ServingCluster,l,k,L_Macro,L_Pico,Num_TxAnt_Pico);
                    if (m ~= l) || (kk ~=k)
                        Kz_DL = Kz_DL + H_k' * (beam*beam') * H_k * ant_gain;
                    else
                        H = H_k'*beam;
                    end
                end
            end
        end
        try
            invKz = inv(Kz_DL);
            
            sinr = real(H'*invKz*H)*ant_gain;
            User(l,k).inst_rate = log2(1 + sinr/gap);
            
            if abs(User(l,k).inst_rate - User(l,k).prev_inst_rate)/User(l,k).prev_inst_rate > max_rate_diff
                max_rate_diff = abs(User(l,k).inst_rate - User(l,k).prev_inst_rate)/User(l,k).prev_inst_rate;
            end
            
            User(l,k).prev_inst_rate = User(l,k).inst_rate;
            User(l,k).wgt_rate = User(l,k).wgt_wsr*User(l,k).inst_rate;
            
            sumrate = sumrate + User(l,k).inst_rate;
        end
        wgt_sumrate = wgt_sumrate + User(l,k).wgt_rate;
        %       SumBackhaul = SumBackhaul + User(l,k).inst_rate*length(User(l,k).ServingCluster);
        
        for iBS = 1:length(User(l,k).ServingCluster)
            curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
            curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
            
            Per_BS_Bkhaul(curr_cell,curr_BS) = Per_BS_Bkhaul(curr_cell,curr_BS) + User(l,k).inst_rate;
            
        end
        
        if User(l,k).inst_rate < 1e-2
            user_index_to_delete = [user_index_to_delete ik];
        end
        
    end
    
    Cells(l).Scheduled_User(user_index_to_delete) = [];
    curr_user_num = curr_user_num + length(Cells(l).Scheduled_User);
end


