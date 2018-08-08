function [User, max_diff,min_norm,max_norm]= BF_wmmse(L,User,Cells,Chn,lambda,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico)

%%% --------
%%% This function computes transmit beamformer for fixed weight and
%%% receiver beamformer
%%% --------
max_diff = 0;
ant_gain = 1;
Total_BS = L_Macro + L_Pico;
min_norm = inf;
max_norm = 0;
for l = 1:L
    for ik = 1:length(Cells(l).Scheduled_User)
        k = Cells(l).Scheduled_User(ik);
        
%         if length(lambda) == 1  %total power constraint
%             Kz_UL_lambda = eye(User(l,k).ServingAnt)*lambda;
%         elseif length(lambda) > 2   %per BS power constraint         
%             lambda_diag = [];
%             for iBS = 1:length(User(l,k).ServingCluster)
%                 curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
%                 curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
%                 if curr_BS == 1 %Macro BS
%                     lambda_diag = [lambda_diag;ones(Num_TxAnt_Macro,1)*lambda(User(l,k).ServingCluster(iBS))];
%                 else    %Pico BS
%                     lambda_diag = [lambda_diag;ones(Num_TxAnt_Pico,1)*lambda(User(l,k).ServingCluster(iBS))];
%                 end
%             end
%             
%             Kz_UL_lambda = diag(lambda_diag);
%             
%         end
        
%         Kz_UL_intf = zeros(size(Kz_UL_lambda));
%         for m = 1:L
%             for ikk = 1:length(Cells(m).Scheduled_User)
%                 kk = Cells(m).Scheduled_User(ikk);
%                 
%                 if (m ~= l) || (kk ~=k)
%                     h = (User(m,kk).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, m, kk, L_Macro,L_Pico,Num_TxAnt_Pico)';
%                     Kz_UL_intf = Kz_UL_intf + User(m,kk).wgt_mse*(h'*h)*ant_gain;
%                 else
%                     H = (User(m,kk).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, m, kk, L_Macro,L_Pico,Num_TxAnt_Pico)';
%                     Kz_UL_intf = Kz_UL_intf + User(m,kk).wgt_mse*(H'*H)*ant_gain;
%                 end
%                 
%             end
%             
%         end
%         
%         curr_norm = norm(Kz_UL_intf);
%         if curr_norm < min_norm
%             min_norm = curr_norm;
%         end
%         if curr_norm > max_norm
%             max_norm = curr_norm;
%         end
%         
%         Kz_UL = Kz_UL_lambda + Kz_UL_intf;
        H = (User(l,k).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)';
        invKz = pinv(User(l,k).Kz_UL);
        
        previous_beam = User(l,k).beam_tx;
        
        User(l,k).beam_tx= invKz*H'*sqrt(ant_gain)*User(l,k).wgt_mse;
        
        diff = norm(User(l,k).beam_tx - previous_beam)/norm(previous_beam);
        
        if diff > max_diff
            max_diff = diff;
        end
        
    end
    
end