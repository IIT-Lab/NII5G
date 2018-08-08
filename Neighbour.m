function [User Asso_User] = Neighbour(Num_Cell, K, path_loss_dB, per_BS_power_constraint_matrix,User,M,noise_dBm)

%%% -----------------------------------------------------------------------
%%% This function finds the nearest M neighbouring BS's for every user
%%% accoring to the pathloss and BS trasnmit power budget
%%% -----------------------------------------------------------------------
Total_BS = 4;
%User = struct([]);
% Tx_Power = [P_max_Macro_dBm*ones(1,L_Macro) P_max_Pico_dBm*ones(1,L_Pico)];
% Tx_Power_Matrix = kron(Tx_Power,ones(Num_Cell,1));
Asso_User = zeros(M,Num_Cell,Total_BS);
%num of users claiming each BS as the strongest BS
for iCell = 1:Num_Cell
    for k = 1:K
        
        path_loss = squeeze(path_loss_dB(:,:,iCell,k));
        Rcv_Power = dec2dB(per_BS_power_constraint_matrix) - path_loss;
        Rcv_Power_vec =reshape(Rcv_Power',[],1);
        
        [~,SortedBS] = sort(Rcv_Power_vec,'descend');
        User(iCell,k).NeighbourCluster = SortedBS(1:M);
        
        
        for m = 1:M
            Cell_index = ceil(SortedBS(m)/Total_BS);
            BS_index = mod(SortedBS(m)-1,Total_BS)+1;
            User(iCell,k).Neighbour(m).Cell_BS = [Cell_index BS_index];
            
            Asso_User(m, Cell_index, BS_index) = Asso_User(m, Cell_index, BS_index) + 1;
            
        end
        
    end
end
