function [User, Cells] = Direct_Beamformers(Num_Cell,K,Cells,User,Chn,noise_dBm,L_Macro,...
    L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,per_BS_power_constraint_matrix,Num_Serving_User)

leakage = 1;
Power_Alloc = 1;
Total_BS = L_Pico+L_Macro;
for l = 1:Num_Cell
    
    for iBS = 1:Total_BS
        
        Num_Users_Served = Num_Serving_User(l,iBS); % Number of Users Served by this BS
        if Num_Users_Served == 0
            continue
        end
        BS_cluster = zeros(Num_Users_Served,2); % List of Users Served by this BS
        count = 1;
        for iCell = 1:Num_Cell
            for ik = 1:length(Cells(iCell).Scheduled_User)
                k = Cells(iCell).Scheduled_User(ik);
                if User(iCell,k).ServingCluster == Total_BS*(l-1)+iBS
                    BS_cluster(count,1) = iCell;
                    BS_cluster(count,2) = k;
                    count = count + 1;
                end
            end
        end
        
        for k_index = 1:Num_Users_Served
            iCell = BS_cluster(k_index,1);
            k = BS_cluster(k_index,2);
            Hconcat = [];
            
            % Local leakage
            if leakage == 1
                for kk_index = 1:Num_Users_Served
                    if kk_index ~= k_index
                        iC = BS_cluster(kk_index,1);
                        kk = BS_cluster(kk_index,2);    % User identification: (iCell,k)
                        if iBS == 1
                            Hconcat = [Hconcat Get_Chn(Chn, User(iC,kk).ServingCluster, iC, kk, L_Macro,L_Pico,Num_TxAnt_Macro)];    % Channel Matrix Between User (iCell,k) and its Serving eRRH (l,iAnt)
                        else
                            Hconcat = [Hconcat Get_Chn(Chn, User(iC,kk).ServingCluster, iC, kk, L_Macro,L_Pico,Num_TxAnt_Pico)];    % Channel Matrix Between User (iCell,k) and its Serving eRRH (l,iAnt)
                        end
                    end
                end
                
                if iBS == 1
                    Hi = Get_Chn(Chn, User(iCell,k).ServingCluster, iCell, k, L_Macro,L_Pico,Num_TxAnt_Macro);
                else
                    Hi = Get_Chn(Chn, User(iCell,k).ServingCluster, iCell, k, L_Macro,L_Pico,Num_TxAnt_Pico);
                end
                
                if Num_Users_Served == 1
                    Hconcat = zeros(size(Hi,1),size(Hi,1));
                end
            
            % Global leakage
            elseif leakage == 2
                for iC = 1:Num_Cell
                    for ikk = 1:length(Cells(iC).Scheduled_User)
                        kk = Cells(iC).Scheduled_User(ikk);
                        if kk ~= k || iC ~= iCell
                            if iBS == 1
                                Hconcat = [Hconcat Get_Chn(Chn, User(iCell,k).ServingCluster, iC, kk, L_Macro,L_Pico,Num_TxAnt_Macro)];    % Channel Matrix Between User (iC,kk) and eRRH (l,iAnt)
                            else
                                Hconcat = [Hconcat Get_Chn(Chn, User(iCell,k).ServingCluster, iC, kk, L_Macro,L_Pico,Num_TxAnt_Pico)];    % Channel Matrix Between User (iC,kk) and eRRH (l,iAnt)
                            end
                        end
                    end
                end
                
                if iBS == 1
                    Hi = Get_Chn(Chn, User(iCell,k).ServingCluster, iCell, k, L_Macro,L_Pico,Num_TxAnt_Macro);
                else
                    Hi = Get_Chn(Chn, User(iCell,k).ServingCluster, iCell, k, L_Macro,L_Pico,Num_TxAnt_Pico);
                end
                
                if K == 1 && Num_Cell == 1
                    Hconcat = zeros(size(Hi,1),size(Hi,1));
                end
                
            end
            
            Hconcat = Hconcat * (Hconcat');
            if Power_Alloc == 1
                Mconcat = inv(10^(noise_dBm/10)*eye(size(Hconcat,1))/(User(iCell,k).power_sqrt^2) + Hconcat)*(Hi*(Hi'));
            else
                Mconcat = inv(10^(noise_dBm/10)*eye(size(Hconcat,1))/(per_BS_power_constraint_matrix(l,iBS)*Num_Users_Served) + Hconcat)*(Hi*(Hi'));
            end
            [V,D] = eig(Mconcat);
            [~,permutation] = sort( abs(diag(D)), 'descend' );
            beam = V(:,permutation(1));
            if Power_Alloc == 1
                User(iCell,k).beam_tx = User(iCell,k).power_sqrt * beam;
            else
                User(iCell,k).beam_tx = sqrt(per_BS_power_constraint_matrix(l,iBS)/Num_Users_Served)*beam;
            end
        end
    end
end