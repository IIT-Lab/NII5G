function [User, Cells, Asso_User, Num_Serving_User, Log_Data] = Pool_clustering(User,Cells,Num_Cell, K, Chn, S, sigmaList, Total_BS, path_loss_dB, snr_gap_dB, intf_2nd, ...
    per_BS_power_constraint_matrix,per_BS_power_constraint,per_BS_bkhaul_constraint,noise_dBm,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,L_Macro,L_Pico,P_max_total_dBm,iter,tau_zero,lambda,posNumber,Log_Data,Mode,iK,Num_Iter)

% Path Loss Based Clustering
if Mode == 1
    %every scheduled user connects with its best BS (lowest PL)
    [User, Asso_User] = Neighbour(Num_Cell, K, path_loss_dB, per_BS_power_constraint_matrix,User,1,noise_dBm);
    Asso_User = squeeze(Asso_User);
    Num_Serving_User = zeros(Num_Cell,Total_BS); %Number of users each BS supports
    for iCell = 1:Num_Cell
        for k = 1:K
            User(iCell,k).ServingCluster = User(iCell,k).NeighbourCluster; 
            User(iCell,k).ServingAnt = 0;
            for iBS = 1:length(User(iCell,k).ServingCluster)
                curr_cell = User(iCell,k).Neighbour(iBS).Cell_BS(1);
                curr_BS = User(iCell,k).Neighbour(iBS).Cell_BS(2);
                User(iCell,k).bkhaul_wgt(curr_cell,curr_BS) = 0.01;
                Num_Serving_User(curr_cell,curr_BS) = Num_Serving_User(curr_cell,curr_BS) + 1;
                if curr_BS == 1 %Macro BS
                    User(iCell,k).ServingAnt = User(iCell,k).ServingAnt + Num_TxAnt_Macro;
                else
                    User(iCell,k).ServingAnt = User(iCell,k).ServingAnt + Num_TxAnt_Pico;
                end
            end
        end
    end
    for iCell = 1:Num_Cell
        Cells(iCell).Scheduled_User = 1:K ;    % Need to change to respect fronthaul capacity
    end
    
elseif Mode == 1.5
    M = 1;
    Num_Scheduled_User = K;
    [User, Asso_User] = Neighbour(Num_Cell, K, path_loss_dB, per_BS_power_constraint_matrix,User,M,noise_dBm);
    
    for iCell = 1:Num_Cell
        [Cells(iCell).Scheduled_User, Cells(iCell).UnScheduled_User]= Round_Robin(Cells(iCell).UnScheduled_User,Num_Scheduled_User);
        if isempty(Cells(iCell).UnScheduled_User)
            Cells(iCell).UnScheduled_User = 1:K;
        end
    end
    
    Num_Serving_User = zeros(Num_Cell,Total_BS); %num of users each BS supports at the begining
    for iCell = 1:Num_Cell
        for ik = 1:length(Cells(iCell).Scheduled_User)
            k = Cells(iCell).Scheduled_User(ik);
            User(iCell,k).ServingCluster = User(iCell,k).NeighbourCluster; 
            User(iCell,k).ServingAnt = 0;
            for iBS = 1:length(User(iCell,k).ServingCluster)
                
                curr_cell = User(iCell,k).Neighbour(iBS).Cell_BS(1);
                curr_BS = User(iCell,k).Neighbour(iBS).Cell_BS(2);
                User(iCell,k).bkhaul_wgt(curr_cell,curr_BS) = 0.01;
                Num_Serving_User(curr_cell,curr_BS) = Num_Serving_User(curr_cell,curr_BS) + 1;
                if curr_BS == 1 %Macro BS
                    User(iCell,k).ServingAnt = User(iCell,k).ServingAnt + Num_TxAnt_Macro;
                else
                    User(iCell,k).ServingAnt = User(iCell,k).ServingAnt + Num_TxAnt_Pico;
                end
            end
        end
    end
    
    %Initialize wgts
    max_wgt = 0;
    for iCell = 1:Num_Cell
        for ik = 1:length(Cells(iCell).Scheduled_User)
            k = Cells(iCell).Scheduled_User(ik);
            User(iCell,k).wgt_wsr = 1/User(iCell,k).long_avg_rate;
            User(iCell,k).inst_rate = 0.01;
            User(iCell,k).prev_inst_rate = 0.01;
            if User(iCell,k).wgt_wsr > max_wgt
                max_wgt = User(iCell,k).wgt_wsr;
            end
        end
    end
    
    %normalize user wgts by max_wgt
    for iCell = 1:Num_Cell
        for ik = 1:length(Cells(iCell).Scheduled_User)
            k = Cells(iCell).Scheduled_User(ik);
            User(iCell,k).wgt_wsr = User(iCell,k).wgt_wsr/max_wgt;
            
        end
    end
    
    % initialize TX beamformers and sparse Matrix
    for l = 1:Num_Cell
        for kk = 1:length(Cells(l).Scheduled_User)
            k = Cells(l).Scheduled_User(kk); %scheduled user in cell l
            Num_ServingAnt = User(l,k).ServingAnt;
            temp_tx = randn(Num_ServingAnt,1) + 1i*randn(Num_ServingAnt,1);
            
            ant_head = 1;
            temp_tx2 = [];
            for m = 1:length(User(l,k).ServingCluster)
                
                curr_cell = ceil(User(l,k).ServingCluster(m)/Total_BS);
                curr_BS = User(l,k).ServingCluster(m) - (curr_cell-1)*Total_BS;
                if curr_BS == 1
                    Num_Tx_Ant = Num_TxAnt_Macro;
                else
                    Num_Tx_Ant = Num_TxAnt_Pico;
                end
                ant_tail = ant_head + Num_Tx_Ant - 1;
                norm_tx = temp_tx(ant_head:ant_tail)/norm(temp_tx(ant_head:ant_tail));
                temp_tx2 = [temp_tx2;sqrt(per_BS_power_constraint_matrix(curr_cell,curr_BS)/Num_Serving_User(curr_cell,curr_BS))*norm_tx];
                ant_head = ant_tail+1;
                User(l,k).BS_Power(m) = per_BS_power_constraint_matrix(curr_cell,curr_BS)/Num_Serving_User(curr_cell,curr_BS);
            end

            User(l,k).beam_tx = temp_tx2;
            if length(User(l,k).beam_tx) ~= User(l,k).ServingAnt
                fprintf('$$$ Tx Initialization Not Right $$$ \n');
            end
            User(l,k).SparseMatrix = eye(User(l,k).ServingAnt);
        end
    end
    
    %check per BS power constraint
    
    [sumpower, per_BS_power, User] = Sum_Power(Num_Cell,Cells,User,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
    
    [User,sumrate, wgt_sumrate, Per_BS_Bkhaul, Cells] = WSR_Solver_Zeronorm(Num_Cell,K,Cells,User,P_max_total_dBm,Chn,noise_dBm,snr_gap_dB,intf_2nd,...
        L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,2,per_BS_power_constraint,per_BS_power_constraint_matrix,per_BS_bkhaul_constraint,iter,tau_zero,lambda); 
    %Log_Data(posNumber,iter+1).Clustering_sumrate = sumrate;
    %Log_Data(posNumber,iter+1).Clustering_wgt_sumrate = wgt_sumrate;
    for l = 1:Num_Cell
        for k = 1:K
            User(l,k).power_sqrt = norm(User(l,k).beam_tx);
        end
    end 
    
    Num_Serving_User = zeros(Num_Cell,Total_BS);
    for l = 1:Num_Cell
        for k = 1:K
            iBS = mod(User(l,k).ServingCluster,Total_BS);
            if iBS == 0
                iBS = 4;
            end
            iCell = ((User(l,k).ServingCluster - iBS)/Total_BS) + 1;
            Num_Serving_User(iCell,iBS) =  + 1 ;
        end
    end

% WSR Optimization Based Clustering
elseif Mode == 2
    M = 8;
    Num_Scheduled_User = K;
    [User, Asso_User] = Neighbour(Num_Cell, K, path_loss_dB, per_BS_power_constraint_matrix,User,M,noise_dBm);
    
    for iCell = 1:Num_Cell
        [Cells(iCell).Scheduled_User, Cells(iCell).UnScheduled_User]= Round_Robin(Cells(iCell).UnScheduled_User,Num_Scheduled_User);
        if isempty(Cells(iCell).UnScheduled_User)
            Cells(iCell).UnScheduled_User = 1:K;
        end
    end
    
    Num_Serving_User = zeros(Num_Cell,Total_BS); %num of users each BS supports at the begining
    for iCell = 1:Num_Cell
        for ik = 1:length(Cells(iCell).Scheduled_User)
            k = Cells(iCell).Scheduled_User(ik);
            User(iCell,k).ServingCluster = User(iCell,k).NeighbourCluster; 
            User(iCell,k).ServingAnt = 0;
            for iBS = 1:length(User(iCell,k).ServingCluster)
                
                curr_cell = User(iCell,k).Neighbour(iBS).Cell_BS(1);
                curr_BS = User(iCell,k).Neighbour(iBS).Cell_BS(2);
                User(iCell,k).bkhaul_wgt(curr_cell,curr_BS) = 0.01;
                Num_Serving_User(curr_cell,curr_BS) = Num_Serving_User(curr_cell,curr_BS) + 1;
                if curr_BS == 1 %Macro BS
                    User(iCell,k).ServingAnt = User(iCell,k).ServingAnt + Num_TxAnt_Macro;
                else
                    User(iCell,k).ServingAnt = User(iCell,k).ServingAnt + Num_TxAnt_Pico;
                end
            end
        end
    end
    
    %Initialize wgts
    max_wgt = 0;
    for iCell = 1:Num_Cell
        for ik = 1:length(Cells(iCell).Scheduled_User)
            k = Cells(iCell).Scheduled_User(ik);
            User(iCell,k).wgt_wsr = 1/User(iCell,k).long_avg_rate;
            User(iCell,k).inst_rate = 0.01;
            User(iCell,k).prev_inst_rate = 0.01;
            if User(iCell,k).wgt_wsr > max_wgt
                max_wgt = User(iCell,k).wgt_wsr;
            end
        end
    end
    
    %normalize user wgts by max_wgt
    for iCell = 1:Num_Cell
        for ik = 1:length(Cells(iCell).Scheduled_User)
            k = Cells(iCell).Scheduled_User(ik);
            User(iCell,k).wgt_wsr = User(iCell,k).wgt_wsr/max_wgt;
            
        end
    end
    
    % initialize TX beamformers and sparse Matrix
    for l = 1:Num_Cell
        for kk = 1:length(Cells(l).Scheduled_User)
            k = Cells(l).Scheduled_User(kk); %scheduled user in cell l
            Num_ServingAnt = User(l,k).ServingAnt;
            temp_tx = randn(Num_ServingAnt,1) + 1i*randn(Num_ServingAnt,1);
            
            ant_head = 1;
            temp_tx2 = [];
            for m = 1:length(User(l,k).ServingCluster)
                
                curr_cell = ceil(User(l,k).ServingCluster(m)/Total_BS);
                curr_BS = User(l,k).ServingCluster(m) - (curr_cell-1)*Total_BS;
                if curr_BS == 1
                    Num_Tx_Ant = Num_TxAnt_Macro;
                else
                    Num_Tx_Ant = Num_TxAnt_Pico;
                end
                ant_tail = ant_head + Num_Tx_Ant - 1;
                norm_tx = temp_tx(ant_head:ant_tail)/norm(temp_tx(ant_head:ant_tail));
                temp_tx2 = [temp_tx2;sqrt(per_BS_power_constraint_matrix(curr_cell,curr_BS)/Num_Serving_User(curr_cell,curr_BS))*norm_tx];
                ant_head = ant_tail+1;
                User(l,k).BS_Power(m) = per_BS_power_constraint_matrix(curr_cell,curr_BS)/Num_Serving_User(curr_cell,curr_BS);
            end

            User(l,k).beam_tx = temp_tx2;
            if length(User(l,k).beam_tx) ~= User(l,k).ServingAnt
                fprintf('$$$ Tx Initialization Not Right $$$ \n');
            end
            User(l,k).SparseMatrix = eye(User(l,k).ServingAnt);
        end
    end
    
    %check per BS power constraint
    
    [~, ~, User] = Sum_Power(Num_Cell,Cells,User,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
    
    [User, ~, ~, ~, Cells] = WSR_Solver_Zeronorm(Num_Cell,K,Cells,User,P_max_total_dBm,Chn,noise_dBm,snr_gap_dB,intf_2nd,...
        L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,2,per_BS_power_constraint,per_BS_power_constraint_matrix,per_BS_bkhaul_constraint,iter,tau_zero,lambda); 
    
    PowerList = struct([]);
    for iCell = 1:3
        for k = 1:K
            ant_head = 1;
            power = [];
            for index = 1:length(User(iCell,k).ServingCluster)
                clust = User(iCell,k).ServingCluster(index);
                if mod(clust,Total_BS) == 1
                    ant_tail = ant_head + Num_TxAnt_Macro - 1;
                else
                    ant_tail = ant_head + Num_TxAnt_Pico - 1;
                end
                beam = User(iCell,k).beam_tx(ant_head:ant_tail);
                power = [power  beam'*beam];
                ant_head = ant_tail + 1;
            end
            PowerList(iCell,k).power = power;
        end
    end

    % Extract the implicit clustering
    for l = 1:Num_Cell
        for k = 1:K
            if any(Cells(l).Scheduled_User == k)
                % Associate User(l,k) to its serving eRRH r
                [~,Perm] = sort(PowerList(l,k).power,'descend');
                ServingBS = User(l,k).ServingCluster(Perm(1));
                ant_head = 1;
                for iCell = 1:Perm(1)-1
                    if mod(User(l,k).ServingCluster(iCell),Total_BS) == 1
                        ant_head = ant_head + Num_TxAnt_Macro;
                    else
                        ant_head = ant_head + Num_TxAnt_Pico;
                    end
                end
                if mod(ServingBS,Total_BS) == 1
                    ant_tail = ant_head + Num_TxAnt_Macro - 1;
                else
                    ant_tail = ant_head + Num_TxAnt_Pico - 1;
                end
                User(l,k).beam_tx = User(l,k).beam_tx(ant_head:ant_tail);
                User(l,k).ServingCluster = ServingBS;
                User(l,k).ServingAnt = ant_tail - ant_head + 1;

            else
                User(l,k).ServingCluster = [];
                User(l,k).ServingAnt = 0;
                User(l,k).inst_rate = 0;
                User(l,k).beam_tx = [];
            end
                
        end
    end
    
    [User,~, ~, ~, Cells] = WSR_Solver_Zeronorm_Recomputation(Num_Cell,K,Cells,User,P_max_total_dBm,Chn,noise_dBm,snr_gap_dB,...
        intf_2nd,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,2,per_BS_power_constraint,per_BS_power_constraint_matrix,per_BS_bkhaul_constraint,iter,tau_zero,lambda);
    
    for l = 1:Num_Cell
        for ik = 1:length(Cells(l).Scheduled_User)
            k = Cells(l).Scheduled_User(ik);
            User(l,k).power_sqrt = norm(User(l,k).beam_tx);
        end
    end
    
    % Update the number of Users scheduled at each eRRH
    Num_Serving_User = zeros(Num_Cell,Total_BS);
    for l = 1:Num_Cell
        for ik = 1:length(Cells(l).Scheduled_User)
            k = Cells(l).Scheduled_User(ik) ;
            iBS = mod(User(l,k).ServingCluster,Total_BS);
            if iBS == 0
                iBS = 4;
            end
            iCell = ((User(l,k).ServingCluster - iBS)/Total_BS) + 1;
            Num_Serving_User(iCell,iBS) = Num_Serving_User(iCell,iBS) + 1;         
        end
    end
    
    [~,sumrate,~] = rate_DL_iCSI(Num_Cell,Num_Rx_Ant,User,Cells,Chn,noise_dBm,intf_2nd,L_Macro,L_Pico,Num_TxAnt_Pico);
    Log_Data(iK,posNumber,iter+1).Clustering_sumrate = sumrate; 
    
    for l = 1:Num_Cell
        for k = 1:K
            Log_Data(iK,posNumber,iter+1).Clustering_User(l,k).inst_rate = User(l,k).inst_rate ;
        end
    end
    if iter == Num_Iter
        for l = 1:Num_Cell
            for k = 1:K
                s = 0;
                for ite = 0:Num_Iter
                    s = s + Log_Data(iK,posNumber,ite+1).Clustering_User(l,k).inst_rate;
                end
                Log_Data(iK,posNumber,iter+1).Clustering_User(l,k).long_avg_rate = s/(Num_Iter+1);
            end
        end
    end
    
    for sigma = 1:length(sigmaList)
        [User2,sumrate_iCSI,~] = rate_DL_iCSI(Num_Cell,Num_Rx_Ant,User,Cells,S(sigma).iCSIBank(posNumber,iter+1).channel,noise_dBm,intf_2nd,L_Macro,L_Pico,Num_TxAnt_Pico);
        Log_Data(iK,posNumber,iter+1).Clustering_sumrate_iCSI(sigma).sumrate = sumrate_iCSI;
        Log_Data(iK,posNumber,iter+1).Clustering_sumrate_iCSI(sigma).sigma = sigmaList(sigma);
        for l = 1:Num_Cell
            for k = 1:K
                Log_Data(iK,posNumber,iter+1).Clustering_sumrate_iCSI(sigma).User(l,k).inst_rate = User2(l,k).inst_rate ;
            end
        end
        if iter == Num_Iter
            for l = 1:Num_Cell
                for k = 1:K
                    s = 0;
                    for ite = 0:Num_Iter
                        s = s + Log_Data(iK,posNumber,ite+1).Clustering_sumrate_iCSI(sigma).User(l,k).inst_rate;
                    end
                    Log_Data(iK,posNumber,iter+1).Clustering_sumrate_iCSI(sigma).User(l,k).long_avg_rate = s/(Num_Iter+1);
                end
            end
        end
    end
    
% SLNR Optimization Based Clustering
elseif Mode == 3
    
    
end

