
%initialization
for iCell = 1:Num_Cell
    for k = 1:K
        User(iCell,k).long_avg_rate = 0.1;
        User(iCell,k).avg_rate = 0;
    end
    Cells(iCell).UnScheduled_User = 1:K;
end

M = 8; %maximum number of serving BSs per User

Chn = ChannelBank(posNumber,1).CSI;
path_loss_dB = ChannelBank(posNumber,1).PL;
[User Asso_User] = Neighbour(Num_Cell, K, path_loss_dB, per_BS_power_constraint_matrix,User,M,noise_dBm);

out_iter = 1;
previous_utility = 1;
current_utility = 2;

while abs((current_utility - previous_utility)/current_utility) > 1e-3 || out_iter <=50
%while out_iter <= 1
    % each time schedule 10 users, 5 iterations is needed for all users
    % being served. we wish every user being served at least 10 times

    %if out_iter ~= 1
    %    if mod(out_iter,T) == 0
    %        Chn = GenHetNetChn(dist,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant);
    %    else
    %        Chn = GenNewCSI(Num_Cell,1,K,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant);
    %    end
    %end

    Chn = ChannelBank(posNumber,out_iter).CSI;
    path_loss_dB = ChannelBank(posNumber,out_iter).PL;
    
    previous_utility = current_utility;
    %timerTot = tic;
    
    
    fprintf('$$$$$$ User Scheduling iter = %d $$$$$$$ \n', out_iter);
    %User scheduler
    for iCell = 1:Num_Cell
        [Cells(iCell).Scheduled_User, Cells(iCell).UnScheduled_User]= Round_Robin(Cells(iCell).UnScheduled_User,Num_Scheduled_User);
       
%         Log_Data(out_iter).Cells(iCell).Scheduled_User = Cells(iCell).Scheduled_User;
%         Log_Data(out_iter).Cells(iCell).UnScheduled_User = Cells(iCell).UnScheduled_User; % log user scheduling
        if isempty(Cells(iCell).UnScheduled_User)
            Cells(iCell).UnScheduled_User = 1:K;
        end
  %      Cells(iCell).initial_lambda = zeros(1,Total_BS);
    end
    
    %initial Clustering 
    %every scheduled user connects with all the candidate BS's
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
    
    [User,sumrate, wgt_sumrate, Per_BS_Bkhaul, Cells] = WSR_Solver_WMMSE(Num_Cell,K,Cells,User,P_max_total_dBm,Chn,noise_dBm,snr_gap_dB,intf_2nd,...
        L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,Power_Mode,per_BS_power_constraint,per_BS_power_constraint_matrix,per_BS_bkhaul_constraint); 

    Log_Data(posNumber,out_iter).Per_BS_Bkhaul = Per_BS_Bkhaul;
    
    [Log_Utility,User,Log_Data] = Compute_Utility(Num_Cell,K,User,Cells,out_iter,posNumber,Log_Data);
    
    
    for sigma = 1:11
        [~,sumrate_iCSI,~] = rate_DL_iCSI(Num_Cell,Num_Rx_Ant,User,Log_Data(posNumber,out_iter).Cells,S(sigma).iCSIBank(posNumber,out_iter).channel,noise_dBm,intf_2nd,L_Macro,L_Pico,Num_TxAnt_Pico);
        Log_Data(posNumber,out_iter).iCSI(sigma).sumrate = sumrate_iCSI;
        Log_Data(posNumber,out_iter).iCSI(sigma).sigma = sigmaList(sigma);
    end
    
    
    fprintf('$$$$ Current Log Utility = %f $$$$ \n\n\n\n',Log_Utility);
    
    current_utility = Log_Utility;
    out_iter = out_iter + 1;
    
    save Log_Data_proposedTest.mat Log_Data User;
    %timerTot2 = toc;
    %fprintf('$$$$$$$$ Total Elapsed Time = %f $$$$$$$ \n\n\n', timerTot2-timerTot);
    
end