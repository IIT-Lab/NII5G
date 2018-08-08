clearvars -except ChannelBank iCSIBank dist Num_User Num_BS BSLoc MULoc Log_Data;
%clear;
clc;

%% System Parameter
Num_Cell = 3;
L_Macro = 1;
L_Pico = 3;
ISD = 0.8;%inter sector distance in km
Total_BS = L_Macro + L_Pico;
Num_TxAnt_Macro = 4;
Num_TxAnt_Pico = 2;
Num_Rx_Ant = 2;
K = 30; %number of users per cell
BW = 10e6;%bandwidth
P_max_Macro_dBm = 43 - dec2dB(BW); %dBm/Hz
P_max_Pico_dBm = 30 - dec2dB(BW); %dBm/Hz
P_max_total_dBm = dec2dB((dB2dec(P_max_Macro_dBm)*L_Macro + dB2dec(P_max_Pico_dBm)*L_Pico)*Num_Cell); %dBm/Hz
noise_dBm = -169; %dBm/Hz
snr_gap_dB = 0;
intf_2nd = zeros(Num_Rx_Ant,Num_Rx_Ant);    %2nd tier interference
epsilon = 10;    % Interference temperature threshold (dBm)
tau_zero = 10^(-1);
lambda = 10^(-6/50); % Zero norm approximation parameters

macro_const = [26.5075   68.3125  118.5453  181.3552  268];
pico_const = [7.1892   10.6471   17.1229   29.4793  73];

backhaul_index = 5;
Macro_Bkhaul = macro_const(backhaul_index);%bps/Hz
Pico_Bkhaul = pico_const(backhaul_index);

per_BS_bkhaul_constraint = zeros(Num_Cell,Total_BS);
per_BS_power_constraint_matrix = zeros(Num_Cell,Total_BS);
for l = 1:Num_Cell
    for iBS = 1:Total_BS
        if iBS == 1 %Macro BS
            per_BS_power_constraint_matrix(l,iBS) = dB2dec(P_max_Macro_dBm);
            per_BS_bkhaul_constraint(l,iBS) = Macro_Bkhaul;
        else
            per_BS_power_constraint_matrix(l,iBS) = dB2dec(P_max_Pico_dBm);
            per_BS_bkhaul_constraint(l,iBS) = Pico_Bkhaul;
        end
    end
end
per_BS_power_constraint = reshape(per_BS_power_constraint_matrix',[],1);
%per BS power constraint in dBm/Hz

%% Load Channels
%load('ChannelBank_30UEs_3Cells_20iteUpdate.mat');
%load('iCSI_30UEs_3Cells.mat')

%% Simulation parameter
Clustering_Mode = 2;
Beamforming_Mode = 2;
Power_Mode = 2;
%1:sum power constraint ;
%2: per BS power constraint using subgraident method to find optimal dual
%variable
%21: per BS power constraint using gradient descent method to find optimal dual
%variable
%3: two total power constraint, one for Macro BS's, the other for Pico BS's

Cells = struct([]);
User = struct([]);
Res = struct([]);
%SumBackhaulConstraint = 80;

%%
%initialization
for iCell = 1:Num_Cell
    for k = 1:K
        User(iCell,k).long_avg_rate = 0.1;
        User(iCell,k).avg_rate = 0;
    end
end

T = 20; % Iterations between shadowing changes
T2 = 10; % Iterations between Dynamic Scheduling updates

%%
for posNumber = 1:6
    fprintf('Position Distribution Number : %d\n',posNumber);
    for iter = 0:29
        fprintf('Iteration Number : %d\n',iter+1);
        Chn = ChannelBank(posNumber,iter+1).CSI;
        path_loss_dB = ChannelBank(posNumber,iter+1).PL;
        for iCell = 1:Num_Cell
            Cells(iCell).Scheduled_User = [];
            Cells(iCell).UnScheduled_User = [];
            for k = 1:K
                if isempty(Log_Data(posNumber,iter+1).User(iCell,k).Cluster)
                    Cells(iCell).UnScheduled_User = [Cells(iCell).UnScheduled_User k];
                else
                    Cells(iCell).Scheduled_User = [Cells(iCell).Scheduled_User k];
                end
            end
        end
        
        Num_Serving_User = zeros(Num_Cell,Total_BS); %num of users each BS supports at the begining
        for iCell = 1:Num_Cell
            for ik = 1:length(Cells(iCell).Scheduled_User)
                k = Cells(iCell).Scheduled_User(ik);
                User(iCell,k).ServingCluster = Log_Data(posNumber,iter+1).User(iCell,k).Cluster;
                User(iCell,k).ServingAnt = length(Log_Data(posNumber,iter+1).User(iCell,k).beam_tx);
                curr_BS = mod(User(iCell,k).ServingCluster,Total_BS);
                if curr_BS == 0
                    curr_BS = 4;
                end
                curr_cell = ((User(iCell,k).ServingCluster - curr_BS)/Total_BS) + 1;
                User(iCell,k).bkhaul_wgt(curr_cell,curr_BS) = 0.01;
                Num_Serving_User(curr_cell,curr_BS) = Num_Serving_User(curr_cell,curr_BS) + 1;
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
        
        Res(posNumber,iter+1).Clustering_sumrate = sumrate;
        for sigma = 1:15
            [User, sumrate, Cells] = rate_DL_iCSI(Num_Cell,Num_Rx_Ant,User,Cells,S(sigma).iCSIBank(posNumber,iter+1).channel,noise_dBm,intf_2nd,L_Macro,L_Pico,Num_TxAnt_Pico);
            Res(posNumber,iter+1).iCSI(sigma) = sumrate;
        end
        save Log_Data_proposedTest.mat Res User;
    end
end