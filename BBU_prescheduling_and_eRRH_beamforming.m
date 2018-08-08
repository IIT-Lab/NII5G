%clearvars -except ChannelBank iCSIBank dist Num_User Num_BS BSLoc MULoc S;
clear;
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
%K = 30; %number of users per cell
K_List = [5 10 15 20 25 30];
BW = 10e6;%bandwidth
P_max_Macro_dBm = 43 - dec2dB(BW); %dBm/Hz
P_max_Pico_dBm = 30 - dec2dB(BW); %dBm/Hz
P_max_total_dBm = dec2dB((dB2dec(P_max_Macro_dBm)*L_Macro + dB2dec(P_max_Pico_dBm)*L_Pico)*Num_Cell); %dBm/Hz
noise_dBm = -169; %dBm/Hz
P_noise = dB2dec(noise_dBm) ; %(mW/Hz)
snr_gap_dB = 0;
intf_2nd = zeros(Num_Rx_Ant,Num_Rx_Ant);    %2nd tier interference
epsilon = 10;    % Interference temperature threshold (dBm)
tau_zero = 10^(-1);
lambda = 10^(-6/50); % Zero norm approximation parameters
T = 20; % Iterations between shadowing changes
T2 = 5; % Iterations between Dynamic Scheduling updates

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
%S(1) = load('iCSI_(SSF)_30UEs_3Cells_sigma_1.mat');
%S(2) = load('iCSI_(SSF)_30UEs_3Cells_sigma_sqrt_1_e-1.mat');
%S(3) = load('iCSI_(SSF)_30UEs_3Cells_sigma_1_e-1.mat');

sigmaList = [1 sqrt(1e-1) 1e-1];

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

Log_Data = struct([]);
Num_Iter = 29;
Num_Pos = 1;
%SumBackhaulConstraint = 80;

%%
%for iK = 1:length(K_List)
for iK = 6
    clearvars ChannelBank S dist Num_User Num_BS BSLoc MULoc
    K = K_List(iK) ;
    load(strcat('ChannelBank_',num2str(K),'UEs_3Cells_20iteUpdate.mat'));
    S(1) = load(strcat('iCSI_(SSF)_',num2str(K),'UEs_3Cells_sigma_1.mat'));
    S(2) = load(strcat('iCSI_(SSF)_',num2str(K),'UEs_3Cells_sigma_sqrt_1_e-1.mat'));
    S(3) = load(strcat('iCSI_(SSF)_',num2str(K),'UEs_3Cells_sigma_1_e-1.mat'));
    
    Cells = struct([]);
    User = struct([]);
    
    %initialization
    for iCell = 1:Num_Cell
        for k = 1:K
            User(iCell,k).long_avg_rate = 0.1;
            User(iCell,k).avg_rate = 0;
        end
        Cells(iCell).UnScheduled_User = 1:K;
    end
    
    Num_Scheduled_User = K;
    
    for posNumber = 1:Num_Pos
        fprintf('Position Distribution Number : %d\n',posNumber);
        for iter = 0:Num_Iter
            fprintf('Iteration Number : %d\n',iter+1);
            Chn = ChannelBank(posNumber,iter+1).CSI;
            path_loss_dB = ChannelBank(posNumber,iter+1).PL;
            iChn = S(1).iCSIBank(posNumber,iter+1).channel;
            % Pre-scheduling
            %if mod(iter,T2) == 0
            [User, Cells, Asso_User, Num_Serving_User, Log_Data] = Pool_clustering(User,Cells,Num_Cell, K, Chn, S, sigmaList, Total_BS, path_loss_dB, snr_gap_dB, intf_2nd, per_BS_power_constraint_matrix,...
                per_BS_power_constraint,per_BS_bkhaul_constraint,noise_dBm,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,L_Macro,L_Pico,P_max_total_dBm,iter,tau_zero,lambda,posNumber,Log_Data,Clustering_Mode,iK,Num_Iter);
            %end
            
            % eRRH Beamforming
            for sigma = 0:length(sigmaList)
                if sigma == 0
                    if Beamforming_Mode == 1
                        [User, Cells, sumrate, wgt_sumrate] = eRRH_beamforming(Num_Cell,K,Cells,User,Chn,noise_dBm,snr_gap_dB,...
                            intf_2nd,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,per_BS_power_constraint_matrix,Num_Serving_User,epsilon);
                    elseif Beamforming_Mode == 2
                        [User, Cells] = Direct_Beamformers(Num_Cell,K,Cells,User,Chn,noise_dBm,L_Macro,...
                            L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,per_BS_power_constraint_matrix,Num_Serving_User);
                    end
                    [User, sumrate, Cells] = rate_DL_iCSI(Num_Cell,Num_Rx_Ant,User,Cells,Chn,noise_dBm,intf_2nd,L_Macro,L_Pico,Num_TxAnt_Pico);
                    Log_Data(iK,posNumber,iter+1).pCSI_sumrate = sumrate;
                    for l = 1:Num_Cell
                        for k = 1:K
                            Log_Data(iK,posNumber,iter+1).pCSI_User(l,k).inst_rate = User(l,k).inst_rate;
                        end
                    end
                    if iter == Num_Iter
                        for l = 1:Num_Cell
                            for k = 1:K
                                s = 0;
                                for ite = 0:Num_Iter
                                    s = s + Log_Data(iK,posNumber,ite+1).pCSI_User(l,k).inst_rate;
                                end
                                Log_Data(iK,posNumber,iter+1).pCSI_User(l,k).long_avg_rate = s/(Num_Iter+1);
                            end
                        end
                    end 
                else
                    if Beamforming_Mode == 1
                        [User, Cells, sumrate, wgt_sumrate] = eRRH_beamforming(Num_Cell,K,Cells,User,S(sigma).iCSIBank(posNumber,iter+1).channel,noise_dBm,snr_gap_dB,...
                            intf_2nd,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,per_BS_power_constraint_matrix,Num_Serving_User,epsilon);
                    elseif Beamforming_Mode == 2
                        [User, Cells] = Direct_Beamformers(Num_Cell,K,Cells,User,S(sigma).iCSIBank(posNumber,iter+1).channel,noise_dBm,L_Macro,...
                            L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,per_BS_power_constraint_matrix,Num_Serving_User);
                    end
                    % fprintf('Sumrate = %d \n',sumrate);
                    [User, sumrate, Cells] = rate_DL_iCSI(Num_Cell,Num_Rx_Ant,User,Cells,S(sigma).iCSIBank(posNumber,iter+1).channel,noise_dBm,intf_2nd,L_Macro,L_Pico,Num_TxAnt_Pico);
                    Log_Data(iK,posNumber,iter+1).BF(sigma).sumrate = sumrate;
                    Log_Data(iK,posNumber,iter+1).BF(sigma).sigma = sigmaList(sigma);
                    for l = 1:Num_Cell
                        for k = 1:K
                            Log_Data(iK,posNumber,iter+1).BF(sigma).User(l,k).inst_rate = User(l,k).inst_rate;
                        end
                    end
                    if iter == Num_Iter
                        for l = 1:Num_Cell
                            for k = 1:K
                                s = 0;
                                for ite = 0:Num_Iter
                                    s = s + Log_Data(iK,posNumber,ite+1).BF(sigma).User(l,k).inst_rate;
                                end
                                Log_Data(iK,posNumber,iter+1).BF(sigma).User(l,k).long_avg_rate = s/(Num_Iter+1);
                            end
                        end
                    end

                end
                %[Log_Utility,User,Log_Data] = Compute_Utility_2(Num_Cell,K,User,Cells,iter+1,posNumber,Log_Data);
            end
            save Log_Data_proposedTest.mat Log_Data User;
        end
    end
    
end