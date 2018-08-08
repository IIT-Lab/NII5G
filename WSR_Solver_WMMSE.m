function [User,sumrate, wgt_sumrate,Per_BS_Bkhaul, Cells] = ...
    WSR_Solver_WMMSE(Num_Cell,K,Cells,User,P_max_total_dBm,Chn,noise_dB,snr_gap_dB,...
    intf_2nd,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,Power_Mode,per_BS_power_constraint,per_BS_power_constraint_matrix,per_BS_bkhaul_constraint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function solves WSR problem using WMMSE approach with
%%% 1. SumPower constraint in P_max_total_dBm (Power_Mode = 1) or
%%%     per BS power constraint (Power_Mode = 2)
%%% 2. Given set of weights for WSR stored in User.wgt_wsr
%%% 3. Given scheduled users stored in Cells.Scheduled_User
%%% 4. Given clustering scheme for each scheduled user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% initialize TX beamformers
for l = 1:Num_Cell
    for kk = 1:length(Cells(l).Scheduled_User)
        k = Cells(l).Scheduled_User(kk); %scheduled user in cell l
        Num_ServingAnt = User(l,k).ServingAnt;
        temp_tx = randn(Num_ServingAnt,1) + 1i*randn(Num_ServingAnt,1);
        User(l,k).beam_tx = sqrt(per_BS_power_constraint(2)/10)*temp_tx/norm(temp_tx);
    end
end
Total_BS = L_Macro+L_Pico;
previous_rate = 1;
current_rate = 2;
max_iteration = 1000;
iteration = 1;
max_diff = 1;
bkhaul_flag = 1;
Log_Data = struct([]);

%%
fprintf('$$$ Solve WSR via WMMSE $$$ \n')
%while (abs(current_rate - previous_rate)/current_rate) > 1e-2 || bkhaul_flag || iteration < 30
while (abs(current_rate - previous_rate)/current_rate) > 1e-3 ||  iteration < 30
    % stop criterion: WSR converges or max_iteration reaches
    tic;
   
    previous_rate = current_rate;
    fprintf('\n $$$$$$$$$ WMMSE iteration= %d   $$$$$$$$$$ \n', iteration);
    
    %fprintf('$$$$$$$$$  Update Rx BF and Weight   $$$$$$$$$$ \n');
    [User] = RXBF_Wgt(Num_Cell,Cells,Num_Rx_Ant,User,Chn,noise_dB,snr_gap_dB,intf_2nd,L_Macro,L_Pico,Num_TxAnt_Pico);   
    
    %fprintf('$$$$$$$$$  Compute CVX parameters   $$$$$$$$$$ \n');
    [User, Total_Ant, Kz, b,Per_BS_Coop,Per_BS_Coop_Bkhaul,Per_BS_ServAnt] = ...
        Compute_Kz_UL(Num_Cell,User,Cells,Num_TxAnt_Macro,Num_TxAnt_Pico,L_Macro,L_Pico,Chn);
    
    N_k = max(Total_BS - floor(iteration/10),2);
    %fprintf('$$$$$$$$$  Compute Tx BF using CVX   $$$$$$$$$$ \n');
    %if Power_Mode == 3
        [User, max_diff] = CVX_QCQP_TxBF(User,Cells,N_k,Num_Cell,Power_Mode,Chn,L_Macro,L_Pico,...
            Num_TxAnt_Pico,Num_TxAnt_Macro,per_BS_power_constraint_matrix,Total_BS,per_BS_bkhaul_constraint,...
            Total_Ant, Kz, b,Per_BS_Coop,Per_BS_Coop_Bkhaul,Per_BS_ServAnt);
    %end
    
    %fprintf('$$$$ Tx BF diff = %f   $$$$ \n',max_diff);
    
    %check power constraint
    [~,per_BS_power,User,max_power_diff,max_dev_wgt] = Sum_Power(Num_Cell,Cells,User,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
    power_gap = (per_BS_power - per_BS_power_constraint_matrix)./per_BS_power_constraint_matrix;
    max_power_gap = max(power_gap);
    
    fprintf('$$$ max deviated wgt = %f $$$\n',max_dev_wgt);
    
    fprintf('$$$ max power gap = %f $$$ \n',max_power_gap);
    fprintf('\n');
    fprintf('$$$ max power diff = %f $$$\n',max_power_diff);
    % sumpower = Sum_Power(Num_Cell,Cells,User,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
    
    [User, sumrate, wgt_sumrate, Per_BS_Bkhaul, curr_user_num, Cells, max_rate_diff] = ...
        rate_DL(Num_Cell,Num_Rx_Ant,User,Cells,Chn,noise_dB,intf_2nd,L_Macro,L_Pico,Num_TxAnt_Pico);
    
    bkhaul_gap = (Per_BS_Bkhaul - per_BS_bkhaul_constraint)./per_BS_bkhaul_constraint;
    max_bkhaul_gap = max(bkhaul_gap);
    fprintf('$$$ max rate diff = %f $$$ \n',max_rate_diff);
    fprintf('$$$ max bkhaul gap = %f $$$ \n',max_bkhaul_gap);
    
    if max(max_bkhaul_gap) < 0.01
        bkhaul_flag = 0;
    else
        bkhaul_flag = 1;
    end
    
    current_rate = wgt_sumrate;
    
    fprintf('$$$ Current Survived User = %d  $$$ \n',curr_user_num);
    fprintf('$$$$$$$$$  Updated Sum Rate= %f   $$$$$$$$$$ \n', sumrate);
    
%     if (abs(current_rate - previous_rate)/current_rate) < 1e-2
%         (abs(current_rate - previous_rate)/current_rate)
%     end
    
    iteration = iteration + 1;
    
    if iteration > max_iteration
        break;
    end
%     if max_power_diff <1e-1
%         iteration
%     end
    Log_Data = [Log_Data User];
    toc;
%    save Log_Data_proposed_New_10.mat Log_Data User;
   
   
end

%fprintf('$$$$$$$$$  Final Sum Rate= %f  $$$$$$$$$$ \n', sumrate);