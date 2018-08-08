function [User,sumrate, wgt_sumrate,Per_BS_Bkhaul, Cells] = ...
    WSR_Solver_Zeronorm_Recomputation(Num_Cell,K,Cells,User,P_max_total_dBm,Chn,noise_dB,snr_gap_dB,...
    intf_2nd,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,Power_Mode,per_BS_power_constraint,per_BS_power_constraint_matrix,per_BS_bkhaul_constraint,iter,tau_zero,lambda)


Total_BS = L_Macro+L_Pico;
previous_rate = 1;
current_rate = 2;
max_iteration = 100;
iteration = 1;

% Initialize Zero Norm approximation weights
for l = 1:Num_Cell
    for ik = 1:length(Cells(l).Scheduled_User)
        k = Cells(l).Scheduled_User(ik);
        User(l,k).beta = 1/tau_zero;
    end
end
%%
fprintf('$$$ Solve WSR via WMMSE $$$ \n')
%while (abs(current_rate - previous_rate)/current_rate) > 1e-2 || bkhaul_flag || iteration < 30
while (abs(current_rate - previous_rate)/current_rate) > 1e-3 ||  iteration <= 30
    tic;
   
    previous_rate = current_rate;
    fprintf('\n $$$$$$$$$ WMMSE iteration= %d   $$$$$$$$$$ \n', iteration);
    
    %fprintf('$$$$$$$$$  Update Rx BF and Weight   $$$$$$$$$$ \n');
    [User] = RXBF_Wgt(Num_Cell,Cells,Num_Rx_Ant,User,Chn,noise_dB,snr_gap_dB,intf_2nd,L_Macro,L_Pico,Num_TxAnt_Pico);   
    
    %fprintf('$$$$$$$$$  Compute CVX parameters   $$$$$$$$$$ \n');
    [User, Total_Ant, Kz, b,Per_BS_Coop,Per_BS_Coop_Bkhaul,Per_BS_ServAnt] = ...
        Compute_Kz_UL(Num_Cell,User,Cells,Num_TxAnt_Macro,Num_TxAnt_Pico,L_Macro,L_Pico,Chn);
    
    N_k = max(Total_BS - floor(iteration/10),2);

    [User, ~] = CVX_QCQP_TxBF_Zero(User,Cells,N_k,Num_Cell,Power_Mode,Chn,L_Macro,L_Pico,...
        Num_TxAnt_Pico,Num_TxAnt_Macro,per_BS_power_constraint_matrix,Total_BS,per_BS_bkhaul_constraint,...
        Total_Ant, Kz, b,Per_BS_Coop,Per_BS_Coop_Bkhaul,Per_BS_ServAnt,iter,tau_zero,lambda);
    
    %check power constraint
    [~,~,User,~,~] = Sum_Power(Num_Cell,Cells,User,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico);
    
    [User, sumrate, wgt_sumrate, Per_BS_Bkhaul, curr_user_num, Cells, ~] = ...
        rate_DL(Num_Cell,Num_Rx_Ant,User,Cells,Chn,noise_dB,intf_2nd,L_Macro,L_Pico,Num_TxAnt_Pico);
    
    current_rate = wgt_sumrate;
    
    fprintf('$$$ Current Survived User = %d  $$$ \n',curr_user_num);
    fprintf('$$$$$$$$$  Updated Sum Rate= %f   $$$$$$$$$$ \n', sumrate);
    
    iteration = iteration + 1;
    
    if iteration > max_iteration
        break;
    end
   
end