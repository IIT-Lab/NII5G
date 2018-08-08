function [User, Cells, sumrate, wgt_sumrate] = eRRH_beamforming(Num_Cell,K,Cells,User,Chn,noise_dB,snr_gap_dB,...
    intf_2nd,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,per_BS_power_constraint_matrix,Num_Serving_User,epsilon)


Total_BS = L_Pico+L_Macro;
max_iteration = 100;
ant_gain = 1;
snr_gap = dB2dec(snr_gap_dB);

inter_Temp = 0;
inter_Range = 1;

% initialize TX beamformers
for l = 1:Num_Cell
    for ik = 1:length(Cells(l).Scheduled_User)
        k = Cells(l).Scheduled_User(ik);
        Num_ServingAnt = User(l,k).ServingAnt;
        temp_tx = randn(Num_ServingAnt,1) + 1i*randn(Num_ServingAnt,1);
        User(l,k).beam_tx = sqrt(per_BS_power_constraint_matrix(1,2)/10)*temp_tx/norm(temp_tx);
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

for l = 1:Num_Cell
    for iBS = 1:Total_BS
        Num_Users_Served = Num_Serving_User(l,iBS); % Number of Users Served by this BS
        if Num_Users_Served == 0
            continue
        elseif Num_Users_Served == 2
            1;
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
        
        previous_rate = 1;
        current_rate = 2;
        iteration = 1;
        
        fprintf('$$$ Solve WSR for BS (%d,%d) $$$ \n',l,iBS)
        while (abs(current_rate - previous_rate)/current_rate) > 1e-3 || iteration <= 30    % stop criterion: WSR converges or max_iteration reaches
            tic;
            previous_rate = current_rate;
            fprintf('\n $$$$$$$$$ WMMSE iteration= %d   $$$$$$$$$$ \n', iteration); 
            for k_index = 1:Num_Users_Served
                iCell = BS_cluster(k_index,1);
                k = BS_cluster(k_index,2);    % User identification: (iCell,k)
                H = Get_Chn(Chn, User(iCell,k).ServingCluster, iCell,k, L_Macro,L_Pico,Num_TxAnt_Pico);    % Channel Matrix Between User (iCell,k) and its Serving eRRH (l,iAnt)
        
                Kz_DL = eye(Num_Rx_Ant)* 10^(noise_dB/10) + intf_2nd;    %total received signal covariance matrix for user (iCell,k)
                
                if inter_Range == 1
                    % We Consider only Interference from Users under the same eRRH (Local)
                    for kk_index = 1:Num_Users_Served
                            Cell_ID = BS_cluster(kk_index,1);
                            kk = BS_cluster(kk_index,2);
                            beam = User(Cell_ID,kk).beam_tx;
                            % Hkk = Get_Chn(Chn, User(Cell_ID,kk).ServingCluster, iCell,k, L_Macro,L_Pico,Num_TxAnt_Pico)
                            Kz_DL = Kz_DL + ant_gain*H'*(beam*beam')*H;
                    end
                
                elseif inter_Range == 2
                    % We Consider Interference from all Users (Global)
                    for Cell_ID = 1:Num_Cell
                        for kk = 1:K
                                beam = User(Cell_ID,kk).beam_tx;
                                Hkk = Get_Chn(Chn, User(Cell_ID,kk).ServingCluster, iCell,k, L_Macro,L_Pico,Num_TxAnt_Pico);
                                Kz_DL = Kz_DL + ant_gain*Hkk'*(beam*beam')*Hkk;
                        end
                    end
                    
                end
        
                beam_rx = Kz_DL\(H'*User(iCell,k).beam_tx)*sqrt(ant_gain);
        
                User(iCell,k).beam_rx = beam_rx;
        
                User(iCell,k).mse = 1 - real( (User(iCell,k).beam_rx)' * H' * User(iCell,k).beam_tx * sqrt(ant_gain) );
                %User(iCell,k).mse = User(iCell,k).beam_rx'*H'*User(iCell,k).beam_tx + 1 - 2*real( (User(iCell,k).beam_rx)' * H' * User(iCell,k).beam_tx * sqrt(ant_gain) );
                
                %wgt_mse(l,s,k,n) = 1/R_avg(l,s,k)/mse(l,s,k,n);
                lambda = 1/(1 + (1/(User(iCell,k).mse) - 1)/snr_gap );
                User(iCell,k).wgt_mse = lambda*User(iCell,k).wgt_wsr/(User(iCell,k).mse)^2/snr_gap;
       
            end

            %Compute Parameters for CVX Solver to Compute Tx BF
            Total_Ant = 0;
            Kz = [];
            b = [];
            Per_BS_Coop(l,iBS).Coop_Matrix = [];
            Per_BS_ServAnt(l,iBS).index = [];
            %BS_intf = 0;

            for k_index = 1:Num_Users_Served
                iCell = BS_cluster(k_index,1);
                k = BS_cluster(k_index,2);
                ant_head = 1;
                %BS_intf_k = 0;
                if iBS == 1 %Macro BS
                    Per_BS_ServAnt(l,iBS).index = [Per_BS_ServAnt(l,iBS).index (Total_Ant+ant_head):(Total_Ant+ant_head+Num_TxAnt_Macro-1)];
                    Per_BS_Coop(l,iBS).Coop_Matrix = blkdiag(Per_BS_Coop(l,iBS).Coop_Matrix, eye(Num_TxAnt_Macro));
                    ant_head = ant_head + Num_TxAnt_Macro;
                else    %Pico BS
                    Per_BS_ServAnt(l,iBS).index = [Per_BS_ServAnt(l,iBS).index (Total_Ant+ant_head):(Total_Ant+ant_head+Num_TxAnt_Pico-1)];
                    Per_BS_Coop(l,iBS).Coop_Matrix = blkdiag(Per_BS_Coop(l,iBS).Coop_Matrix, eye(Num_TxAnt_Pico));
                    ant_head = ant_head + Num_TxAnt_Pico;
                end
            
                Kz_UL_intf = zeros(User(iCell,k).ServingAnt,User(iCell,k).ServingAnt);
                
                if inter_Range == 1    % We Consider only the Interference from Users Served by the Same BS (Local)
                    for kk_index = 1:Num_Users_Served 
                        Cell_ID = BS_cluster(kk_index,1);
                        kk = BS_cluster(kk_index,2);
                        h = (User(Cell_ID,kk).beam_rx)' * Get_Chn(Chn, User(iCell,k).ServingCluster, Cell_ID, kk, L_Macro,L_Pico,Num_TxAnt_Pico)';
                        Kz_UL_intf = Kz_UL_intf + User(Cell_ID,kk).wgt_mse*(h'*h)*ant_gain;
                        %BS_intf_k = BS_intf_k + (User(iCell,k).beam_rx)' * Get_Chn(Chn, User(iCell,k).ServingCluster, Cell_ID, kk, L_Macro,L_Pico,Num_TxAnt_Pico)';
                    end
                    
                elseif inter_Range == 2    % We Consider all Interference from Users (Global)
                    for Cell_ID = 1:Num_Cell
                        for kk = 1:K
                            h = (User(Cell_ID,kk).beam_rx)' * Get_Chn(Chn, User(Cell_ID,kk).ServingCluster, iCell, k, L_Macro,L_Pico,Num_TxAnt_Pico)';
                            Kz_UL_intf = Kz_UL_intf + User(Cell_ID,kk).wgt_mse*(h'*h)*ant_gain;
                        end
                    end
                end
        
                User(iCell,k).Kz_UL_intf = Kz_UL_intf;
                Kz = blkdiag(Kz,sqrtm(Kz_UL_intf)); %matrix in the obj
                %BS_intf = BS_intf + BS_intf_k'*BS_intf_k*ant_gain;
                b_k = pinv(sqrtm(Kz_UL_intf))*User(iCell,k).wgt_mse*Get_Chn(Chn, User(iCell,k).ServingCluster, iCell, k, L_Macro,L_Pico,Num_TxAnt_Pico)*User(iCell,k).beam_rx;
                b = [b;b_k];    %vector in the obj
                Total_Ant = Total_Ant + User(iCell,k).ServingAnt;        
            end 
            N_k = max(Total_BS - floor(iteration/10),2);
            N = Total_Ant;

            cvx_begin quiet
            variable V(N,1) complex;

            minimize norm( Kz*V - b)
            subject to
        
            if ~isempty(Per_BS_Coop(l,iBS).Coop_Matrix) && ~isempty(Per_BS_ServAnt(l,iBS).index)
                norm(V) <= sqrt(per_BS_power_constraint_matrix(l,iBS)); %per BS power constraint
            end
            
            if inter_Temp == 1
                
            end
            %BS_intf <= dB2dec(epsilon);
            
            cvx_end

            cvx_status

            %assign cvx solution to beam_tx
            ant_head = 1;
            for k_index = Num_Users_Served
                iCell = BS_cluster(k_index,1);
                k = BS_cluster(k_index,2); 
        
                ant_tail = ant_head + User(iCell,k).ServingAnt -1;
                User(iCell,k).beam_tx = V(ant_head:ant_tail);
                ant_head = ant_tail+1;
   
            end

            gap = 1;
            sumrate = 0;
            wgt_sumrate = 0;
            for ik = 1:Num_Users_Served
                iCell = BS_cluster(ik,1);
                k = BS_cluster(ik,2);
                Kz_DL = eye(Num_Rx_Ant)* 10^(noise_dB/10) + intf_2nd;
                for ikk = 1:Num_Users_Served
                    iC = BS_cluster(ikk,1);
                    kk = BS_cluster(ikk,2);
                    beam = User(iC,kk).beam_tx;
                    H_k = Get_Chn(Chn,User(iC,kk).ServingCluster,iCell,k,L_Macro,L_Pico,Num_TxAnt_Pico);
                    if (iC ~= iCell) || (kk ~=k)
                        Kz_DL = Kz_DL + H_k' * (beam*beam') * H_k * ant_gain;
                    else
                        H = H_k'*beam;
                    end
                end
                
                invKz = inv(Kz_DL);
                sinr = real(H'*invKz*H)*ant_gain;
                User(iCell,k).inst_rate = log2(1 + sinr/gap);
                User(iCell,k).wgt_rate = User(iCell,k).wgt_wsr*User(iCell,k).inst_rate;
                
                sumrate = sumrate + User(iCell,k).inst_rate;
                wgt_sumrate = wgt_sumrate + User(iCell,k).wgt_rate;
            end
        
            current_rate = wgt_sumrate;
    
            fprintf('$$$$$$$$$  Updated Sum Rate= %f   $$$$$$$$$$ \n', sumrate);
    
            iteration = iteration + 1;
            
            if iteration > max_iteration
                if current_rate < previous_rate*0.95
                    iteration = max_iteration ;
                end
                break;
            end
            toc;
        end
    end
end
        