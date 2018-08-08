function [sumpower, per_BS_power, User, max_diff, max_dev_wgt] = Sum_Power(L,Cells,User,L_Macro,L_Pico,Num_TxAnt_Macro,Num_TxAnt_Pico)

sumpower = 0;
Total_BS = L_Macro + L_Pico;
per_BS_power = zeros(L,Total_BS);

delta = 1e-12;
max_diff = 0;
max_dev = 0; % max deviation of the backhaul weight compared with the real bf power 
max_dev_wgt = 0;
for l = 1:L
    for ik = 1:length(Cells(l).Scheduled_User)
        k = Cells(l).Scheduled_User(ik);
        beam = User(l,k).beam_tx;
        sumpower = sumpower + real(beam'*beam);
        
        ant_head = 1;
        null_entry = [];
        null_BS = [];
        for iBS = 1:length(User(l,k).ServingCluster)
            curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
            curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
            
            if curr_BS == 1 %Macro BS
                ant_tail = ant_head + Num_TxAnt_Macro - 1;
            else
                ant_tail = ant_head + Num_TxAnt_Pico - 1;
            end
            
            curr_BS_power = real((beam(ant_head : ant_tail))'*beam(ant_head : ant_tail));
            
            per_BS_power(curr_cell,curr_BS) = per_BS_power(curr_cell,curr_BS) + curr_BS_power;
            
            
            if abs(User(l,k).bkhaul_wgt(curr_cell,curr_BS)*curr_BS_power - 1) > max_dev
                max_dev = abs(User(l,k).bkhaul_wgt(curr_cell,curr_BS)*curr_BS_power - 1);
                max_dev_wgt = User(l,k).bkhaul_wgt(curr_cell,curr_BS)*curr_BS_power;
            end
            
            
            User(l,k).bkhaul_wgt(curr_cell,curr_BS) = 1/(curr_BS_power + delta);
                       
            if abs(curr_BS_power - User(l,k).BS_Power(iBS))/User(l,k).BS_Power(iBS) > max_diff
                max_diff = abs(curr_BS_power - User(l,k).BS_Power(iBS))/User(l,k).BS_Power(iBS);
            end
                      
            User(l,k).BS_Power(iBS) = curr_BS_power;
            
            if User(l,k).BS_Power(iBS) < 1e-10
                curr_null_entry = [ant_head:ant_tail];
                null_entry = [null_entry curr_null_entry];
                null_BS = [null_BS iBS];
            end
            ant_head = ant_tail + 1;
        end
        
        User(l,k).ServingCluster(null_BS) = [];
        User(l,k).beam_tx(null_entry) = [];
        User(l,k).ServingAnt = User(l,k).ServingAnt - length(null_entry);
        %         SparseVec2 = diag(User(l,k).SparseMatrix);
        %         SparseVec = diag(User(l,k).SparseMatrix);
        %         SparseVec(null_entry) = [];
        if isempty(User(l,k).ServingCluster)
            User(l,k).ServingCluster = curr_BS;
            User(l,k).beam_tx = beam(curr_null_entry);
            User(l,k).ServingAnt = length(User(l,k).beam_tx);
            %             SparseVec = SparseVec2(curr_null_entry);
        end
        %         User(l,k).SparseMatrix = diag(SparseVec);
        
    end
end

