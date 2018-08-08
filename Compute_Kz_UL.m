function [User, Total_Ant, Kz, b,Per_BS_Coop,Per_BS_Coop_Bkhaul,Per_BS_ServAnt]= Compute_Kz_UL(Num_Cell,User,Cells, Num_TxAnt_Macro,Num_TxAnt_Pico,L_Macro,L_Pico,Chn)
%%
Total_BS = L_Macro+L_Pico;
ant_gain = 1;
Total_Ant = 0;

%compute parameters for cvx solver to compute Tx BF
Kz = [];
b = [];
for m = 1:Num_Cell
    for l = 1:Total_BS
        Per_BS_Coop(m,l).Coop_Matrix = [];
        Per_BS_Coop_Bkhaul(m,l).Coop_Matrix = [];
        Per_BS_ServAnt(m,l).index = [];
    end
end

for l = 1:Num_Cell
    for ik = 1:length(Cells(l).Scheduled_User)
        k = Cells(l).Scheduled_User(ik);

        ant_head = 1;
%         Pico_Coop_Matrix = zeros(User(l,k).ServingAnt,User(l,k).ServingAnt);
%         Macro_Coop_Matrix = zeros(User(l,k).ServingAnt,User(l,k).ServingAnt);
%         User(l,k).BS_Coop = zeros(User(l,k).ServingAnt,User(l,k).ServingAnt,Num_Cell,Total_BS);
%         User(l,k).BS_Coop_bkhaul = zeros(User(l,k).ServingAnt,User(l,k).ServingAnt,Num_Cell,Total_BS);
        for iBS = 1:length(User(l,k).ServingCluster)
           % temp = zeros(User(l,k).ServingAnt,1);
           % temp2 = zeros(User(l,k).ServingAnt,1);
            curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
            curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
            
            if curr_BS == 1 %Macro BS
              %  temp(ant_head:(ant_head+Num_TxAnt_Macro-1)) = ones(Num_TxAnt_Macro,1);
              %  temp2(ant_head:(ant_head+Num_TxAnt_Macro-1)) = ones(Num_TxAnt_Macro,1)*sqrt(User(l,k).prev_inst_rate*User(l,k).bkhaul_wgt(curr_cell,curr_BS));
                temp2 = ones(Num_TxAnt_Macro,1)*sqrt(User(l,k).prev_inst_rate*User(l,k).bkhaul_wgt(curr_cell,curr_BS));
                Per_BS_ServAnt(curr_cell,curr_BS).index = [Per_BS_ServAnt(curr_cell,curr_BS).index (Total_Ant+ant_head):(Total_Ant+ant_head+Num_TxAnt_Macro-1)];
                Per_BS_Coop(curr_cell,curr_BS).Coop_Matrix = blkdiag(Per_BS_Coop(curr_cell,curr_BS).Coop_Matrix, eye(Num_TxAnt_Macro));
                Per_BS_Coop_Bkhaul(curr_cell,curr_BS).Coop_Matrix = blkdiag(Per_BS_Coop_Bkhaul(curr_cell,curr_BS).Coop_Matrix, diag(temp2));
                
                ant_head = ant_head + Num_TxAnt_Macro;
              %  Macro_Coop_Matrix = Macro_Coop_Matrix + diag(temp);
            else    %Pico BS
              %  temp(ant_head:(ant_head+Num_TxAnt_Pico-1)) = ones(Num_TxAnt_Pico,1);
              %  temp2(ant_head:(ant_head+Num_TxAnt_Pico-1)) = ones(Num_TxAnt_Pico,1)*sqrt(User(l,k).prev_inst_rate*User(l,k).bkhaul_wgt(curr_cell,curr_BS));
                temp2 = ones(Num_TxAnt_Pico,1)*sqrt(User(l,k).prev_inst_rate*User(l,k).bkhaul_wgt(curr_cell,curr_BS));
                Per_BS_ServAnt(curr_cell,curr_BS).index = [Per_BS_ServAnt(curr_cell,curr_BS).index (Total_Ant+ant_head):(Total_Ant+ant_head+Num_TxAnt_Pico-1)];
                Per_BS_Coop(curr_cell,curr_BS).Coop_Matrix = blkdiag(Per_BS_Coop(curr_cell,curr_BS).Coop_Matrix, eye(Num_TxAnt_Pico));
                Per_BS_Coop_Bkhaul(curr_cell,curr_BS).Coop_Matrix = blkdiag(Per_BS_Coop_Bkhaul(curr_cell,curr_BS).Coop_Matrix, diag(temp2));
                 
                ant_head = ant_head + Num_TxAnt_Pico;
              %  Pico_Coop_Matrix = Pico_Coop_Matrix + diag(temp);
            end
%             User(l,k).ServeBS(iBS).CoopMatrix = diag(temp);
%             User(l,k).BS_Coop(:,:,curr_cell,curr_BS) = diag(temp);
%             User(l,k).BS_Coop_bkhaul(:,:,curr_cell,curr_BS) = diag(temp2);
        end
        
%         User(l,k).Pico_Coop_Matrix = Pico_Coop_Matrix;
%         User(l,k).Macro_Coop_Matrix = Macro_Coop_Matrix;

        Kz_UL_intf = zeros(User(l,k).ServingAnt,User(l,k).ServingAnt);
        for m = 1:Num_Cell
            for ikk = 1:length(Cells(m).Scheduled_User)
                kk = Cells(m).Scheduled_User(ikk);

                    h = (User(m,kk).beam_rx)' * Get_Chn(Chn, User(l,k).ServingCluster, m, kk, L_Macro,L_Pico,Num_TxAnt_Pico)';
                    Kz_UL_intf = Kz_UL_intf + User(m,kk).wgt_mse*(h'*h)*ant_gain;   
            end  
        end
        
        User(l,k).Kz_UL_intf = Kz_UL_intf;
        
        Kz = blkdiag(Kz,sqrtm(Kz_UL_intf)); %matrix in the obj
        b_k = pinv(sqrtm(Kz_UL_intf))*User(l,k).wgt_mse*Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)*User(l,k).beam_rx;
        b = [b;b_k];    %vector in the obj
        
        Total_Ant = Total_Ant + User(l,k).ServingAnt;
        
    end 
end