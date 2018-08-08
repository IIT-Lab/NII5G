function [User, max_diff] = CVX_QCQP_TxBF(User,Cells,N_k,Num_Cell,Power_Mode,Chn,L_Macro,L_Pico,...
    Num_TxAnt_Pico,Num_TxAnt_Macro,per_BS_power_constraint_matrix,Total_BS,per_BS_bkhaul_constraint,...
    Total_Ant, Kz, b,Per_BS_Coop,Per_BS_Coop_Bkhaul,Per_BS_ServAnt)

%%%% This function finds the optimal Tx BF using cvx conditioned on
%%%% different power constraints and sparse constraints


% K = 0; %number of scheduled users
% for l=1:Num_Cell
%     K = K + length(Cells(l).Scheduled_User);
% end

%N = User(1,Cells(1).Scheduled_User(1)).ServingAnt;

%compute parameters
% Kz = [];
% b = [];
% Macro_Coop = [];
% Pico_Coop = [];
% for m = 1:Num_Cell
%     for l = 1:Total_BS
%         Per_BS_Coop(m,l).Coop_Matrix = [];
%         Per_BS_Coop_Bkhaul(m,l).Coop_Matrix = [];
%     end
% end

%N=0;
% for l = 1:Num_Cell
%     for ik = 1:length(Cells(l).Scheduled_User)
% %         k = Cells(l).Scheduled_User(ik);
% %         N = N + User(l,k).ServingAnt;
% %         Kz_k = User(l,k).Kz_UL_intf;
% %         Kz = blkdiag(Kz,sqrtm(Kz_k)); %matrix in the obj
% %
% %         b_k = pinv(sqrtm(Kz_k))*User(l,k).wgt_mse*Get_Chn(Chn, User(l,k).ServingCluster, l, k, L_Macro,L_Pico,Num_TxAnt_Pico)*User(l,k).beam_rx;
% %         b = [b;b_k];    %vector in the obj
% %
% %         Macro_Coop = blkdiag(Macro_Coop,User(l,k).Macro_Coop_Matrix);
% %         Pico_Coop = blkdiag(Pico_Coop,User(l,k).Pico_Coop_Matrix);
%
%         for c = 1:Num_Cell
%             for m = 1:Total_BS
%                 Per_BS_Coop(c,m).Coop_Matrix = blkdiag(Per_BS_Coop(c,m).Coop_Matrix,User(l,k).BS_Coop(:,:,c,m));
%                 Per_BS_Coop_Bkhaul(c,m).Coop_Matrix = blkdiag(Per_BS_Coop_Bkhaul(c,m).Coop_Matrix,User(l,k).BS_Coop_bkhaul(:,:,c,m));
%             end
%         end
%     end
% end

N = Total_Ant;

%fprintf( 'call cvx solver \n');
cvx_begin quiet
variable V(N,1) complex; %Tx BF

minimize norm( Kz*V - b)
subject to

%norm(Macro_Coop*V) <= sqrt(per_BS_power_constraint_matrix(1,1));  %Macro Total power constraint

%norm(Pico_Coop*V) <= sqrt(sum(per_BS_power_constraint_matrix(1,2:end)));  %Pico Total power constraint
for m = 1:Num_Cell
    for l = 1:Total_BS
        
        %         if size(Per_BS_Coop(m,l).Coop_Matrix,2) ~= size(Per_BS_ServAnt(m,l).index,2)
        %             size(Per_BS_Coop(m,l).Coop_Matrix,2)
        %             size(Per_BS_ServAnt(m,l).index,2)
        %         end
        
        if ~isempty(Per_BS_Coop(m,l).Coop_Matrix) && ~isempty(Per_BS_ServAnt(m,l).index)
            %             size(Per_BS_Coop(m,l).Coop_Matrix,2)
            %             size(Per_BS_ServAnt(m,l).index,2)
            
            norm(Per_BS_Coop(m,l).Coop_Matrix * V(Per_BS_ServAnt(m,l).index)) <= sqrt(per_BS_power_constraint_matrix(m,l)); %per BS power constraint
            
            norm(Per_BS_Coop_Bkhaul(m,l).Coop_Matrix * V(Per_BS_ServAnt(m,l).index)) <= sqrt(per_BS_bkhaul_constraint(m,l)); %per BS backhaul constraint
            
        end
        
    end
end

% ant_head = 1;
% for l = 1:Num_Cell
%     for ik = 1:length(Cells(l).Scheduled_User)
%         k = Cells(l).Scheduled_User(ik);
%         ant_tail = ant_head + User(l,k).ServingAnt -1;
%         if length(User(l,k).ServingCluster) > 2
%             norm(sqrtm(User(l,k).SparseMatrix)*V(ant_head : ant_tail)) <= sqrt(N_k);     %sparse constraint
%         end
%
%         ant_head = ant_tail + 1;
%
%     end
% end

cvx_end

cvx_status

%assign cvx solution to beam_tx and find the max beam_tx difference for convergence check
max_diff = 0;
ant_head = 1;
for l = 1:Num_Cell
    for ik = 1:length(Cells(l).Scheduled_User)
        k = Cells(l).Scheduled_User(ik);
        previous_beam = User(l,k).beam_tx;
        
        ant_tail = ant_head + User(l,k).ServingAnt -1;
        
        User(l,k).beam_tx = V(ant_head:ant_tail);
        diff = norm(User(l,k).beam_tx - previous_beam)/norm(previous_beam);
        
        if diff > max_diff
            max_diff = diff;
        end
        
        ant_head = ant_tail+1;
        
    end
end


%update sparse matrix according to l2 norm approximation

% delta = 1e-10;
%
% for l = 1:Num_Cell
%     for ik = 1:length(Cells(l).Scheduled_User)
%         k = Cells(l).Scheduled_User(ik);
%         SparseMatrix = [];
%
%         ant_head = 1;
%         for iBS = 1:length(User(l,k).ServingCluster)
%
%             curr_cell = ceil(User(l,k).ServingCluster(iBS)/Total_BS);
%             curr_BS = User(l,k).ServingCluster(iBS) - (curr_cell-1)*Total_BS;
%             if curr_BS == 1 %Macro BS
%                 ant_tail = ant_head + Num_TxAnt_Macro - 1;
%                 WgtMatrix = 1/(norm(User(l,k).beam_tx(ant_head:ant_tail))^2 + delta)*eye(Num_TxAnt_Macro);
%                 SparseMatrix = blkdiag(SparseMatrix,WgtMatrix);
%                 ant_head = ant_tail+1;
%             else    %Pico BS
%                 ant_tail = ant_head + Num_TxAnt_Pico - 1;
%                 WgtMatrix = 1/(norm(User(l,k).beam_tx(ant_head:ant_tail))^2 + delta)*eye(Num_TxAnt_Pico);
%                 SparseMatrix = blkdiag(SparseMatrix,WgtMatrix);
%                 ant_head = ant_tail+1;
%             end
%
%         end
%
%         User(l,k).SparseMatrix = SparseMatrix;
%
%     end
% end