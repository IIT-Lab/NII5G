function [User, max_diff] = CVX_QCQP_TxBF_Zero(User,Cells,N_k,Num_Cell,Power_Mode,Chn,L_Macro,L_Pico,...
    Num_TxAnt_Pico,Num_TxAnt_Macro,per_BS_power_constraint_matrix,Total_BS,per_BS_bkhaul_constraint,...
    Total_Ant, Kz, b,Per_BS_Coop,Per_BS_Coop_Bkhaul,Per_BS_ServAnt,iter,tau_zero,lambda)

N = Total_Ant;

%fprintf( 'call cvx solver \n');
cvx_begin quiet
variable V(N,1) complex; %Tx BF

minimize norm( Kz*V - b)
subject to

for m = 1:Num_Cell
    for l = 1:Total_BS        
        if ~isempty(Per_BS_Coop(m,l).Coop_Matrix) && ~isempty(Per_BS_ServAnt(m,l).index)
            
            norm(Per_BS_Coop(m,l).Coop_Matrix * V(Per_BS_ServAnt(m,l).index)) <= sqrt(per_BS_power_constraint_matrix(m,l)); %per BS power constraint
            
            norm(Per_BS_Coop_Bkhaul(m,l).Coop_Matrix * V(Per_BS_ServAnt(m,l).index)) <= sqrt(per_BS_bkhaul_constraint(m,l)); %per BS backhaul constraint            
        end       
    end
end

for l = 1:Num_Cell
    for ik = 1:length(Cells(l).Scheduled_User)
        ZeroNorm_Constraint(l,ik,User,Cells,Total_BS,Num_TxAnt_Macro,Num_TxAnt_Pico,V) <= 1;
    end
end

cvx_end

cvx_status

%assign cvx solution to beam_tx and find the max beam_tx difference for convergence check
max_diff = 0;
ant_head = 1;
tau = tau_zero*lambda^iter;
for l = 1:Num_Cell
    for ik = 1:length(Cells(l).Scheduled_User)
        k = Cells(l).Scheduled_User(ik);
        previous_beam = User(l,k).beam_tx;
        
        ant_tail = ant_head + User(l,k).ServingAnt -1;
        User(l,k).beta = 1/(sum_square_abs(V(ant_head:ant_tail)) + tau);
        User(l,k).beam_tx = V(ant_head:ant_tail);
        diff = norm(User(l,k).beam_tx - previous_beam)/norm(previous_beam);
        
        if diff > max_diff
            max_diff = diff;
        end
        
        ant_head = ant_tail+1;
        
    end
end
