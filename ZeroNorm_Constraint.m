function [Norm2_k] = ZeroNorm_Constraint(l,ik,User,Cells,Total_BS,Num_TxAnt_Macro,Num_TxAnt_Pico,V)

k = Cells(l).Scheduled_User(ik);
clust = User(l,k).ServingCluster;
ant_head=1;
serv_head = 1;
for ikk = 1:ik-1
    kk = Cells(l).Scheduled_User(ikk);
    serv_head = serv_head + User(l,kk).ServingAnt;
end
serv_tail = serv_head + User(l,k).ServingAnt;
beam = V(serv_head:serv_tail);

Norm2_k = 0;
for iBS = 1:length(clust)
    BS = clust(iBS);
    if mod(BS,Total_BS) == 1
        ant_tail = ant_head + Num_TxAnt_Macro;
    else
        ant_tail = ant_head + Num_TxAnt_Pico;
    end
    Norm2_k = Norm2_k + User(l,k).beta * sum_square_abs(beam(ant_head:ant_tail-1)) ;
    ant_head = ant_tail;
end