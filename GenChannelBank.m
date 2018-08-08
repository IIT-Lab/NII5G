K = 10;
T = 20;
ChannelBank = struct([]);
Num_Cell = 3;
Num_User = K;

for pos = 1:50
    % Initialize Users' Distribution
    [Chn, path_loss_dB, dist, MULoc] = GenHetNetChn7Cell(7,1,K,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant, ISD);
    if Num_Cell == 3
        Chn = GenHetNetChn(dist,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant);
    end
    
    % Generate CSIs w/ Positions Remaining Fixed
    for run = 1:50
        if mod(run-1,T) == 0
            Chn = GenHetNetChn(dist,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant);
        else
            Chn = GenNewCSI(Num_Cell,1,K,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,path_loss_dB);
        end
        ChannelBank(pos,run).CSI = Chn;
        ChannelBank(pos,run).PL = path_loss_dB(1:Num_Cell,:,:,:);
    end
end
save ('ChannelBank_10UEs_3Cells_20iteUpdate.mat', 'ChannelBank', 'BSLoc','MULoc','dist','Num_BS','Num_User','-v7.3');