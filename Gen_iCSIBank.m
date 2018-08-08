iCSIBank = struct([]);
sigma = sqrt(1e-1);

for pos = 1:50
    for run = 1:50
        path_loss_dB = ChannelBank(pos,run).PL ;
        iCSIBank(pos,run).channel = GenNew_iCSI(Num_Cell, 1, Num_User, Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,path_loss_dB,ChannelBank(pos,run).CSI,sigma);
    end
end
save ('iCSI_(SSF)_30UEs_3Cells_sigma_sqrt_1_e-1.mat', 'iCSIBank', '-v7.3');