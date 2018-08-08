%InitializeParameters;
L = Num_Cell;
Q = Num_Rx_Ant;
ant_gain = 1;
gap = 1;
SR_List2 = [];
for posNumber = 1:1
    for iter = 1:30
        Chn = iCSIBank(posNumber,iter).channel;
        sumrate = 0;
        Total_BS = L_Macro+L_Pico;
        for l = 1:L
            for k = 1:K
                Kz_DL = eye(Q)* 10^(noise_dBm/10) + intf_2nd;
                for m = 1:L
                    for kk = 1:K
                        %beam = User(m,kk).beam_tx;  %!!!!!!!!!!!!!!!!! ikk !!!!!!!!
                        beam = Log_Data(posNumber,iter).sumrate_iCSI(m,kk).beam_tx;
                        H_k = Get_Chn(Chn,Log_Data(posNumber,iter).sumrate_iCSI(m,kk).ServingCluster,l,k,L_Macro,L_Pico,Num_TxAnt_Pico);
                        %H_k = Get_Chn(Chn,User(m,kk).ServingCluster,l,k,L_Macro,L_Pico,Num_TxAnt_Pico);
                        %H_k = Get_Chn(Chn,User(m,ikk).Cluster,l,k,L_Macro,L_Pico,Num_TxAnt_Pico);
                        if (m ~= l) || (kk ~=k)
                            %if ~isempty(H_k)
                            Kz_DL = Kz_DL + H_k' * (beam*beam') * H_k * ant_gain;
                            %end
                        else
                            H = H_k'*beam;
                        end
                    end
                end
            end
            try
                invKz = inv(Kz_DL);
                
                sinr = real(H'*invKz*H)*ant_gain;
                Rate = log2(1 + sinr/gap);
                
                sumrate = sumrate + Rate;
            end
        end
        SR_List2 = [SR_List2 sumrate];
    end
end
