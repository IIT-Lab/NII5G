function [User, sumrate,Cells] = ...
    rate_DL_iCSI(L,Q,User,Cells,Chn,noise_dB,intf_2nd,L_Macro,L_Pico,Num_TxAnt_Pico)

%%% -----------------------------------------------------------------------
%%% This function computes rate for given transmit beamformers
%%% -----------------------------------------------------------------------
ant_gain = 1;
gap = 1;
sumrate = 0;
Total_BS = L_Macro+L_Pico;
for l = 1:L
    for ik = 1:length(Cells(l).Scheduled_User)
        k = Cells(l).Scheduled_User(ik);
        Kz_DL = eye(Q)* 10^(noise_dB/10) + intf_2nd;
        for m = 1:L
            for ikk = 1:length(Cells(m).Scheduled_User)
                kk = Cells(m).Scheduled_User(ikk);
                beam = User(m,kk).beam_tx;  %!!!!!!!!!!!!!!!!! ikk !!!!!!!!
                H_k = Get_Chn(Chn,User(m,kk).ServingCluster,l,k,L_Macro,L_Pico,Num_TxAnt_Pico);
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
        try
        invKz = inv(Kz_DL);
        
        sinr = real(H'*invKz*H)*ant_gain;
        User(l,k).inst_rate = log2(1 + sinr/gap);
        
        sumrate = sumrate + User(l,k).inst_rate;
        end
    end
end