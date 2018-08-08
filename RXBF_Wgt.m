function [User] = RXBF_Wgt(Num_Cell,Cells,Q,User,Chn,noise_dB,snr_gap_dB,intf_2nd,L_Macro,L_Pico,Num_TxAnt_Pico)

%%% ----------------------------------------------------------------------
%%% This function computes
%%% 1. MMSE receive beamforming: User.beam_rx
%%% 2. MSE: User.mse
%%% 3. weight for WMMSE: User.wgt_mse
%%% for given transmit beamformers.
%%% ----------------------------------------------------------------------

ant_gain = 1;
snr_gap = dB2dec(snr_gap_dB);
for l = 1:Num_Cell
    for k_index = 1:length(Cells(l).Scheduled_User)
        k = Cells(l).Scheduled_User(k_index);
        H = Get_Chn(Chn, User(l,k).ServingCluster, l,k, L_Macro,L_Pico,Num_TxAnt_Pico);
        
        Kz_DL = eye(Q)* 10^(noise_dB/10) + intf_2nd;
        %total received signal covariance matrix for user k
        
        for m = 1:Num_Cell
            for kk_index = 1:length(Cells(m).Scheduled_User)
                kk = Cells(m).Scheduled_User(kk_index);
                beam = User(m,kk).beam_tx;
                Hkk = Get_Chn(Chn, User(m,kk).ServingCluster, l,k, L_Macro,L_Pico,Num_TxAnt_Pico);
                Kz_DL = Kz_DL + ant_gain*Hkk'*(beam*beam')*Hkk;
                
            end
        end
        
        beam_rx = Kz_DL\(H'*User(l,k).beam_tx)*sqrt(ant_gain);
        
        User(l,k).beam_rx = beam_rx;
        
        User(l,k).mse = 1 - real( (User(l,k).beam_rx)' * H' * User(l,k).beam_tx * sqrt(ant_gain) );
        %wgt_mse(l,s,k,n) = 1/R_avg(l,s,k)/mse(l,s,k,n);
        
        lambda = 1/(1 + (1/(User(l,k).mse) - 1)/snr_gap );
        User(l,k).wgt_mse = lambda*User(l,k).wgt_wsr/(User(l,k).mse)^2/snr_gap;
        
    end
end

