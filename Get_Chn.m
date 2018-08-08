  function H = Get_Chn(Chn, BS_set, l,k, L_Macro,L_Pico,Num_TxAnt_Pico)

%%%% -------------------------------------------------------------------------
%%% This function finds the channel matrix from the BSs stored in vector BS_set to user
%%% (l,k)
%%%% -------------------------------------------------------------------------

H = [];
Total_BS = L_Macro + L_Pico;

for iBS = 1:length(BS_set)
    curr_cell = ceil(BS_set(iBS)/Total_BS);
    curr_BS = BS_set(iBS) - (curr_cell-1)*Total_BS;
    if curr_BS == 1  %Macro BS
        H_temp = squeeze(Chn(curr_cell,curr_BS,:,l,k,:));
    else  %Pico BS
        H_temp = squeeze(Chn(curr_cell,curr_BS,1:Num_TxAnt_Pico,l,k,:));
    end
    H = [H;H_temp]; 
    
end

%H is the size of 
% row: number of antennas serving user (l,k)
% column: number of receive antennas