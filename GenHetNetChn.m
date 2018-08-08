function [ChnNew] = GenHetNetChn(Dist,Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function generates channel for HetNet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input
% Dist(l,n,m,k): Distance from l'th cell, n'th BS to k'th user in m'th
% cell

%% Output
%Chn(l,n,a,m,k,b): Channel from a Tx antenna of n BS in l cell to b Rx
%antenna of k user in m cell

%% Initialization
%Num_Cell = size(Dist,1);
Num_Cell = 3;
Num_BS = size(Dist,2);
Num_User = size(Dist,4);
ChnNew = zeros(Num_Cell,Num_BS,Num_TxAnt_Macro,Num_Cell,Num_User,Num_Rx_Ant);

%%
path_loss_dB = zeros(Num_Cell,Num_BS,Num_Cell,Num_User);

shadow_SD_Macro = 8;
shadow_SD_Pico = 10;

ant_gain_Macro = 8;%in dBi
ant_gain_Pico = 5; %in dBi

%for each mobile
for m = 1:Num_Cell
    for k = 1:Num_User
        
        for l = 1:Num_Cell
            for n = 1:Num_BS
                
                if n == 1 %Macro BS
                    ant_gain = ant_gain_Macro;
                    shadowing_dB = shadow_SD_Macro*randn;
                    %fprintf('l=%d, n=%d, m=%d, k=%d \n',l,n,m,k)
                    path_loss_dB(l,n,m,k) = 128.1 + 37.6*log10(Dist(l,n,m,k)) + shadowing_dB - ant_gain;
                    
                    for iTx = 1:Num_TxAnt_Macro
                        for iRx = 1:Num_Rx_Ant
                            %%%-------------------------------------
                            %%% Small Scale Fading: Rayleigh Fading
                            %%% independent fading for each pair of anteanns
                            %%%-------------------------------------
                            
                            fd = (randn + 1i*randn)*sqrt(1/2); % 0dB
                            ChnNew(l,n,iTx,m,k,iRx) = 10^(-path_loss_dB(l,n,m,k)/20) * fd;
                        end
                    end

                else
                    ant_gain = ant_gain_Pico;
                    shadowing_dB = shadow_SD_Pico*randn;
                    path_loss_dB(l,n,m,k) = 140.7 + 36.7*log10(Dist(l,n,m,k)) + shadowing_dB - ant_gain;
                    for iTx = 1:Num_TxAnt_Pico
                        for iRx = 1:Num_Rx_Ant
                            %%%-------------------------------------
                            %%% Small Scale Fading: Rayleigh Fading
                            %%% independent fading for each pair of anteanns
                            %%%-------------------------------------
                            
                            fd = (randn + 1i*randn)*sqrt(1/2); % 0dB
                            ChnNew(l,n,iTx,m,k,iRx) = 10^(-path_loss_dB(l,n,m,k)/20) * fd;
                        end
                    end
                end

            end
        end
        
    end
end
save Chn5.mat ChnNew path_loss_dB;