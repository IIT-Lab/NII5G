function [ChnNew] = GenNewCSI(Num_Cell, S, Num_User, Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,path_loss_dB)

%%%---generate CSI from each BS to each user
Num_BS = 4; %1 Macro BS, 3 Pico BS's per Cell
ChnNew = zeros(Num_Cell,Num_BS,Num_TxAnt_Macro,Num_Cell,Num_User,Num_Rx_Ant);

neighbor = [1 2 3 4 5 6 7;
    2  9 10  3  1  7  8;
    3 10 11 12  4  1  2;
    4  3 12 13 14  5  1;
    5  1  4 14 15 16  6;
    6  7  1  5 16 17 18;
    7  8  2  1  6 18 19];

mapping = [1; 2; 3; 4; 5; 6; 7; 4; 6; 5; 7; 6; 2; 7; 3; 2; 4; 3; 5];

%%
% for each mobile
for m=1:Num_Cell
    for t=1:S
        for k=1:Num_User
            
            % shadowing from mobile to each of the Macro BS (correlated)
            %shadowing_dB_Macro = shadow_SD_Macro*randn(Num_Cell,1); % independent log-norm shadowing from Macro BS with 8dB stadard deviation
            %shadowing_dB_Macro = cov_factor'*shadowing_dB_Macro; % 0.5 correlation between BS
            
            % path loss from the mobile to each of the BS/sector
            for l=1:Num_Cell
                
                vl = neighbor(m,l);   % equivalent neighbour in the 19-cell configuration
                vb = mapping(vl);     % index of the wrapped around neighbour
                
                for n = 1:Num_BS
                    
                    if n == 1 %Macro BS
                        
                        for p=1:Num_TxAnt_Macro
                            for q=1:Num_Rx_Ant
                                
                                %%%-------------------------------------
                                %%% Small Scale Fading: Rayleigh Fading
                                %%% independent fading for each pair of anteanns
                                %%%-------------------------------------
                                fd = (randn + j*randn)*sqrt(1/2);     %%% 0dB variance
                                
                                ChnNew(vb,n,p,m,k,q) = 10^(-path_loss_dB(vb,n,m,k)/20) * fd;
                                
                            end;
                        end;
                    else  %Pico BS
                        
                        for p=1:Num_TxAnt_Pico
                            for q=1:Num_Rx_Ant
                                
                                %%%-------------------------------------
                                %%% Small Scale Fading: Rayleigh Fading
                                %%% independent fading for each pair of anteanns
                                %%%-------------------------------------
                                fd = (randn + j*randn)*sqrt(1/2);     %%% 0dB variance
                                
                                ChnNew(vb,n,p,m,k,q) = 10^(-path_loss_dB(vb,n,m,k)/20) * fd;
                                
                            end;
                        end;
                    end
                end;
            end
        end;
    end;
end;

%save('ChnNew.mat','ChnNew', 'BSLoc','MULoc','dist','path_loss_dB','Num_BS','Num_User');
return;
