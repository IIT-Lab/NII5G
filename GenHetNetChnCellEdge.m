function [Chn path_loss_dB UE_Loc] = GenHetNetChnCellEdge(Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant,L_Macro,L_Pico,ISD,K,power_gap)

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
Num_Cell = 1;
Num_BS = L_Macro+L_Pico;
Num_User = K;
Chn = zeros(Num_Cell,Num_BS,Num_TxAnt_Macro,Num_Cell,Num_User,Num_Rx_Ant);
UE_Loc = zeros(Num_Cell,K);
%%
path_loss_dB = zeros(Num_Cell,Num_BS,Num_Cell,Num_User);

shadow_SD_Macro = 8;
shadow_SD_Pico = 10;

ant_gain_Macro = 8;%in dBi
ant_gain_Pico = 5; %in dBi

%BS location
Macro_Loc = 0;
Pico_Loc = zeros(1,L_Pico);
temp_x = ISD/sqrt(3)/2/2;
temp_y = ISD/sqrt(3)/2/2*sqrt(3);
Pico_Loc(1) = temp_x + 1i*temp_y;
Pico_Loc(2) = -temp_x + 1i*temp_y;
Pico_Loc(3) = -temp_x - 1i*temp_y;
Pico_Loc(4) = temp_x - 1i*temp_y;

ICO = 3;

%for each mobile
for m = 1:Num_Cell
    for k = 1:Num_User
        k
        while 1>0
            %generate user location until meets ICO condition
            temp_x = (rand*2 - 1)*ISD/sqrt(3);
            temp_y = (rand*2 - 1)*ISD/2;
            if (temp_y < -sqrt(3)*temp_x + ISD) ...    %right boundary
                    && (temp_y > sqrt(3)*temp_x - ISD)...  %down_right boundary
                    && (temp_y > -sqrt(3)*temp_x - ISD) ...  %down_left boundary
                    && (temp_y < sqrt(3)*temp_x + ISD)  ...     %left boundary
                    && sqrt(temp_x^2 + temp_y^2) >= 0.035 %user cannot be too close to Macro BS
                
                temp_UE = temp_x + 1i*temp_y;
                UE2Macro = abs(temp_UE);
                UE2Pico = abs(temp_UE - Pico_Loc);
                path_loss_dB_Macro = 0;
                path_loss_dB_Pico = zeros(Num_Cell,L_Pico);
                for l = 1:Num_Cell
                    for n = 1:Num_BS
                        
                        if n == 1 %Macro BS
                            ant_gain = ant_gain_Macro;
                            shadowing_dB = shadow_SD_Macro*randn;
                            path_loss_dB_Macro = 128.1 + 37.6*log10(UE2Macro) + shadowing_dB - ant_gain;
                            
                        else
                            ant_gain = ant_gain_Pico;
                            shadowing_dB = shadow_SD_Pico*randn;
                            path_loss_dB_Pico(n-1) = 140.7 + 36.7*log10(UE2Pico(n-1)) + shadowing_dB - ant_gain;
                            
                        end
                        
                    end
                end
                RcvPower = [-path_loss_dB_Macro -power_gap-path_loss_dB_Pico];
                SortRcvPower = sort(RcvPower,'descend');
                if SortRcvPower(1) - SortRcvPower(3) <= ICO
                    UE_Loc(m,k) = temp_UE;
                    
                    for l=1:Num_Cell
                        for n = 1:Num_BS
                            
                            if n == 1
                                path_loss_dB(l,n,m,k) = path_loss_dB_Macro;
                                for iTx = 1:Num_TxAnt_Macro
                                    for iRx = 1:Num_Rx_Ant
                                        
                                        fd = (randn + 1i*randn)*sqrt(1/2); % 0dB
                                        Chn(l,n,iTx,m,k,iRx) = 10^(-path_loss_dB(l,n,m,k)/20) * fd;
                                        
                                    end
                                end
                                
                            else
                                path_loss_dB(l,n,m,k) = path_loss_dB_Pico(n-1);
                                for iTx = 1:Num_TxAnt_Pico
                                    for iRx = 1:Num_Rx_Ant
                                        
                                        fd = (randn + 1i*randn)*sqrt(1/2); % 0dB
                                        Chn(l,n,iTx,m,k,iRx) = 10^(-path_loss_dB(l,n,m,k)/20) * fd;
                                        
                                    end
                                end
                            end
                            
                        end
                    end
                    break;
                end
                
            end
        end
        
    end
end
figure;
hold on;
plot(real(Macro_Loc),imag(Macro_Loc),'r*');
plot(real(Pico_Loc),imag(Pico_Loc),'m*');
plot(real(UE_Loc(1,:)),imag(UE_Loc(1,:)),'bo');
grid on;
saveas(gcf,'Topology7_1','fig');
save Channel7_1.mat Chn path_loss_dB UE_Loc Num_TxAnt_Macro Num_TxAnt_Pico Num_Rx_Ant L_Macro L_Pico ISD K power_gap ICO;