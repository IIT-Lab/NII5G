function [Chn, path_loss_dB, dist, MULoc] = GenHetNetChn7Cell(Num_Cell, S, Num_User, Num_TxAnt_Macro,Num_TxAnt_Pico,Num_Rx_Ant, B2Bdist)

%%%%------------------------------------------
%%% This function generates channel matrix in frequency domain
%%% for frequency-selective fading channels based on a 7-cell wrap around topology
%%% No sectorization
%%% Two Tier Het Net
%%%%---------------------------------------

Num_BS = 4; %1 Macro BS, 3 Pico BS's per Cell
Chn = zeros(Num_Cell,Num_BS,Num_TxAnt_Macro,Num_Cell,Num_User,Num_Rx_Ant);

%%%-----------------------------------------------------------------
%%% Path loss: depending on distance between the BS and each of the mobiles
%%% Log-normal Shadowing: 8dB standard deviation. Correlation 0.5 between
%%%     BSs (same shadowing among sections for the same BS.
%%%-----------------------------------------------------------------

%%% Base station fixed locations for the 19 cells
BSLoc = 0.5*B2Bdist* [0; sqrt(3)+j; 2*j; -sqrt(3)+j; -sqrt(3)-j; -2*j; sqrt(3)-j;
    sqrt(12); sqrt(12)+2*j; sqrt(3)+3*j; 4*j; -sqrt(3)+3*j; -sqrt(12)+2*j; -sqrt(12);
    -sqrt(12)-2*j; -sqrt(3)-3*j; -4*j; sqrt(3)-3*j; sqrt(12)-2*j];

PicoLoc = 0.5*B2Bdist*[1/sqrt(3); (cos(2*pi/3)+j*sin(2*pi/3))/sqrt(3); (cos(2*pi/3)+j*sin(-2*pi/3))/sqrt(3)];

dist = zeros(Num_Cell,Num_BS,Num_Cell,Num_User);
%%
%%%---generate users

MULoc = zeros(Num_Cell,S,Num_User);
if S == 1 %no sectorization
    for l=1:Num_Cell,
        k = 1;
        while k <= Num_User
            temp_x = (rand*2 - 1)*B2Bdist/sqrt(3);
            %shift uniform distribution among (0,1) to (-B2Bdist/sqrt(3),B2Bdist/sqrt(3))
            %where x can possibly locate.
            temp_y = (rand*2 - 1)*B2Bdist/2;
            %shift uniform distribution among (0,1) to
            %(-B2Bdist/2,B2Bdist/2) where y can possibly locate.
            
            if (temp_y < -sqrt(3)*temp_x + B2Bdist) ...    %right boundary
                    && (temp_y > sqrt(3)*temp_x - B2Bdist)...  %down_right boundary
                    && (temp_y > -sqrt(3)*temp_x - B2Bdist) ...  %down_left boundary
                    && (temp_y < sqrt(3)*temp_x + B2Bdist)  ...     %left boundary
                    && sqrt(temp_x^2 + temp_y^2) >= 0.035       %min dist from Macro to UE: 35m
                
                
                minUE2PicoDist = min(abs(temp_x + j*temp_y - PicoLoc));
                
                if minUE2PicoDist >= 3/1000         %min dist from Pico to UE: 3m
                    MULoc(l,S,k) = temp_x + 1i*temp_y + BSLoc(l);
                    %if (x,y) is in the hexagon, set it as a user location
                    
                    k = k+1;%user number increase by 1
                end
            end
        end
    end;
end;

%%
%%%---generate CSI from each BS to each user

path_loss_dB = zeros(Num_Cell,Num_BS,Num_Cell,Num_User);
% path_loss = zeros(L,L,K);
cov = 0.5*ones(Num_Cell,Num_Cell)+0.5*eye(Num_Cell); % creating shadowing correlation of 0.5
cov_factor = chol(cov);

%%%-----------------------------------
%%% Multipath Time Delay Profile
%%%-----------------------------------
% h_dB = [0 -9.7 -19.2  -Inf -22.8];      %%% ITU-R M.1225 Pedestrain A model
% Tap = length(h_dB);
% h = 10.^(h_dB/20);

%Tap = 1;

%%%----
%%% Shadowing
%%%---
%shadow_SD = 8;

%%%--------------
%%% Antenna pattern parameters
%%%--------------
%A_min =  0; %dB
%theta_3dB = 2*pi/360*70;

shadow_SD_Macro = 8;
shadow_SD_Pico = 8;

ant_gain_Macro = 15;%in dBi
ant_gain_Pico = 15; %in dBi

%%-- wrap around pattern
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
for m=1:Num_Cell,
    for t=1:S,
        for k=1:Num_User,
            
            % shadowing from mobile to each of the Macro BS (correlated)
            shadowing_dB_Macro = shadow_SD_Macro*randn(Num_Cell,1); % independent log-norm shadowing from Macro BS with 8dB stadard deviation
            shadowing_dB_Macro = cov_factor'*shadowing_dB_Macro; % 0.5 correlation between BS
            
            % path loss from the mobile to each of the BS/sector
            for l=1:Num_Cell,
                
                vl = neighbor(m,l);   % equivalent neighbour in the 19-cell configuration
                vb = mapping(vl);     % index of the wrapped around neighbour
                
                for n = 1:Num_BS
                    
                    if n == 1 %Macro BS
                        
                        dist(vb,n,m,k) = abs(BSLoc(vl) - MULoc(m,t,k));
                        ant_gain = ant_gain_Macro;
                        
                        path_loss_dB(vb,n,m,k) = 128.1 + 37.6*log10(dist(vb,n,m,k)) + shadowing_dB_Macro(vb) - ant_gain;
                        
                        for p=1:Num_TxAnt_Macro,
                            for q=1:Num_Rx_Ant,
                                
                                %%%-------------------------------------
                                %%% Small Scale Fading: Rayleigh Fading
                                %%% independent fading for each pair of anteanns
                                %%%-------------------------------------
                                fd = (randn + j*randn)*sqrt(1/2);     %%% 0dB variance
                                
                                Chn(vb,n,p,m,k,q) = 10^(-path_loss_dB(vb,n,m,k)/20) * fd;
                                
                            end;
                        end;
                    else  %Pico BS
                        dist(vb,n,m,k) = abs(BSLoc(vl) + PicoLoc(n-1) - MULoc(m,t,k));
                        ant_gain = ant_gain_Pico;
                        shadowing_dB_Pico = shadow_SD_Pico*randn;
                        path_loss_dB(vb,n,m,k) = 140.7 + 36.7*log10(dist(vb,n,m,k)) + shadowing_dB_Pico - ant_gain;
                        
                        for p=1:Num_TxAnt_Pico,
                            for q=1:Num_Rx_Ant,
                                
                                %%%-------------------------------------
                                %%% Small Scale Fading: Rayleigh Fading
                                %%% independent fading for each pair of anteanns
                                %%%-------------------------------------
                                fd = (randn + j*randn)*sqrt(1/2);     %%% 0dB variance
                                
                                Chn(vb,n,p,m,k,q) = 10^(-path_loss_dB(vb,n,m,k)/20) * fd;
                                
                            end;
                        end;
                    end
                end;
            end
        end;
    end;
end;

AllBSLoc = zeros(Num_Cell,Num_BS);

for l = 1:Num_Cell
    for n = 1:Num_BS
        if n == 1
            AllBSLoc(l,n) = BSLoc(l);
        else
            AllBSLoc(l,n) = BSLoc(l) + PicoLoc(n-1);
        end
        
    end
end


%figure;
%hold on;
%for l = 1:Num_Cell
%    for n = 1:Num_BS
%        
%        if n == 1
%            plot(real(AllBSLoc(l,n)),imag(AllBSLoc(l,n)),'rs');           
%        else
%            plot(real(AllBSLoc(l,n)),imag(AllBSLoc(l,n)),'r*');
%        end
%
%    end
%end



%for l=1:Num_Cell,
%for l=1:3
%    for s=1:S,
%        loc = zeros(Num_User,1);
%        for k=1:Num_User;
%            loc(k) = MULoc(l,s,k);
%        end;
%        plot(real(loc),imag(loc),'bo');
%    end;
%end;
%grid on;
%legend('Macro Base Station', 'Pico Base Station','Mobile User')
%saveas(gcf,'topology.fig','fig')

% figure;
% ch_gain = zeros(N,1);
% for k=1:K,
%   subplot(K,1,k)
%   for n=1:N,
%     ch_gain(n) = dec2dB(norm(Chn(1,1,1,1,k,1,1,n))^2);
%   end;
%   plot(ch_gain);
% end;
% xlabel('frequency tones');
% saveas(gcf,'channel.fig','fig')

%save('Chn5.mat','Chn', 'BSLoc','MULoc','dist','path_loss_dB','Num_BS','Num_User');
return;
