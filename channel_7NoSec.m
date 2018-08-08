function [Chn, path_loss_dB, dist] = channel_7NoSec(L, S, K, P, Q, N, B2Bdist)

%%%%------------------------------------------
%%% This function generates channel matrix in frequency domain     
%%% for frequency-selective fading channels based on a 7-cell wrap around topology
%%% No sectorization
%%%%---------------------------------------

%%% Input format
% L: base-stations
% S: sectors per base-station
% K: users per sector
% N: frequency tones
% P: base-station antennas
% Q: mobile antennas
% B2Bdist: Base-station to base-station distance

%%% Output format (complex gain)
% Chn(l, s, m, t, k, p, q, n)
% complex gain from the l-th BS, s-th sector, p-th antenna
% to the kth user in the m-th BS, t-th sector, q-th antenna
% in n-th frequency tone

%%%------------------------------------
%%% Initialize parameters
%%%-----------------------------------
Chn = zeros(L,S,L,S,K,P,Q,N);  % Channel matrix from the (l,s)-th BS to 
                               % (m,t,k)-th mobile. (p,q) antenna, tone n.
                               
%%%-----------------------------------------------------------------
%%% Path loss: depending on distance between the BS and each of the mobiles
%%% Log-normal Shadowing: 8dB standard deviation. Correlation 0.5 between
%%%     BSs (same shadowing among sections for the same BS.
%%%-----------------------------------------------------------------

%%% Base station fixed locations for the 19 cells
BSLoc = 0.5*B2Bdist* [0; sqrt(3)+j; 2*j; -sqrt(3)+j; -sqrt(3)-j; -2*j; sqrt(3)-j; 
	sqrt(12); sqrt(12)+2*j; sqrt(3)+3*j; 4*j; -sqrt(3)+3*j; -sqrt(12)+2*j; -sqrt(12); 
	-sqrt(12)-2*j; -sqrt(3)-3*j; -4*j; sqrt(3)-3*j; sqrt(12)-2*j];

dist = zeros(L,L,K);
%%
%%%---generate users

MULoc = zeros(L,S,K);
if S == 1 %no sectorization
    for l=1:L,
     k = 1;
        while k <= K
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
                    && sqrt(temp_x^2 + temp_y^2) >= 0.035

              MULoc(l,S,k) = temp_x + 1i*temp_y + BSLoc(l);   
              %if (x,y) is in the hexagon, set it as a user location
              
              k = k+1;%user number increase by 1

            end
        end        
    end;  
end;

%%
%%%---generate CSI from each BS to each user

path_loss_dB = zeros(L,S,L,S,K);
path_loss = zeros(L,L,K);
cov = 0.5*ones(L,L)+0.5*eye(L); % creating shadowing correlation of 0.5
cov_factor = chol(cov);

%%%-----------------------------------
%%% Multipath Time Delay Profile
%%%-----------------------------------
% h_dB = [0 -9.7 -19.2  -Inf -22.8];      %%% ITU-R M.1225 Pedestrain A model 
% Tap = length(h_dB);
% h = 10.^(h_dB/20);

Tap = 1;

%%%----
%%% Shadowing
%%%---
shadow_SD = 8;

%%%--------------
%%% Antenna pattern parameters
%%%--------------
A_min =  0; %dB
%theta_3dB = 2*pi/360*70;


%%-- wrap around pattern
neighbor = [1 2 3 4 5 6 7;
    2  9 10  3  1  7  8;
    3 10 11 12  4  1  2;
    4  3 12 13 14  5  1;
    5  1  4 14 15 16  6;
    6  7  1  5 16 17 18;
    7  8  2  1  6 18 19];
    
mapping = [1; 2; 3; 4; 5; 6; 7; 4; 6; 5; 7; 6; 2; 7; 3; 2; 4; 3; 5];


% for each mobile
for m=1:L,
  for t=1:S,
    for k=1:K,
        
        % shadowing from mobile to each of the BS (correlated)    
        shadowing_dB = shadow_SD*randn(L,1); % independent log-norm shadowing from L BS with 0dB stadard deviation
        shadowing_dB = cov_factor'*shadowing_dB; % 0.5 correlation between BS
        
        % path loss from the mobile to each of the BS/sector
        for l=1:L,
            
          vl = neighbor(m,l);   % equivalent neighbour in the 19-cell configuration
          vb = mapping(vl);     % index of the wrapped around neighbour
          
          dist(vb,m,k) = abs(BSLoc(vl) - MULoc(m,t,k));
  %        theta = angle(MULoc(m,t,k) - BSLoc(vl));
                      
          for s=1:S,
             
  %          ref_theta = 2*pi/12 + 2*pi/S*(s-1);
   
            % [l, s, dist, theta/(2*pi)*360, ref_theta/(2*pi)*360]
            % angle(exp(j*(theta - ref_theta)))/(2*pi)*360
            antenna_pattern = A_min;
            
          	path_loss_dB(vb,s,m,t,k) = 128.1 + 37.6*log10(dist(vb,m,k)) + shadowing_dB(vb) + antenna_pattern;
        
            path_loss(vb,m,k) = path_loss_dB(vb,s,m,t,k);
            for p=1:P,
                for q=1:Q,
                    
                    %%%-------------------------------------
                    %%% Small Scale Fading: Rayleigh Fading 
                    %%% independent fading for each pair of anteanns
                    %%%-------------------------------------
                    fd = (randn(1,Tap) + j*randn(1,Tap))*sqrt(1/2);     %%% 0dB variance
                
                    Chn(vb,s,m,t,k,p,q,:) = 10^(-path_loss_dB(vb,s,m,t,k)/20) * fd;
         
                end;
            end;   
          end;
        end; 
    end;
  end;
end;

figure;
hold on;
plot(real(BSLoc(1:L)),imag(BSLoc(1:L)),'r*');
for l=1:L,
    for s=1:S,
        loc = zeros(K,1);
        for k=1:K;
            loc(k) = MULoc(l,s,k);
        end;
        plot(real(loc),imag(loc),'bo'); 
    end;
end;
grid on; 
legend('Base Station', 'Mobile User')
saveas(gcf,'topology.fig','fig') 

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

 save('Chn_forPdist.mat','Chn', 'BSLoc','MULoc','dist','path_loss');
 return;
