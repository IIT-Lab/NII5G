function [Dist Macro_Loc Pico_Loc UE_Loc]= GenHetNetTopology(L_Macro,L_Pico,ISD,K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  This function generates the topology for a HetNet single cell %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input
%L_Macro: Number of Macro BS's (also the number of cells)
%L_Pico: Number of Pico BS's
%ISD: inter-sector distance
%K: Number of users per cell

%% Output
%Dist(l,n,m,k): the distance from the n'th BS in l'th cell to the k'th user
%in m'th cell

%% Initialization
Dist = zeros(L_Macro,L_Macro+L_Pico,L_Macro,K);
K_hotspot = 10; %number of users around a Pico BS
radius_hotspot = 40/1000; % K_hotspot users will be placed around radius_hotspot in km around Pico BS's
K_remain = K - K_hotspot*L_Pico; %K_remain users will be placed uniformly in the cell
%%
if L_Macro == 1
    Macro_Loc = 0;
    Pico_Loc = zeros(1,L_Pico);
    Pico_index = 1;
    %generate Pico BS locations
    while Pico_index <= L_Pico
        temp_x = (rand*2 - 1)*ISD/sqrt(3);
        temp_y = (rand*2 - 1)*ISD/2;
        if (temp_y < -sqrt(3)*temp_x + ISD) ...    %right boundary
                && (temp_y > sqrt(3)*temp_x - ISD)...  %down_right boundary
                && (temp_y > -sqrt(3)*temp_x - ISD) ...  %down_left boundary
                && (temp_y < sqrt(3)*temp_x + ISD)  ...     %left boundary
                && sqrt(temp_x^2 + temp_y^2) >= 75/1000 %Pico BS cannot be too close to Macro:75m
            
            temp_pico = temp_x + 1i*temp_y;
            if Pico_index > 1
                P2PDist = abs((temp_pico - Pico_Loc(1:Pico_index-1)));
                min_P2PDist = min(P2PDist);
                if min_P2PDist >= 40/1000   %Pico to Pico cannot be too close:40m
                    Pico_Loc(Pico_index) = temp_pico;
                    Pico_index = Pico_index + 1;
                end
            else
                Pico_Loc(Pico_index) = temp_pico;
                Pico_index = Pico_index + 1;
            end
        end
    end
    
    UE_Loc = zeros(1,K);
    user_index = 1;
    %generate hotspot users
    for Pico_index = 1:L_Pico
        hotspot_user_index = 1;
        while hotspot_user_index <= K_hotspot
            temp_x = (rand*2 - 1)*radius_hotspot + real(Pico_Loc(Pico_index));
            temp_y = (rand*2 - 1)*radius_hotspot + imag(Pico_Loc(Pico_index));
            if (temp_y < -sqrt(3)*temp_x + ISD) ...    %right boundary
                    && (temp_y > sqrt(3)*temp_x - ISD)...  %down_right boundary
                    && (temp_y > -sqrt(3)*temp_x - ISD) ...  %down_left boundary
                    && (temp_y < sqrt(3)*temp_x + ISD)  ...     %left boundary
                    && sqrt(temp_x^2 + temp_y^2) >= 35/1000     %UE cannot be too close to Macro:35m
                
                temp_UE = temp_x + 1i*temp_y;
                if min(abs(temp_UE - Pico_Loc)) >= 3/1000; %UE cannot be too close to Pico:3m
                    UE_Loc(user_index) = temp_x + 1i*temp_y;
                    hotspot_user_index = hotspot_user_index + 1;
                    user_index = user_index + 1;
                end
            end
        end
    end
    
    %generate the rest of the users
    while user_index <= K
        temp_x = (rand*2 - 1)*ISD/sqrt(3);
        temp_y = (rand*2 - 1)*ISD/2;
        if (temp_y < -sqrt(3)*temp_x + ISD) ...    %right boundary
                && (temp_y > sqrt(3)*temp_x - ISD)...  %down_right boundary
                && (temp_y > -sqrt(3)*temp_x - ISD) ...  %down_left boundary
                && (temp_y < sqrt(3)*temp_x + ISD)  ...     %left boundary
                && sqrt(temp_x^2 + temp_y^2) >= 0.035 %user cannot be too close to Macro BS
            
            temp_UE = temp_x + 1i*temp_y;
            UE2Pico = abs(temp_UE - Pico_Loc);
            min_UE2Pico = min(UE2Pico);
            if min_UE2Pico >= 3/1000
                UE_Loc(user_index) = temp_UE;
                user_index = user_index + 1;
            end
        end
    end
    
    %calculate distance
    for k = 1:K
        Dist(L_Macro,1,L_Macro,k) = abs(UE_Loc(k) - Macro_Loc);%distance between Macro BS to user k
        for Pico_index = 1:L_Pico
            Dist(L_Macro,Pico_index+1,L_Macro,k) = abs(UE_Loc(k) - Pico_Loc(Pico_index));
        end
    end
    
    figure;
    hold on;
    plot(real(Macro_Loc),imag(Macro_Loc),'r*');
    plot(real(Pico_Loc),imag(Pico_Loc),'m*');
    plot(real(UE_Loc(1,:)),imag(UE_Loc(1,:)),'bo');
    grid on;
    saveas(gcf,'Topology','fig');
else
    fprintf('$$$ HetNet Topology has not generalized to multiple cells /n $$$');
end
save topology.mat Dist Macro_Loc Pico_Loc UE_Loc;
