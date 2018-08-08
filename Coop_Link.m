function E = Coop_Link(L,L_Macro,L_Pico,K,Num_TxAnt_Macro,Num_TxAnt_Pico,User)

%%% ---------------
%%% This function finds the cooperation link matrix E for each user
%%% according to the Clustering Scheme
%%% ---------------



E = zeros(P*S*L,P*S*L,L,K);

for l = 1:L
    for s = 1:S
        for k = 1:K
            for n = 1:N
                nbr_cell = squeeze(nbr_order(l,s,k,:));
                temp = zeros(P*S*L,1);
                
                for m = 1:M
                    index = nbr_cell(m);
                    temp( (index-1)*P*S+1 : index*P*S  ) = ones(P*S,1);
                end
                
                E(:,:,l,s,k,n) = diag(temp);
                
            end
        end
    end
end