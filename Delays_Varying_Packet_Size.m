PacketSize = [1e-4 1.78e-4 3.16e-4 5.62e-4 1e-3 1.78e-3 3.16e-3 5.62e-3 1e-2 1.78e-2 3.16e-2 5.62e-2 1e-1];
Tavg = zeros(4,length(PacketSize));

for size = 1:length(PacketSize)
    res = Delays(300e-6,PacketSize(size),Simulations,1e-3,20,3,10);
    s_pBF = 0;
    s_BF2 = 0;
    s_pClust = 0;
    s_Clust2 = 0;
    
    for pos = 4:27
        if pos == 10
            continue;
        end
        for i = 1:3
            for j = 1:20
                s_pBF = s_pBF + res(pos,i,j).pCSI_BF;
            end
        end
    end
    Tavg(1,size) = s_pBF / (23*60);
    
    for pos = 4:27
        if pos == 10
            continue;
        end
        for i = 1:3
            for j = 1:20
                s_BF2 = s_BF2 + res(pos,i,j).iCSI_BF_2;
            end
        end
    end
    Tavg(2,size) = s_BF2 / (23*60);
    
    for pos = 4:27
        if pos == 10
            continue;
        end
        for i = 1:3
            for j = 1:20
                s_pClust = s_pClust + res(pos,i,j).pCSI_Clust;
            end
        end
    end
    Tavg(3,size) = s_pClust / (23*60);
    
    for pos = 4:27
        if pos == 10
            continue;
        end
        for i = 1:3
            for j = 1:20
                s_Clust2 = s_Clust2 + res(pos,i,j).iCSI_Clust_2;
            end
        end
    end
    Tavg(4,size) = s_Clust2 / (23*60);
end

