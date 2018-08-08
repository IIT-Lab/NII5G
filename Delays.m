function [res] = Delays(Delta, PacketLength, Simulations, FrameLength, K, Num_Cell, ClusteringPeriod)

res = struct([]);
for sim = 1:size(Simulations,2)
    load(Simulations(sim).File);
    Num_Ite = size(Log_Data,3);
    for pos = Simulations(sim).Start:Simulations(sim).Stop
        for l = 1:Num_Cell
            for k = 1:K
                
                data = 0;
                for ite = 1:Num_Ite
                    data = data + Log_Data(4,pos,ite).pCSI_User(l,k).inst_rate * FrameLength;
                    if data >= PacketLength
                        res(pos,l,k).pCSI_BF = ite * (FrameLength + Delta/ClusteringPeriod);
                        break;
                    end
                    res(pos,l,k).pCSI_BF = (Num_Ite + 1) * (FrameLength + Delta/ClusteringPeriod);
                end
                
                data = 0;
                for ite = 1:Num_Ite
                    data = data + Log_Data(4,pos,ite).BF(1).User(l,k).inst_rate * FrameLength;
                    if data >= PacketLength
                        res(pos,l,k).iCSI_BF_1 = ite * (FrameLength + Delta/ClusteringPeriod);
                        break;
                    end
                    res(pos,l,k).iCSI_BF_1 = (Num_Ite + 1) * (FrameLength + Delta/ClusteringPeriod);
                end
                
                data = 0;
                for ite = 1:Num_Ite
                    data = data + Log_Data(4,pos,ite).BF(2).User(l,k).inst_rate * FrameLength;
                    if data >= PacketLength
                        res(pos,l,k).iCSI_BF_2 = ite * (FrameLength + Delta/ClusteringPeriod);
                        break;
                    end
                    res(pos,l,k).iCSI_BF_2 = (Num_Ite + 1) * (FrameLength + Delta/ClusteringPeriod);
                end
                
                data = 0;
                for ite = 1:Num_Ite
                    data = data + Log_Data(4,pos,ite).BF(3).User(l,k).inst_rate * FrameLength;
                    if data >= PacketLength
                        res(pos,l,k).iCSI_BF_3 = ite * (FrameLength + Delta/ClusteringPeriod);
                        break;
                    end
                    res(pos,l,k).iCSI_BF_3 = (Num_Ite + 1) * (FrameLength + Delta/ClusteringPeriod);
                end
                
                data = 0;
                for ite = 1:Num_Ite
                    data = data + Log_Data(4,pos,ite).Clustering_User(l,k).inst_rate * FrameLength;
                    if data >= PacketLength
                        res(pos,l,k).pCSI_Clust = ite * (FrameLength + Delta);
                        break;
                    end
                    res(pos,l,k).pCSI_Clust = (Num_Ite + 1) * (FrameLength + Delta);
                end
                
                data = 0;
                for ite = 1:Num_Ite
                    data = data + Log_Data(4,pos,ite).Clustering_sumrate_iCSI(1).User(l,k).inst_rate * FrameLength;
                    if data >= PacketLength
                        res(pos,l,k).iCSI_Clust_1 = ite * (FrameLength + Delta);
                        break;
                    end
                    res(pos,l,k).iCSI_Clust_1 = (Num_Ite + 1) * (FrameLength + Delta);
                end
                
                data = 0;
                for ite = 1:Num_Ite                   
                    data = data + Log_Data(4,pos,ite).Clustering_sumrate_iCSI(2).User(l,k).inst_rate * FrameLength;
                    if data >= PacketLength
                        res(pos,l,k).iCSI_Clust_2 = ite * (FrameLength + Delta);
                        break;
                    end
                    res(pos,l,k).iCSI_Clust_2 = (Num_Ite + 1) * (FrameLength + Delta);
                end
                
                data = 0;
                for ite = 1:Num_Ite                   
                    data = data + Log_Data(4,pos,ite).Clustering_sumrate_iCSI(3).User(l,k).inst_rate * FrameLength;
                    if data >= PacketLength
                        res(pos,l,k).iCSI_Clust_3 = ite * (FrameLength + Delta);
                        break;
                    end
                    res(pos,l,k).iCSI_Clust_3 = (Num_Ite + 1) * (FrameLength + Delta);
                end
            end
        end
    end
end