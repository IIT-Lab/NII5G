y = [];

for k = 1:30
    for iCell = 1:3
        y = [y User(iCell,k).long_avg_rate];
    end
end

res2 = [];
for temp = 0:0.05:20
    f = 0.00;
    for x = 1:90
        f = f + (double(temp >= y(x)))/90;
    end
    res2 = [res2 f];
end

%% BF / iCSI

load('Log_Data_proposedTest1_3.mat')
BF_iCSI_1_cdf = [];
for pos = 1:3
    for i = 1:3
        for j = 1:20
            BF_iCSI_1_cdf = [BF_iCSI_1_cdf Log_Data(4,pos,30).BF(1).User(i,j).long_avg_rate];
        end
    end
end

load('Log_Data_proposedTest4_9.mat');
for pos = 4:9
    for i = 1:3
        for j = 1:20
            BF_iCSI_1_cdf = [BF_iCSI_1_cdf Log_Data(4,pos,50).BF(1).User(i,j).long_avg_rate];
        end
    end
end

load('Log_Data_proposedTest11_16.mat');
for pos = 11:16
    for i = 1:3
        for j = 1:20
            BF_iCSI_1_cdf = [BF_iCSI_1_cdf Log_Data(4,pos,50).BF(1).User(i,j).long_avg_rate];
        end
    end
end

BF_iCSI_1_cdf = sort(BF_iCSI_1_cdf);

res_BF_iCSI_1_cdf = [];
for temp = 0:0.05:10
    f = 0.0;
    for x = 1:900
        if temp < BF_iCSI_1_cdf(x)
            break;
        else
            f = f + 1;
        end
    end
    res_BF_iCSI_1_cdf = [res_BF_iCSI_1_cdf f/900];
end

%% BF / pCSI

load('Log_Data_proposedTest1_3.mat')
BF_pCSI_cdf = [];
for pos = 1:3
    for i = 1:3
        for j = 1:20
            BF_pCSI_cdf = [BF_pCSI_cdf Log_Data(4,pos,30).pCSI_User(i,j).long_avg_rate];
        end
    end
end

load('Log_Data_proposedTest4_9.mat');
for pos = 4:9
    for i = 1:3
        for j = 1:20
            BF_pCSI_cdf = [BF_pCSI_cdf Log_Data(4,pos,50).pCSI_User(i,j).long_avg_rate];
        end
    end
end

load('Log_Data_proposedTest11_16.mat');
for pos = 11:16
    for i = 1:3
        for j = 1:20
            BF_pCSI_cdf = [BF_pCSI_cdf Log_Data(4,pos,50).pCSI_User(i,j).long_avg_rate];
        end
    end
end

BF_pCSI_cdf = sort(BF_pCSI_cdf);

res_BF_pCSI_cdf = [];
for temp = 0:0.05:10
    f = 0.0;
    for x = 1:900
        if temp < BF_pCSI_cdf(x)
            break;
        else
            f = f + 1;
        end
    end
    res_BF_pCSI_cdf = [res_BF_pCSI_cdf f/900];
end

%% Clustering / iCSI

load('Log_Data_proposedTest1_3.mat')
pClust_iCSI_1_cdf = [];
for pos = 1:3
    for i = 1:3
        for j = 1:20
            pClust_iCSI_1_cdf = [pClust_iCSI_1_cdf Log_Data(4,pos,30).Clustering_sumrate_iCSI(1).User(i,j).long_avg_rate];
        end
    end
end

load('Log_Data_proposedTest4_9.mat');
for pos = 4:9
    for i = 1:3
        for j = 1:20
            pClust_iCSI_1_cdf = [pClust_iCSI_1_cdf Log_Data(4,pos,50).Clustering_sumrate_iCSI(1).User(i,j).long_avg_rate];
        end
    end
end

load('Log_Data_proposedTest11_16.mat');
for pos = 11:16
    for i = 1:3
        for j = 1:20
            pClust_iCSI_1_cdf = [pClust_iCSI_1_cdf Log_Data(4,pos,50).Clustering_sumrate_iCSI(1).User(i,j).long_avg_rate];
        end
    end
end

pClust_iCSI_1_cdf = sort(pClust_iCSI_1_cdf);

res_Clust_iCSI_1_cdf = [];
for temp = 0:0.05:10
    f = 0.0;
    for x = 1:900
        if temp < pClust_iCSI_1_cdf(x)
            break;
        else
            f = f + 1;
        end
    end
    res_Clust_iCSI_1_cdf = [res_Clust_iCSI_1_cdf f/900];
end

%% Clustering / pCSI

load('Log_Data_proposedTest1_3.mat')
pClust_cdf = [];
for pos = 1:3
    for i = 1:3
        for j = 1:20
            pClust_cdf = [pClust_cdf Log_Data(4,pos,30).Clustering_User(i,j).long_avg_rate];
        end
    end
end

load('Log_Data_proposedTest4.mat')
for i = 1:3
    for j = 1:20
        pClust_cdf = [pClust_cdf Log_Data(4,4,30).Clustering_User(i,j).long_avg_rate];
    end
end
pClust_cdf = sort(pClust_cdf);

res_pClust_cdf = [];
for temp = 0:0.05:10
    f = 0.0;
    for x = 1:240
        if temp < pClust_cdf(x)
            break;
        else
            f = f + 1;
        end
    end
    res_pClust_cdf = [res_pClust_cdf f/240];
end

%%

figure; plot(0:0.05:10,res_pClust_cdf,'k-',0:0.05:10,res_Clust_iCSI_3_cdf,'k--',0:0.05:10,res_Clust_iCSI_2_cdf,'k-.',0:0.05:10,res_Clust_iCSI_1_cdf,'k:',...
    0:0.05:10,res_BF_pCSI_cdf,'r-',0:0.05:10,res_BF_iCSI_3_cdf,'r--',0:0.05:10,res_BF_iCSI_2_cdf,'r-.',0:0.05:10,res_BF_iCSI_1_cdf,'r:');

%% Average Delays

s_pBF = 0;
s_BF1 = 0;
s_BF2 = 0;
s_BF3 = 0;
s_pClust = 0;
s_Clust1 = 0;
s_Clust2 = 0;
s_Clust3 = 0;
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
s_pBF = s_pBF / (23*60);

for pos = 4:27
    if pos == 10
        continue;
    end
    for i = 1:3
        for j = 1:20
            s_BF1 = s_BF1 + res(pos,i,j).iCSI_BF_1;
        end
    end
end
s_BF1 = s_BF1 / (23*60);

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
s_BF2 = s_BF2 / (23*60);

for pos = 4:27
    if pos == 10
        continue;
    end
    for i = 1:3
        for j = 1:20
            s_BF3 = s_BF3 + res(pos,i,j).iCSI_BF_3;
        end
    end
end
s_BF3 = s_BF3 / (23*60);

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
s_pClust = s_pClust / (23*60);

for pos = 4:27
    if pos == 10
        continue;
    end
    for i = 1:3
        for j = 1:20
            s_Clust1 = s_Clust1 + res(pos,i,j).iCSI_Clust_1;
        end
    end
end
s_Clust1 = s_Clust1 / (23*60);

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
s_Clust2 = s_Clust2 / (23*60);

for pos = 4:27
    if pos == 10
        continue;
    end
    for i = 1:3
        for j = 1:20
            s_Clust3 = s_Clust3 + res(pos,i,j).iCSI_Clust_3;
        end
    end
end
s_Clust3 = s_Clust3 / (23*60);

%% BF / iCSI Delay

BF_iCSI_1_cdf = [];
for pos = 4:27
    if pos == 10
        continue;
    end
    for i = 1:3
        for j = 1:20
            BF_iCSI_1_cdf = [BF_iCSI_1_cdf res(pos,i,j).iCSI_BF_1];
        end
    end
end

BF_iCSI_1_cdf = sort(BF_iCSI_1_cdf);

res_BF_iCSI_1_cdf = [];
for temp = 0:0.0001:0.07
    f = 0.0;
    for x = 1:1380
        if temp < BF_iCSI_1_cdf(x)
            break;
        else
            f = f + 1;
        end
    end
    res_BF_iCSI_1_cdf = [res_BF_iCSI_1_cdf f/1380];
end

BF_iCSI_2_cdf = [];
for pos = 4:27
    if pos == 10
        continue;
    end
    for i = 1:3
        for j = 1:20
            BF_iCSI_2_cdf = [BF_iCSI_2_cdf res(pos,i,j).iCSI_BF_2];
        end
    end
end

BF_iCSI_2_cdf = sort(BF_iCSI_2_cdf);

res_BF_iCSI_2_cdf = [];
for temp = 0:0.0001:0.07
    f = 0.0;
    for x = 1:1380
        if temp < BF_iCSI_2_cdf(x)
            break;
        else
            f = f + 1;
        end
    end
    res_BF_iCSI_2_cdf = [res_BF_iCSI_2_cdf f/1380];
end

BF_iCSI_3_cdf = [];
for pos = 4:27
    if pos == 10
        continue;
    end
    for i = 1:3
        for j = 1:20
            BF_iCSI_3_cdf = [BF_iCSI_3_cdf res(pos,i,j).iCSI_BF_3];
        end
    end
end

BF_iCSI_3_cdf = sort(BF_iCSI_3_cdf);

res_BF_iCSI_3_cdf = [];
for temp = 0:0.0001:0.07
    f = 0.0;
    for x = 1:1380
        if temp < BF_iCSI_3_cdf(x)
            break;
        else
            f = f + 1;
        end
    end
    res_BF_iCSI_3_cdf = [res_BF_iCSI_3_cdf f/1380];
end

%% BF / pCSI Delay

BF_pCSI_cdf = [];
for pos = 4:27
    if pos == 10
        continue;
    end
    for i = 1:3
        for j = 1:20
            BF_pCSI_cdf = [BF_pCSI_cdf res(pos,i,j).pCSI_BF];
        end
    end
end

BF_pCSI_cdf = sort(BF_pCSI_cdf);

res_BF_pCSI_cdf = [];
for temp = 0:0.0001:0.07
    f = 0.0;
    for x = 1:1380
        if temp < BF_pCSI_cdf(x)
            break;
        else
            f = f + 1;
        end
    end
    res_BF_pCSI_cdf = [res_BF_pCSI_cdf f/1380];
end

%% Clustering / iCSI Delay

Clust_iCSI_1_cdf = [];
for pos = 4:27
    if pos == 10
        continue;
    end
    for i = 1:3
        for j = 1:20
            Clust_iCSI_1_cdf = [Clust_iCSI_1_cdf res(pos,i,j).iCSI_Clust_1];
        end
    end
end

Clust_iCSI_1_cdf = sort(Clust_iCSI_1_cdf);

res_Clust_iCSI_1_cdf = [];
for temp = 0:0.0001:0.07
    f = 0.0;
    for x = 1:1380
        if temp < Clust_iCSI_1_cdf(x)
            break;
        else
            f = f + 1;
        end
    end
    res_Clust_iCSI_1_cdf = [res_Clust_iCSI_1_cdf f/1380];
end

Clust_iCSI_2_cdf = [];
for pos = 4:27
    if pos == 10
        continue;
    end
    for i = 1:3
        for j = 1:20
            Clust_iCSI_2_cdf = [Clust_iCSI_2_cdf res(pos,i,j).iCSI_Clust_2];
        end
    end
end

Clust_iCSI_2_cdf = sort(Clust_iCSI_2_cdf);

res_Clust_iCSI_2_cdf = [];
for temp = 0:0.0001:0.07
    f = 0.0;
    for x = 1:1380
        if temp < Clust_iCSI_2_cdf(x)
            break;
        else
            f = f + 1;
        end
    end
    res_Clust_iCSI_2_cdf = [res_Clust_iCSI_2_cdf f/1380];
end

Clust_iCSI_3_cdf = [];
for pos = 4:27
    if pos == 10
        continue;
    end
    for i = 1:3
        for j = 1:20
            Clust_iCSI_3_cdf = [Clust_iCSI_3_cdf res(pos,i,j).iCSI_Clust_3];
        end
    end
end

Clust_iCSI_3_cdf = sort(Clust_iCSI_3_cdf);

res_Clust_iCSI_3_cdf = [];
for temp = 0:0.0001:0.07
    f = 0.0;
    for x = 1:1380
        if temp < Clust_iCSI_3_cdf(x)
            break;
        else
            f = f + 1;
        end
    end
    res_Clust_iCSI_3_cdf = [res_Clust_iCSI_3_cdf f/1380];
end

%% Clustering / pCSI Delay

Clust_pCSI_cdf = [];
for pos = 4:27
    if pos == 10
        continue;
    end
    for i = 1:3
        for j = 1:20
            Clust_pCSI_cdf = [Clust_pCSI_cdf res(pos,i,j).pCSI_Clust];
        end
    end
end

Clust_pCSI_cdf = sort(Clust_pCSI_cdf);

res_Clust_pCSI_cdf = [];
for temp = 0:0.0001:0.07
    f = 0.0;
    for x = 1:1380
        if temp < Clust_pCSI_cdf(x)
            break;
        else
            f = f + 1;
        end
    end
    res_Clust_pCSI_cdf = [res_Clust_pCSI_cdf f/1380];
end

%%

figure; plot(0:0.0001:0.07,res_Clust_pCSI_cdf,'k-',0:0.0001:0.07,res_Clust_iCSI_3_cdf,'k--',0:0.0001:0.07,res_Clust_iCSI_2_cdf,'k-.',0:0.0001:0.07,res_Clust_iCSI_1_cdf,'k:',...
    0:0.0001:0.07,res_BF_pCSI_cdf,'r-',0:0.0001:0.07,res_BF_iCSI_3_cdf,'r--',0:0.0001:0.07,res_BF_iCSI_2_cdf,'r-.',0:0.0001:0.07,res_BF_iCSI_1_cdf,'r:');