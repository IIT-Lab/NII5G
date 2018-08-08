function [iChn] = Gen_iCSI(Chn,sigma_iCSI)

iChn = Chn + (randn(size(Chn)) + 1i*randn(size(Chn))) * (sigma_iCSI/sqrt(2));

return;

