% convert from dB to decimal

function dec = dB2dec(dB);

dec = 10.^(0.1*dB);