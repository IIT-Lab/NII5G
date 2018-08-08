% convert from dB to decimal

function dB = dec2dB(dec);

dB = 10*log10(dec);