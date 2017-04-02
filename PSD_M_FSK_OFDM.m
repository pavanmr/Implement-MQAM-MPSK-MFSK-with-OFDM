close all
clear
clc

M = 16; % Constellation Size
k = log2(M); % No. of Bits
N = 10^5; % No. of Symbols
fs = 100; % Sampling Freq
fc = 10; % Carrier Freq
Tb = 1; % Symbol Duration
fd = 1/Tb; % Freq Deviation
t = 0:1/fs:Tb-1/fs; % Sampling Instants

ip = rand(1,N*k)>0.5; % Generating 0,1 with equal probability
ipBitReshape =  reshape(ip,k,N).'; % Grouping to N symbols having k bits each
bin2DecMatrix = ones(N,1)*(2.^((k-1):-1:0)) ; % LUT for Conversion
ipDec = sum(ipBitReshape.*bin2DecMatrix,2).'; % Decimal to binary conversion

x = [];
for j = 1:length(ipDec)
    if mod(ipDec(j),2)
        sinfreq = fc + fd*(ipDec(j)-1)/2; % Map constellation to freq
        finsin = repelem(sinfreq,fs);
        xs = (sqrt(2/fs))*sin(2*pi*finsin.*t); % Generating the sin modulated signal
        x(:,j) = xs;
    else
        cosfreq = fc + fd*ipDec(j)/2; % Map constellation to freq
        fincos = repelem(cosfreq,fs);
        xc = (sqrt(2/fs))*cos(2*pi*fincos.*t); % Generating the cos modulated signal
        x(:,j) = xc;
    end
end
x = x(:); % Combine the signals

p = abs(fft((x)))/fs; % DFT
psd = (p(1:fs*N/2).^2)'; % First half
psd1 = [0 psd];
psd2 = [0 0 psd];
psd3 = [0 0 0 psd];
psd4 = [0 0 0 0 psd];

% Moving average filter
psdn = [0.2*psd 0 0 0 0] + [0.2*psd1 0 0 0] + [0.2*psd2 0 0] + [0.2*psd3 0] + 0.2*psd4;

f = (0:fs*N/2+3)/N;
semilogy(f, psdn);
grid on
title('Power spectal density from FFT')
xlabel('Frequency')
ylabel('Power in dB')