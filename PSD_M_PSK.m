close all
clear
clc

N = 10^4; % No.of symbols
M = 16; % Contellation size
k = log2(M); % No. of bits

T = 1; % Time observed
fs = 5*10^2; % Sampling freq
t = 0:1/fs:(T-1/fs);
fc = 200; % Carrier Frequency

% Mapping for binary <--> Gray code conversion
ref = 0:M-1;
map = bitxor(ref,floor(ref/2));
[~,ind] = sort(map);

ipBit = rand(1,N*k)>0.5; % random 1's and 0's
bin2DecMatrix = ones(N,1)*(2.^((k-1):-1:0)) ; % conversion from binary to decimal
ipBitReshape =  reshape(ipBit,k,N).'; % grouping to N symbols having k bits each
ipGray = sum(ipBitReshape.*bin2DecMatrix,2).'; % decimal to binary

% Gray coded constellation mapping
ipDec = ind(ipGray+1)-1; % bit group to constellation point
ipPhase = (ipDec*2*pi/M)'; % conversion to phase

sm = cos(ipPhase)*cos(2*pi*fc*t) - sin(ipPhase)*sin(2*pi*fc*t);
smt = sm';
fsm = smt(:);
msm = mean(sm);

p = abs(fft((fsm)))*2/fs;
psd = (p(1:fs*N/2).^2)';
psd1 = [0 psd];
psd2 = [0 0 psd];
psd3 = [0 0 0 psd];
psd4 = [0 0 0 0 psd];

% moving average
psdn = [0.2*psd 0 0 0 0] + [0.2*psd1 0 0 0] + [0.2*psd2 0 0] + [0.2*psd3 0] + 0.2*psd4;

f = (0:fs*N/2+3)/N;
plot(f, psdn);
grid on
% axis([190 210 0.5 10^5])
title('Power spectal density from FFT')
xlabel('Frequency')
ylabel('Power in dB')
