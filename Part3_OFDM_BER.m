close all
clear
clc

M = 4; % Constellation Size
k = log2(M); % No. of Bits
N = 10^5; % No. of Symbols
fs = 100; % Sampling Freq
fc = 20; % Carrier Freq
Tb = 1; % Symbol Duration
fd = 1/Tb; % Freq Deviation
t = 0:1/fs:Tb-1/fs; % Sampling Instants
Eb_N0_dB = 0:15; % Range of SNR values
Es_N0_dB  = Eb_N0_dB + 10*log10(k);

for ii = 1:length(Eb_N0_dB)
    disp(['Running Sim for SNR = ',num2str(ii)])
    
    ip = rand(1,N*k)>0.5; % Generating 0,1 with equal probability
    ipBitReshape =  reshape(ip,k,N).'; % Grouping to N symbols having k bits each
    bin2DecMatrix = ones(N,1)*(2.^((k-1):-1:0)) ; % LUT for Conversion
    ipDec = sum(ipBitReshape.*bin2DecMatrix,2).'; % Decimal to binary conversion
    x = [];
    for j = 1:length(ipDec)
        if mod(ipDec(j),2)
            sinfreq = fc + fd*(ipDec(j)-1)/2; 
            finsin = repelem(sinfreq,fs);
            xs = (sqrt(2/fs))*sin(2*pi*finsin.*t); % Generating the sin modulated signal
            x(:,j) = xs;
        else
            cosfreq = fc + fd*ipDec(j)/2; % Constellation to freq
            fincos = repelem(cosfreq,fs);
            xc = (sqrt(2/fs))*cos(2*pi*fincos.*t); % Generating the cos modulated signal
            x(:,j) = xc;
        end
    end
    x = x(:)'; % Combine
    
    % Adding noise
    n = 1/sqrt(2)*(randn(1,N*fs) + 1i*randn(1,N*fs)); % White guassian noise, 0dB variance1
    y = x + 10^(-Es_N0_dB(ii)/20)*n; % Additive white gaussian noise
    
    % Demodulation
    fop = [];
    for j = 1:M/2
        opc(:,j) = conv(y, sqrt(2/fs)*cos(2*pi*(fc+fd*(j-1))*t)); % correlating with cos
        fop = [fop;opc(fs+1:fs:end,j)]; % Sample at t = T
        ops(:,j) = conv(y, sqrt(2/fs)*sin(2*pi*(fc+fd*(j-1))*t)); % correlating with sin
        fop = [fop;ops(fs+1:fs:end,j)]; % Sample at t = T
    end
    
    fop = reshape(fop,N,M);
    [ff,opDe] = max(abs(fop),[],2); % Choose max
    opDet = (opDe - 1)';
    % Convert Decimal to binary
    opBinHat = dec2bin(opDet);
    opBinHat = opBinHat.';
    opBinHat = opBinHat(1:end).';
    opBinHat = reshape(str2num(opBinHat).',k,N).' ;
    
    nBitErr(ii) = size(find((ipBitReshape - opBinHat)),1); % Counting the number of errors
end

simBer = nBitErr/(k*N) + 1e-100; % Adding small value to prevent -ve infinity
theoryBer = M/2*0.5*erfc(sqrt((10.^(Eb_N0_dB/10)/2))); % Theoretical BER    

figure
semilogy(Eb_N0_dB,theoryBer,'bs-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer,'mx-','LineWidth',2);
axis([0 16 10^-5 0.5])
grid on
legend(['Theoretical BER, M = ',num2str(M)], ['Simulated BER, M = ',num2str(M)]);
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve')