% Bit Error Rate for 16-PSK modulation using Gray modulation mapping
close all
clear
clc

N = 10^5; % number of symbols
M = 2;   % constellation size
k = log2(M); % bits per symbol

Eb_N0_dB  = 0:25; % multiple Es/N0 values
Es_N0_dB  = Eb_N0_dB + 10*log10(k);

% Mapping for binary <--> Gray code conversion
ref = 0:M-1;
map = bitxor(ref,floor(ref/2));
[tt, ind] = sort(map);

estPhase = zeros(1,N);
for ii = 1:length(Eb_N0_dB)
    
    % symbol generation
    ipBit = rand(1,N*k,1)>0.5; % random 1's and 0's
    ipBitReshape =  reshape(ipBit,k,N).'; % grouping to N symbols having k bits each
    bin2DecMatrix = ones(N,1)*(2.^((k-1):-1:0)); % conversion from binary to decimal
    grayRef = sum(ipBitReshape.*bin2DecMatrix,2).'; % decimal to binary
    
    % Gray coded constellation mapping
    ipDec = ind(grayRef+1)-1; % bit group to constellation point
    ipPhase = ipDec*2*pi/M; % conversion to phase
    s = exp(1i*ipPhase);  % modulation
    
    % adding noise
    n = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N)); % white guassian noise, 0dB variance1
    y = s + 10^(-Es_N0_dB(ii)/20)*n; % additive white gaussian noise
    
    % demodulation
    opPhase = angle(y);
    % unwrap phase
    opPhase(opPhase<0) = opPhase(opPhase<0) + 2*pi;
    % rounding the received phase to the closest constellation
    estPhase = 2*pi/M*round(opPhase/(2*pi/M)); % phase equalization
    estPhase(estPhase==2*pi) = 0;
    estDec = round(estPhase*M/(2*pi));
    
    % Decimal to Gray code conversion
    ipGrayHat = map(estDec+1); % converting to decimal
    estBinary = dec2bin(ipGrayHat,k) ; % decimal to binary
    
    % converting binary string to number
    estBinary = estBinary.';
    estBinary = estBinary(1:end).';
    estBinary = str2num(estBinary).' ;
    
    % counting errors
    nBitErr(ii) = size(find(ipBit- estBinary),2); % couting the number of errors
end

simBer = nBitErr/(N*k) + 1e-100; % small term added to prevent -infinty in log

if(M == 8 || M == 16)
    theoryBer = (1/k)*erfc(sqrt(k*10.^(Eb_N0_dB/10))*sin(pi/M));
else
    theoryBer = (1/2)*erfc(sqrt(10.^(Eb_N0_dB/10)));
end

figure
semilogy(Eb_N0_dB,theoryBer,'bs-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer,'mx-','LineWidth',2);
axis([0 20 10^-5 1])
grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for 16-PSK modulation')
