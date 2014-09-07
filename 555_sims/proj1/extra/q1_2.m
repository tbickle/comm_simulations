close all; clear; clc;

K = 10^4; % # of instances (memory limitation)

% select transmission and channel parameters
EbN0dBv=0:8;
EbN0v=10.^(EbN0dBv/10);
N0=1;
Ebv=EbN0v*N0;
Esv=Ebv;

BER = zeros(1,length(Esv));
for i = 1:length(Esv) % test each SNR
    % generate bit stream
    b = (sign(randn(1, K))+1)/2;

    % transmit
    s = sqrt(Esv(i))*(2*b-1);

    % channel
    n = sqrt(N0/2)*randn(1,K);

    % receive
    r = s + n;

    % detection (MAP)
    b_est = (sign(r)+1)/2;

    % BER
    BER(i) = sum(abs(b - b_est))/K;
    
    if (BER(i) < 100/K) break; end
end

figure(1);
semilogy(EbN0dBv,BER,'-o',EbN0dBv,qfunc(sqrt(2*EbN0v)));
