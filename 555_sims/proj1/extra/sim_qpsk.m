close all; clear; clc;

K = 10^6; % # of instances (memory limitation)

% select transmission and channel parameters
EbN0dB=0:18;
EbN0=10.^(EbN0dB/10);
N0=1;
Eb=EbN0*N0;
Es=Eb;

BER = zeros(1,length(Es));
for i = 1:length(Es) % test each SNR
    
    % generate bit stream
    b = (sign(randn(1, K))+1)/2;

    % transmit
    sym_alph = [1+1j, -1+1j, 1-1j, -1-1j]/sqrt(2); % symbol alphabet
    bits = log2(length(sym_alph)); % bits per symbol
    symbols = reshape(b, bits, []);
    sym_sel = 2.^[0:bits-1]*symbols; % select symbol (uniquely identify symbol using binary #s)
    s = sqrt(Es(i))*sym_alph(sym_sel+1);
    % convert to real signal -> mix up -> transmit over channel
    
    % channel
    n = sqrt(N0/2)*(randn(1,K/bits) + 1j*randn(1,K/bits))/sqrt(2);

    % receive -> mix down -> convert to I-Q
    r = s + n;

    % detection (MAP)
    rx_symbols = [real(r)<0; imag(r)<0]*1.0;
    b_est = reshape(rx_symbols, 1, [])*1.0;
    
    % calculate BER
    bit_errors = sum(b ~= b_est)*1.0;
    BER(i) = bit_errors/K;

    if (BER(i) < 100/K) break; end

end

% figure(1);
% scatter(real(r), imag(r)); % final SNR

ber_qpsk = qfunc(sqrt(2*EbN0));
ser_qpsk = 2*ber_qpsk;

figure(2);
semilogy(EbN0dB, ber_qpsk, ...
         EbN0dB, ser_qpsk, ...
         EbN0dB,BER, '-o' ...
         );
axis([0 18 10^-8 10^0]);
legend('BER QPSK theo.','SER QPSK theo.','QPSK sim.');
