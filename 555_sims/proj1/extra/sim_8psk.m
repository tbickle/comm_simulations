% note: qfunc(sqrt(2*Eb/N0)) = 0.5*erfc(sqrt(2*Eb/N0)/sqrt(2))
close all; clear; clc;

K = 3*10^4; % # of instances (memory limitation)

% select transmission and channel parameters
EbN0dB=0:18;
EbN0=10.^(EbN0dB/10);
N0=2;
Eb=EbN0*N0;
Es=Eb;
M = 8; % # of symbols in alpahbet
bits = log2(M); % bits per symbol
a = 1/sqrt(2); % constant

BER = zeros(1,length(Es));
for i = 1:length(Es) % test each SNR
    
    % generate bit stream
    b = (sign(randn(1, K))+1)/2;
    
       %[lsb,x,msb]
%     b = [0,0,0, 1,0,0, 1,0,1, 0,0,1, 0,1,1, ...
%          1,1,1, 1,1,0, 0,1,0, 0,0,0, 0,0,1];

    % transmit
%     sym_alph = [1, a*(1+1j), 1j, a*(-1+1j), -1, a*(-1-1j), -1j, a*(1-1j)]; % regular
    sym_alph = [1, a*(1+1j), a*(1-1j), -1j, a*(-1+1j), 1j, -1, a*(-1-1j)]; % grey code
    symbols = reshape(b, bits, []);
    sym_sel = 2.^[0:bits-1]*symbols; % select symbol (uniquely identify symbol using binary #s)
%     sym_sel+1
    s = sqrt(Es(i))*sym_alph(sym_sel+1);
    % convert to real signal -> mix up -> transmit over channel
    
    % channel
    n = sqrt(N0/2)*(randn(1,K/bits) + 1j*randn(1,K/bits))*a*1;

    % receive -> mix down -> convert to I-Q
    r = s + n;

    % detection (MAP)
%     symbols
    b_est = MAP(r, bits);

    % calculate BER
    bit_errors = sum(b ~= b_est)*1.0;
    BER(i) = bit_errors/K;

    if (BER(i) < 100/K) break; end

end

% figure(1);
% scatter(real(r), imag(r)); % final SNR
% axis([-1 1 -1 1]);

ber_bpsk = qfunc(sqrt(2*EbN0));
ser_8psk = 2*qfunc(sqrt(2*EbN0*log2(M))*sin(pi/M)); % theoretical
ber_8psk = ser_8psk/bits;
test = 0.5*erfc(sqrt(2*EbN0)/sqrt(2));

figure(2);
semilogy(EbN0dB, ber_bpsk, ...
         EbN0dB, ber_8psk, ...
         EbN0dB, ser_8psk, ...
         EbN0dB, BER,'-o', ...
         EbN0dB, test, '-o' ...
         );
axis([0 18 10^-8 10^0]);
legend('BER BPSK theo.','BER 8PSK theo.','SER 8PSK theo.','8PSK sim.');
