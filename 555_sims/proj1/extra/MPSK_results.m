close all; clear; clc;
load mpsk_results.mat

EbN0 = 10.^(EbN0dB/10);

% THEORETICAL CALCULATIONS
ber_bpsk_theo = qfunc(sqrt(2*EbN0));
ser_bpsk_theo = ber_bpsk_theo;
ber_qpsk_theo = qfunc(sqrt(2*EbN0));
ser_qpsk_theo = 1-(1-ber_qpsk_theo).^2;
ser_8psk_theo = 2*qfunc(sqrt(2*EbN0*log2(8))*sin(pi/8)); % calculated
ber_8psk_theo = ser_8psk_theo/log2(8);
ser_16psk_theo = 2*qfunc(sqrt(2*EbN0*log2(16))*sin(pi/16)); % calculated
ber_16psk_theo = ser_16psk_theo/log2(16);

% PLOTS
figure(1);
subplot(1,2,1);
semilogy(EbN0dB*log2(2), ser_bpsk_sim, 'ok', ...
         EbN0dB*log2(2), ser_bpsk_theo, 'k', ...
         EbN0dB*log2(4), ser_qpsk_sim, 'or', ...
         EbN0dB*log2(4), ser_qpsk_theo, 'r', ...
         EbN0dB*log2(8), ser_8psk_sim, 'og', ...
         EbN0dB*log2(8), ser_8psk_theo, 'g', ...
         EbN0dB*log2(16), ser_16psk_sim, 'ob', ...
         EbN0dB*log2(16), ser_16psk_theo, 'b');
axis([min(EbN0dB*log2(2)) max(EbN0dB*log2(2)) 10^-8 10^0]);
% legend('SER simulated','SER theoretical');
ylabel('SER'); xlabel('E_s/N_o (dB)');
title('Symbol Error Rate');
   
subplot(1,2,2);
semilogy(EbN0dB, ber_bpsk_sim, 'ok', ...
         EbN0dB, ber_qpsk_sim, 'or', ...
         EbN0dB, ber_bpsk_theo, 'k', ... % ber_bpsk_theo = ber_qpsk_theo
         EbN0dB, ber_8psk_sim, 'og', ...
         EbN0dB, ber_8psk_theo, 'g', ...
         EbN0dB, ber_16psk_sim, 'ob', ...
         EbN0dB, ber_16psk_theo);
axis([min(EbN0dB) max(EbN0dB) 10^-8 10^0]);
% legend('BER simulated','BER theoretical');
ylabel('BER'); xlabel('E_b/N_o (dB)');
title('Bit Error Rate');
