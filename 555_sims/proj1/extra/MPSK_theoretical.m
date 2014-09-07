clear; clc; close all;
load mpsk_theo.mat



% CALCULATIONS
    EbN0 = 10.^(EbN0dB/10);

    % BPSK
    ber1 = BER_BPSK; % produced by Matlab Communications Toolbox
    ber2 = qfunc(sqrt(2*EbN0)); % calculated
    
    % QPSK
    ber3 = BER_QPSK; % produced by Matlab Communications Toolbox
    ber4 = qfunc(sqrt(2*EbN0)); % calculated
    
    % 8PSK
    ber5 = BER_8PSK; % produced by Matlab Communications Toolbox
    M = 8; bits = log2(M);
    ber6 = 2*qfunc(sqrt(2*EbN0*log2(M))*sin(pi/M))/bits; % calculated
    
    % 16PSK
    ber7 = BER_16PSK; % produced by Matlab Communications Toolbox
    M = 16; bits = log2(M);
    ber8 = 2*qfunc(sqrt(2*EbN0*log2(M))*sin(pi/M))/bits; % calculated
    
    clear M bits;

    
    
if(1)
% PLOTS
    % BPSK
    figure(1);
    subplot(2,2,1);
    semilogy(EbN0dB, ber1, '-o', ...
             EbN0dB, ber2);
    axis([0 18 10^-8 10^0]);
    legend('BPSK Matlab','BPSK');
    ylabel('BER'); xlabel('Eb/N0 (dB)');
    title('BPSK (theoretical): BER vs. Eb/N0');


    % QPSK
    % figure(2);
    subplot(2,2,2);
    semilogy(EbN0dB, ber3, '-o', ...
             EbN0dB, ber4);
    axis([0 18 10^-8 10^0]);
    legend('QPSK Matlab','QPSK');
    ylabel('BER'); xlabel('Eb/N0 (dB)');
    title('QPSK (theoretical): BER vs. Eb/N0');


    % 8PSK
    % figure(3);
    subplot(2,2,3);
    semilogy(EbN0dB, ber5, '-o', ...
             EbN0dB, ber6);
    axis([0 18 10^-8 10^0]);
    legend('8PSK Matlab','8PSK');
    ylabel('BER'); xlabel('Eb/N0 (dB)');
    title('8PSK (theoretical): BER vs. Eb/N0');


    % 16PSK
    % figure(4);
    subplot(2,2,4);
    semilogy(EbN0dB, ber7, '-o', ...
             EbN0dB, ber8);
    axis([0 18 10^-8 10^0]);
    legend('16PSK Matlab','16PSK');
    ylabel('BER'); xlabel('Eb/N0 (dB)');
    title('16PSK (theoretical): BER vs. Eb/N0');


    % All
    figure(2);
    semilogy(EbN0dB, ber2, '-o', ...
             EbN0dB, ber4, ...
             EbN0dB, ber6, ...
             EbN0dB, ber8);
%     semilogy(EbN0dB, ber1, '-o', ...
%              EbN0dB, ber3, ...
%              EbN0dB, ber5, ...
%              EbN0dB, ber7);
    axis([0 18 10^-8 10^0]);
    legend('BPSK','QPSK','8PSK','16PSK');
    ylabel('BER'); xlabel('Eb/N0 (dB)');
    title('BER vs. Eb/N0 (theoretical)');
end

%     % Test
%     figure(3);
%     semilogy(EbN0dB, ber5, '-o', ...
%              EbN0dB, ber6);
%     axis([0 18 10^-8 10^0]);
