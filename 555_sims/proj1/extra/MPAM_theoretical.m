clear; clc; close all;
load mpam_theo.mat

% Pbe = 2*(M-1)/M*qfunc(sqrt(6*log2(M)/(M^2-1)*EbN0)); % calculated
% Pse = 2*(M-1)/M*qfunc(sqrt(6/(M^2-1)*EsN0));

% CALCULATIONS
    EbN0 = 10.^(EbN0dB/10);

    % 2PAM
    ber1 = ber_2pam_theo; % produced by Matlab Communications Toolbox
    M = 2;
    ber2 = 2*(M-1)/M*qfunc(sqrt(6*log2(M)/(M^2-1)*EbN0))/log2(M); % calculated
    
    % 4PAM
    ber3 = ber_4pam_theo; % produced by Matlab Communications Toolbox
    M = 4;
    ber4 = 2*(M-1)/M*qfunc(sqrt(6*log2(M)/(M^2-1)*EbN0))/log2(M); % calculated
    
    % 8PAM
    ber5 = ber_8pam_theo; % produced by Matlab Communications Toolbox
    M = 8;
    ber6 = 2*(M-1)/M*qfunc(sqrt(6*log2(M)/(M^2-1)*EbN0))/log2(M); % calculated
    
    % 16PAM
    ber7 = ber_16pam_theo; % produced by Matlab Communications Toolbox
    M = 16;
    ber8 = 2*(M-1)/M*qfunc(sqrt(6*log2(M)/(M^2-1)*EbN0))/log2(M); % calculated
    
    clear M; %bits;

    
    
if(1)
% PLOTS
    % 2PAM
    figure(1);
    subplot(2,2,1);
    semilogy(EbN0dB, ber1, '-o', ...
             EbN0dB, ber2);
    axis([0 20 10^-8 10^0]);
    legend('2PAM Matlab','2PAM');
    ylabel('BER'); xlabel('Eb/N0 (dB)');
    title('2PAM (theoretical): BER vs. Eb/N0');


    % 4PAM
    % figure(2);
    subplot(2,2,2);
    semilogy(EbN0dB, ber3, '-o', ...
             EbN0dB, ber4);
    axis([0 20 10^-8 10^0]);
    legend('4PAM Matlab','4PAM');
    ylabel('BER'); xlabel('Eb/N0 (dB)');
    title('4PAM (theoretical): BER vs. Eb/N0');


    % 8PAM
    % figure(3);
    subplot(2,2,3);
    semilogy(EbN0dB, ber5, '-o', ...
             EbN0dB, ber6);
    axis([0 20 10^-8 10^0]);
    legend('8PAM Matlab','8PAM');
    ylabel('BER'); xlabel('Eb/N0 (dB)');
    title('8PAM (theoretical): BER vs. Eb/N0');


    % 16PAM
    % figure(4);
    subplot(2,2,4);
    semilogy(EbN0dB, ber7, '-o', ...
             EbN0dB, ber8);
    axis([0 20 10^-8 10^0]);
    legend('16PAM Matlab','16PAM');
    ylabel('BER'); xlabel('Eb/N0 (dB)');
    title('16PAM (theoretical): BER vs. Eb/N0');


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
    axis([0 20 10^-8 10^0]);
    legend('2PAM','4PAM','8PAM','16PAM');
    ylabel('BER'); xlabel('Eb/N0 (dB)');
    title('BER vs. Eb/N0 (theoretical)');
end

%     % Test
%     figure(3);
%     semilogy(EbN0dB, ber5, '-o', ...
%              EbN0dB, ber6);
%     axis([0 20 10^-8 10^0]);
