close all; clear; clc;



% SIMULATION PARAMETERS 
order = 5;
M = 8;              % # of symbols
              
k = log2(M);        % bits per symbol
K = k*10^order;     % # of instances (memory limitation) <--------
num_syms = K/k;    % <---------------------------------------

EbN0dB = 0:18;
EbN0 = 10.^(EbN0dB/10);
N0 = 2;
Eb = EbN0*N0;
Es = k*Eb;
EsN0 = Es/N0;
EsN0dB = 10*log10(EsN0);
SNR_rng = length(EbN0dB);

switch M % [1]
case 2      % BPSK symbols
    sym_map = [1,-1];
case 4      % QPSK symbols
    sym_map = [1,1j,-1,-1j];
case 8      % 8PSK symbols
    sym_map = [1,(1+1j)/sqrt(2),(1-1j)/sqrt(2),-1j,(-1+1j)/sqrt(2),1j,-1,(-1-1j)/sqrt(2)];
end


% SIMULATION
BER = zeros(1,SNR_rng); SER = zeros(1,SNR_rng);
for i = 1:length(Es) % test each SNR
    
    % generate bit stream
    b = (sign(randn(1,K))+1)/2;
    
    % test stream: [lsb,x,msb]
%     b = [0,0,0, 1,0,0, 1,0,1, 0,0,1, 0,1,1, 1,1,1, 1,1,0, 0,1,0, 0,0,0, 0,0,1];

    % transmit
    symbols = reshape(b, k, []);
    sym_enum = 2.^(0:k-1)*symbols;  % select symbol (uniquely identify symbol using binary #s)
    s_unit = sym_map(sym_enum+1);   % unit magnitude signal to be sent
    s = sqrt(Es(i))*s_unit;         % apply tx energy
    % convert to real signal -> mix up -> transmit over channel
    
    % channel
    n = sqrt(N0/2)*(randn(1,num_syms) + 1j*randn(1,num_syms));

    % receive -> mix down -> convert to I-Q
    r = s + n;

    % MAP decision [1]
    sd = zeros(2,num_syms);
    for a = 1:num_syms
        [ee ind]=min(abs(sym_map-r(a)));                      % find minimum distance [ee=min. dist. value, ind=index of min. dist.]
        sd(:,a)= [real(sym_map(ind));imag(sym_map(ind))];     % symbol decision
    end
    
    temp = sign(sd);
    s_est = MAP(temp, k);                       % symbol estimation
    
    % symbol error determination [1]
    sym_err_temp = sum(abs(symbols-s_est));     % a vector of symbol errors (0 = no error, 1+ = error(s))
    sym_err_v = sign(sym_err_temp);             % a vector of symbol errors (0 = no error, 1 = error)
    SER(i) = sum(sym_err_v)/num_syms;           % (# of symbol errors)/(# of symbols)
    
    % bit error determination [1]
    b_est = reshape(s_est, 1, []);              % reshape symbol matrix into bit stream
    b_err_v = abs(b-b_est);                     % a vector of bit errors (0 = no error, 1 = error)
    BER(i) = sum(b_err_v)/(k*num_syms);         % sums all bit errors w/in symbols, then sums up all those values for total bit erro
 
    if (BER(i) < 100/K) break; end

end

% symbols
% s_est
% abs(symbols-s_est)



% THEORETICAL CALCULATIONS
switch M
case 2
   Pbe = qfunc(sqrt(2*EbN0));
   Pse = Pbe;
case 4
   Pbe = qfunc(sqrt(2*EbN0));
   Pse = 1-(1-Pbe).^2;
case 8
   Pse = 2*qfunc(sqrt(2*EbN0*k)*sin(pi/M)); % calculated
   Pbe = Pse/k;
end

%lower and upper bounds
% Pslb = 0.5*erfc(sqrt(2*EsN0)*sin(pi/M)/sqrt(2));
% Psub = 2*Pslb;



% PLOTS [1]
figure(1);
    subplot(1,2,1)
    semilogy(EsN0dB,SER,'ok', ...
             EsN0dB,Pse,'g');%, ...
%              10*log10(EsN0),Pslb,'--b', 10*log10(EsN0),Psub,'--b');
    axis([min(EsN0dB) max(EsN0dB) 10^-8 10^0]);
    legend('SER simulated','SER theoretical');
    ylabel('SER'); xlabel('E_s/N_o (dB)');
    title(['Symbol Error Rate, M=' int2str(M)]);
   
subplot(1,2,2);
semilogy(EbN0dB,BER,'ok', ...
         EbN0dB,Pbe,'-g');
axis([min(EbN0dB) max(EbN0dB) 10^-8 10^0]);
legend('BER simulated','BER theoretical');
ylabel('BER'); xlabel('E_b/N_o (dB)');
title(['Bit Error Rate, M=' int2str(M)]);



% [1] MPSK, by K. Bell (11/22/99)
