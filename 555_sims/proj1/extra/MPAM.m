close all; clear; clc;



% SIMULATION PARAMETERS
order = 3;
M = 2;              % # of symbols: choose 2, 4, 8, or 16
EbN0dB = 0:20;      % select range of Eb/N0
              
k = log2(M);        % bits per symbol
K = k*10^order;     % # of instances (memory limitation)
% K = 8;            % 4PAM test
% K = 30;           % 8PAM test
% K = 64;           % 16PAM test
num_syms = K/k;     % # of symbols

EbN0 = 10.^(EbN0dB/10);
%N0 = 1;
%Eb = EbN0*N0;
%Es = k*Eb;
%EsN0 = Es/N0;
EsN0dB = 10*log10(EbN0*k);
SNR_rng = length(EbN0dB);

switch M
    case 2      % 2ASK symbols
        sym_map = [-1,1];
    case 4      % 4ASK symbols
        sym_map = [-3/sqrt(5) -1/sqrt(5) 3/sqrt(5) 1/sqrt(5)];
%     case 8      % 8ASK symbols
%     case 16     % 16ASK symbols
end


% SIMULATION
BER = zeros(1,SNR_rng); SER = zeros(1,SNR_rng);
for i = 1:length(EbN0) % test each SNR
    
    % generate bit stream
    b = (sign(randn(1,K))+1)/2;
    
    % test bit streams (gray coded)
%     b = [0,0, 1,0, 1,1, 0,1];                             % 4PAM test: [lsb,msb]
%     b = [0,0,0, 1,0,0, 1,0,1, 0,0,1, 0,1,1, ...           % 8PSK test: [lsb,x,msb]
%          1,1,1, 1,1,0, 0,1,0, 0,0,0, 0,0,1];
%     b = [0,0,0,0, 1,0,0,0, 1,1,0,0, 0,1,0,0, ...          % 16PSK test: [lsb,x,x,msb]
%          0,1,1,0, 1,1,1,0, 1,0,1,0, 0,0,1,0, ...
%          0,0,1,1, 1,0,1,1, 1,1,1,1, 0,1,1,1, ...
%          0,1,0,1, 1,1,0,1, 1,0,0,1, 0,0,0,1];

    % transmit
    symbols = reshape(b, k, []);
    sym_enum = 2.^(0:k-1)*symbols+1;  % select symbol (uniquely identify symbol using binary #s)
    sym_enum_xref = mpam_xref(sym_enum,M);   % gray code cross referenced to match symbol map
    s_unit = sym_map(sym_enum_xref);	% unit magnitude signal to be sent
%     s = sqrt(Es(i))*s_unit;             % apply tx energy
    s = s_unit;
    
    No = 10^(-EsN0dB(i)/10);
%     % Random bits
%     B1g=round(rand); B2g=round(rand);       % Gray mapping 
%     % Amplitudes (symbols)             
%     A_gra = S_gra(2*B1g+B2g+1);             % Gray mapping
    
    % convert to real signal -> mix up -> transmit over channel
    
    % channel
    % n = sqrt(N0/2)*randn(1,num_syms)*1;         % 1-dimension
    n = sqrt(No/2)*randn(1,num_syms)*1;         % 1-dimension

    % receive -> mix down -> convert to I-Q
    r = s + n;

    % MAP decision [1]
%     sd = zeros(1,num_syms);
%     for a = 1:num_syms
% %         [ee ind] = min(abs(sym_map*sqrt(Es(i))-r(a)));      % find minimum distance [ee=min. dist. value, ind=index of min. dist.]
%         [ee ind] = min(abs(sym_map-r(a)));      % find minimum distance [ee=min. dist. value, ind=index of min. dist.]
%         sd(1,a) = round(10*sym_map(ind));                             % symbol decision
%     end
    
            % Decision device: Natural
     if (M==2)
         sd = (sign(r)+1)/2;
     else
         sd = zeros(2,num_syms);
         for a = 1:num_syms
            if r(a) < -2/sqrt(5),                        sd(1,a)=0; sd(2,a)=0;
            elseif ((r(a)>-2/sqrt(5)) && (r(a)<=0)),     sd(1,a)=1; sd(2,a)=0;
            elseif ((r(a)>0) && (r(a)<=2/sqrt(5))),      sd(1,a)=0; sd(2,a)=1;
            else sd(1,a)=1; sd(2,a)=1;
            end
         end
     end
     
    s_est = sd;

    
%     s_est = mpam_map(sd,M);
    
    % symbol error determination [1]
    sym_err_v = abs(symbols-s_est);             % a vector of symbol errors (0 = no error, 1+ = error(s) (1=error if M=2))
    if (M~=2)
        sym_err_v = sign(sum(sym_err_v));       % a vector of symbol errors (0 = no error, 1 = error)
    end
    SER(i) = sum(sym_err_v)/num_syms;           % (# of symbol errors)/(# of symbols)
    
    % bit error determination [1]
    b_est = reshape(s_est, 1, []);              % reshape symbol matrix into bit stream
    b_err_v = abs(b-b_est);                     % a vector of bit errors (0 = no error, 1 = error)
    BER(i) = sum(b_err_v)/(k*num_syms);         % sums all bit errors w/in symbols, then sums up all those values for total bit erro
 
%     if (BER(i) < 100/K) break; end

end

if (order==1 || order==2)
    symbols
    s_unit
    sd
    s_est
%     abs(symbols-s_est)
    sym_err_v
end
    No/2
    var(n)


% THEORETICAL CALCULATIONS
% Pb = 2*(M-1)/M*qfunc(sqrt(6*log2(M)/(M^2-1)*EbN0));
% Ps = 2*(M-1)/M*qfunc(sqrt(6/(M^2-1)*EsN0));
% Ps = 2*(M-1)/M*qfunc(sqrt(6*k/(M^2-1)*EbN0));
% Pb = 2*(M-1)/M*qfunc(sqrt(6*k/(M^2-1)*EbN0))/k;
[Pb,Ps] = berawgn(EbN0dB,'pam',M);




% PLOTS [1]
if(1)
figure(1);
    subplot(1,2,1);
    semilogy(EbN0dB,SER,'ok', ...
             EbN0dB,Ps,'g');
    axis([min(EbN0dB) max(EbN0dB) 10^-8 10^0]);
    legend('SER simulated','SER theoretical');
    ylabel('SER'); xlabel('E_b/N_o (dB)');
    title(['Symbol Error Rate, M=' int2str(M)]);
   
subplot(1,2,2);
semilogy(EbN0dB,BER,'ok', ...
         EbN0dB,Pb,'-g');
axis([min(EbN0dB) max(EbN0dB) 10^-8 10^0]);
legend('BER simulated','BER theoretical');
ylabel('BER'); xlabel('E_b/N_o (dB)');
title(['Bit Error Rate, M=' int2str(M)]);

figure(2);
scatter(real(r), imag(r)); % final SNR
end



% [1] MPSK, by K. Bell (11/22/99)
