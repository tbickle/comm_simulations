close all; clear; clc;



% SIMULATION PARAMETERS
order = 4;
M = 2;              % # of symbols: choose 2, 4, 8, or 16
EbN0dB = 0:20;      % select range of Eb/N0
              
k = log2(M);        % bits per symbol
K = k*10^order;     % # of instances (memory limitation)
% K = 8;            % QPSK test
% K = 24;           % 8PSK test
% K = 64;           % 16PSK test
num_syms = K/k;     % # of symbols

EbN0 = 10.^(EbN0dB/10);
EsN0dB = 10*log10(EbN0*k);
rng = length(EbN0dB);

if(0) % 0=PAM, 1=PSK
    switch M
        case 2      % BPSK symbols
            sym_map = [-1, 1];
        case 4      % QPSK symbols
            sym_map = [1, 1j, -1, -1j];
        case 8      % 8PSK symbols
            sym_map = [1, (1+1j)/sqrt(2), 1j, (-1+1j)/sqrt(2), -1, (-1-1j)/sqrt(2), -1j, (1-1j)/sqrt(2)]; % gray code
        case 16     % 16PSK symbols
            sym_map = [cos(0*pi/8)+1j*sin(0*pi/8), cos(1*pi/8)+1j*sin(1*pi/8), cos(2*pi/8)+1j*sin(2*pi/8), cos(3*pi/8)+1j*sin(3*pi/8), ...
                       cos(4*pi/8)+1j*sin(4*pi/8), cos(5*pi/8)+1j*sin(5*pi/8), cos(6*pi/8)+1j*sin(6*pi/8), cos(7*pi/8)+1j*sin(7*pi/8), ...
                       cos(8*pi/8)+1j*sin(8*pi/8), cos(9*pi/8)+1j*sin(9*pi/8), cos(10*pi/8)+1j*sin(10*pi/8), cos(11*pi/8)+1j*sin(11*pi/8), ...
                       cos(12*pi/8)+1j*sin(12*pi/8), cos(13*pi/8)+1j*sin(13*pi/8), cos(14*pi/8)+1j*sin(14*pi/8), cos(15*pi/8)+1j*sin(15*pi/8)];
    end
else
    switch M
    case 2      % 2ASK symbols
        sym_map = [-1,1];
    case 4      % 4ASK symbols
        sym_map = [-3/sqrt(5) -1/sqrt(5) 3/sqrt(5) 1/sqrt(5)];
    end
end



% SIMULATION
BER = zeros(1,rng); SER = zeros(1,rng);
for i = 1:rng % test each SNR
%     EbN0dB(i)
    
    % generate bit stream
    b = (sign(randn(1,K))+1)/2;
%     b = test_data(M);

    % transmit
    symbols = reshape(b, k, []);
    sym_ind = (2.^(0:k-1)*symbols)+1;       % select symbol (uniquely identify symbol using binary #s)
    sym_ind_gray = gray_code(sym_ind,M);    % gray code cross referenced to match symbol map
    s = sym_map(sym_ind_gray);              % unit magnitude signal to be sent
    % convert to real signal -> mix up -> transmit over channel
    
    % channel
    N0 = 10^(-EsN0dB(i)/10);    
    if (M==2), n = sqrt(N0/2)*randn(1,num_syms)*1;                                  % 1-dimension
    else       n = sqrt(N0/2)*(randn(1,num_syms) + 1j*randn(1,num_syms))*1; end     % 2-dimension

    % receive -> mix down -> convert to I-Q
    r = s + n;

    
    
    
    
    % MAP decision [1]
    s_est = zeros(k,num_syms); sd = zeros(1,num_syms);
    if(0) % 0=PAM, 1=PSK
        for a = 1:num_syms
            [ee ind]=min(abs(sym_map-r(a)));        % find minimum distance [ee=min. dist. value, ind=index of min. dist.]
            sd(1,a)= ind;                           % symbol index decision
            s_est(:,a) = MAP(ind,M);
        end
    else
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
    end
    
    
    
    
    
    % symbol error determination [1]
    sym_err_v = abs(symbols-s_est);                     % a vector of symbol errors (0 = no error, 1+ = error(s) (1=error if M=2))
    if (M~=2), sym_err_v = sign(sum(sym_err_v)); end	% a vector of symbol errors (0 = no error, 1 = error)
    SER(i) = sum(sym_err_v)/num_syms;                   % (# of symbol errors)/(# of symbols)
    
    % bit error determination [1]
    b_est = reshape(s_est, 1, []);              % reshape symbol matrix into bit stream
    b_err_v = abs(b-b_est);                     % a vector of bit errors (0 = no error, 1 = error)
    BER(i) = sum(b_err_v)/(k*num_syms);         % sums all bit errors w/in symbols, then sums up all those values for total bit erro
 
%     if (BER(i) < 100/K) break; end

end

if (order==1 || order==2)
    symbols
    sym_ind_gray
    s
    sd
    s_est
%     abs(symbols-s_est)
    sym_err_v
end
    N0/2
    var(n)



% THEORETICAL CALCULATIONS
if(0) % 0=PAM, 1=PSK
    switch M
        case 2
           Pb = qfunc(sqrt(2*EbN0));
           Ps = Pb;
        case 4
           Pb = qfunc(sqrt(2*EbN0));
           Ps = 1-(1-Pb).^2;
        case 8
           Ps = 2*qfunc(sqrt(2*EbN0*k)*sin(pi/M));
           Pb = Ps/k;
        case 16
           Ps = 2*qfunc(sqrt(2*EbN0*k)*sin(pi/M));
           Pb = Ps/k;
    end
else
    [Pb,Ps] = berawgn(EbN0dB,'pam',M);
end



% PLOTS [1]
if(1)
figure(1);
scatter(real(r), imag(r)); % final SNR

figure(2);
    subplot(1,2,1);
    semilogy(EsN0dB,SER,'ok', ...
             EsN0dB,Ps,'g');
    axis([min(EsN0dB) max(EsN0dB) 10^-8 10^0]);
    legend('SER simulated','SER theoretical');
    ylabel('SER'); xlabel('E_s/N_o (dB)');
    title(['Symbol Error Rate, M=' int2str(M)]);
   
subplot(1,2,2);
semilogy(EbN0dB,BER,'ok', ...
         EbN0dB,Pb,'-g');
axis([min(EbN0dB) max(EbN0dB) 10^-8 10^0]);
legend('BER simulated','BER theoretical');
ylabel('BER'); xlabel('E_b/N_o (dB)');
title(['Bit Error Rate, M=' int2str(M)]);
end



% [1] MPSK, by K. Bell (11/22/99)
