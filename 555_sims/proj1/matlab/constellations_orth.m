close all; clear; clc;


% SIMULATION PARAMETERS
order = 4;
M = 2^4;              % # of symbols: choose 2, 4, 8, or 16
sch = 3;            % modulation scheme: 0=PAM, 1=PSK, 2=QAM, 3=Orthogonal
nt = 1;             % normal or test mode: 0=test, 1=normal
plots = 1;          % open plots? 0=no, 1=yes
EbN0dB = 0:25;      % select range of Eb/N0
              
k = log2(M);        % bits per symbol
if(nt)
    K = k*10^order;	% # of instances (memory limitation)
else                % test
    switch M
        case 2,  K = k*10^order;
        case 4,  K = 8;     % 4-ary test
        case 8,  K = 24;	% 8-ary test
        case 16, K = 64;	% 16-ary test
    end
end
num_syms = K/k;     % # of symbols
EbN0 = 10.^(EbN0dB/10);
EsN0dB = 10*log10(EbN0*k);
rng = length(EbN0dB);
sym_tbl_o = sym_table_orth(M);


% SIMULATION
BER = zeros(1,rng); SER = zeros(1,rng);
for i = 1:rng % test each SNR
    
    % generate bit stream
    if (nt), b = (sign(randn(1,K))+1)/2;
    else     b = test_data(M); end

    % transmit
    symbols = reshape(b, k, []);
    sym_ind = (2.^(0:k-1)*symbols)+1;       % select symbol (uniquely identify symbol using binary #s)
    sym_ind_gray = gray_code(sym_ind,M);    % gray code cross referenced to match symbol map
    s = sym_tbl_o(:,sym_ind_gray);              % unit magnitude signal to be sent
    
    % channel
    N0 = 10^(-EsN0dB(i)/10);
    n = sqrt(N0/2)*randn(size(s));

    % receive
    r = s + n;

    % MAP decision [1]
    sym_err = 0; temp = 0;
    s_est = zeros(size(symbols));  test = zeros(size(sym_tbl_o));
    for col = 1:num_syms
        for row = 1:size(r,1); test(row,:) = abs(sym_tbl_o(row,:)-r(row,col)); end
        if (M~=2&&sch~=0||sch==3), test = sum(test); end
        [ee ind] = min(test);
        s_est(:,col) = MAP(ind,M);
        
        % terminate early
        temp = abs(symbols(:,col)-s_est(:,col));
        if(M~=2), temp = sum(temp); end
        if(temp>0), sym_err=sym_err+1; end
        if (sym_err>=100), s_est=s_est(:,1:col); break; end
    end
    
    % symbol error determination [1]
    sym_err_v = abs(symbols(:,1:col)-s_est);  % a vector of symbol errors (0 = no error, 1+ = error(s) (1=error if M=2))
    if (M~=2), sym_err_v = sign(sum(sym_err_v)); end	% a vector of symbol errors (0 = no error, 1 = error)
    SER(i) = sum(sym_err_v)/col;
    
    % bit error determination [1]
    b_est = reshape(s_est,1,[]);                        % reshape symbol matrix into bit stream
    b_err_v = abs(b(1,1:k*col)-b_est);                    % a vector of bit errors (0 = no error, 1 = error)
    BER(i) = sum(b_err_v)/(k*col);
 
%     if (BER(i) < 100/K) break; end
end

if (order==1 || order==2)
    symbols
    sym_ind_gray
    s
    r
    s_est
    sym_err_v % test
end
%     N0/2
%     var(n)


% THEORETICAL CALCULATIONS

%            Pb = qfunc(sqrt(2*EbN0*k/M));
%            Ps = (2^k-1)/2^(k-1)*Pb;
%            Ps = 1-(1-qfunc(sqrt(2*EbN0))).^k;
% %            Pb = Ps/k;
%            Pb = 2^(k-1)/(2^k-1)*Ps;
            [Pb,Ps] = berawgn(EbN0dB,'fsk',M,'coherent');
           
%            Pb_u = M*exp(-EbN0*k/2);
%            [Pb1,Ps1] = berawgn(EbN0dB,'psk',M,'nondiff');
%            Ps1 = 2*exp(-log2(M)*(0.5*EbN0-log(2)));
%            Pb1 = Ps1/k;

% upper bounds
%            Ps_u = exp(-k*(EbN0-2*log(2))/2); % (loose)
           Ps_u = (M-1)*qfunc(sqrt(k*EbN0)); % coherent (tight)
           Pb_u = 2^(k-1)/(2^k-1)*Ps_u;%(M/2)/(M-1)*Ps_u;
%             Ps_u = exp(-k*(EbN0-2*log(2))/2);
%             Pb_u = Ps_u/k;


% PLOTS [1]
if(plots)
%     figure(1);
%     if (M==2||sch==0), scatter(r,zeros(1,length(r)));
%     else scatter(r(1,:),r(2,:)); end
%     axis([-2 2 -2 2]);

    figure(2);
        subplot(1,2,1);
        semilogy(EbN0dB,SER,'-ok', ...
                 EbN0dB,Ps,'g', ...
                 EbN0dB,Ps_u,'r');
                 
        axis([min(EbN0dB) max(EbN0dB) 10^-6 10^-1]);
        legend('SER simulated','SER FSK theoretical','SER upper bound');
        ylabel('SER'); xlabel('E_b/N_o (dB)');
        title(['Symbol Error Rate, M=' int2str(M)]);

    subplot(1,2,2);
    semilogy(EbN0dB,BER,'-ok', ...
             EbN0dB,Pb,'g', ...
             EbN0dB,Pb_u,'r');        

    axis([min(EbN0dB) max(EbN0dB) 10^-6 10^-1]);
    legend('BER simulated','BER FSK theoretical','BER upper bound');
    ylabel('BER'); xlabel('E_b/N_o (dB)');
    title(['Bit Error Rate, M=' int2str(M)]);
end


% [1] MPSK, by K. Bell (11/22/99)
