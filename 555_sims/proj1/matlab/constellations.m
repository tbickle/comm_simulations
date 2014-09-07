close all; clear; clc;


% SIMULATION PARAMETERS
order = 5;
M = 2^4;              % # of symbols: choose 2, 4, 8, or 16
sch = 3;            % modulation scheme: 0=PAM, 1=PSK, 2=QAM
nt = 1;             % normal or test mode: 0=test, 1=normal
plots = 1;          % open plots? 0=no, 1=yes
EbN0dB = 0:25;      % select range of Eb/N0
              
k = log2(M);        % bits per symbol
if(nt)
    K = k*10^order;	% # of instances (memory limitation)
else                % test
    switch M
        case 4,  K = 8;	% 4-ary test
        case 8,  K = 24;	% 8-ary test
        case 16, K = 64;	% 16-ary test
    end
end
num_syms = K/k;     % # of symbols
EbN0 = 10.^(EbN0dB/10);
EsN0dB = 10*log10(EbN0*k);
rng = length(EbN0dB);

sym_tbl = sym_table(M,sch);


% SIMULATION
BER = zeros(1,rng); SER = zeros(1,rng);
for i = 1:rng % test each Eb/N0 increase
    
    % generate bit stream
    if (nt), b = (sign(randn(1,K))+1)/2;
    else     b = test_data(M);
    end

    % transmit
    symbols = reshape(b, k, []);
    sym_ind = (2.^(0:k-1)*symbols)+1;       % select symbol (uniquely identify symbol using binary #s)
    sym_ind_gray = gray_code(sym_ind,M);    % gray code cross referenced to match symbol map
    s = sym_tbl(:,sym_ind_gray);              % unit magnitude signal to be sent
    
    % channel
    N0 = 10^(-EsN0dB(i)/10);
    if (sch~=3),
        if (M==2||sch==0), n = sqrt(N0/2)*randn(1,num_syms)*nt;                                 % 1-dimension
        else               n = sqrt(N0/2)*(randn(1,num_syms) + 1j*randn(1,num_syms))*nt; end	% 2-dimension
    else
        n = sqrt(N0/2)*randn(size(s));
    end

    % receive -> mix down -> convert to I-Q
    r = s + n;

    % MAP decision [1]
    sym_err = 0; temp = 0;
    s_est = zeros(k,num_syms); test = zeros(size(sym_tbl));
    for col = 1:num_syms
        if (sch~=3)           
            [ee ind]=min(abs(sym_tbl-r(col)));        % find minimum distance [ee=min. dist. value, ind=index of min. dist.]
        else
            for row = 1:size(r,1); test(row,:) = abs(sym_tbl(row,:)-r(row,col)); end
            test = sum(test);
            [ee ind] = min(test);                
        end
        s_est(:,col) = MAP(ind,M);

        % terminate early
        temp = abs(symbols(:,col)-s_est(:,col));
        if(M~=2), temp = sum(temp); end
        if(temp>0), sym_err=sym_err+1; end
        if (sym_err>=100), s_est=s_est(:,1:col); break; end
    end
    
    % symbol error determination [1]
    sym_err_v = abs(symbols(:,1:length(s_est))-s_est);  % a vector of symbol errors (0 = no error, 1+ = error(s) (1=error if M=2))
    if (M~=2), sym_err_v = sign(sum(sym_err_v)); end	% a vector of symbol errors (0 = no error, 1 = error)
    SER(i) = sum(sym_err_v)/col;
    
    % bit error determination [1]
    b_est = reshape(s_est,1,[]);                        % reshape symbol matrix into bit stream
    b_err_v = abs(b(1,1:k*col)-b_est);                    % a vector of bit errors (0 = no error, 1 = error)
    BER(i) = sum(b_err_v)/(k*col);
 
    if (BER(i) < 100/K), break; end
end

if (order==1 || order==2)
    symbols
    sym_ind_gray
    s
    s_est
    sym_err_v % test
end
%     N0/2
%     var(n)


% THEORETICAL CALCULATIONS
if(sch == 0)        % PAM
    Ps = 2*(M-1)/M*qfunc(sqrt(6*log2(M)/(M^2-1)*EbN0));
    Pb = Ps/k;
    %     [Pb,Ps] = berawgn(EbN0dB,'pam',M);
elseif (sch == 1)   % PSK
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
%     [Pb,Ps] = berawgn(EbN0dB,'psk',M,'nondiff');
elseif (sch == 2)                % QAM
    switch M
        case 2
           Pb = qfunc(sqrt(2*EbN0));
           Ps = Pb;
        case 4
           Ps = 1-(1-2*(1-1/sqrt(M))*qfunc(sqrt(3/(M-1)*EbN0*k))).^2;
           Pb = Ps/k;
        case 8
           [Pb,Ps] = berawgn(EbN0dB,'qam',M); % need to derive eq.
        case 16
           Ps = 1-(1-2*(1-1/sqrt(M))*qfunc(sqrt(3/(M-1)*EbN0*k))).^2;
           Pb = Ps/k;
    end
% 	[Pb,Ps] = berawgn(EbN0dB,'qam',M);
else
           [Pb,Ps] = berawgn(EbN0dB,'fsk',M,'coherent');
           % upper bounds
%            Ps_u = exp(-k*(EbN0-2*log(2))/2); % (loose)
           Ps_u = (M-1)*qfunc(sqrt(k*EbN0)); % coherent (tight)
           Pb_u = 2^(k-1)/(2^k-1)*Ps_u;%(M/2)/(M-1)*Ps_u;
end


% PLOTS [1]
if(plots)
    if (sch~=3)
        figure(1);
        scatter(real(r), imag(r)); % final SNR
        axis([-2 2 -2 2]);
    end


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
