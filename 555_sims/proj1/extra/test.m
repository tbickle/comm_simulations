% % Probability of a bit error for M-PAM
% % EE 161 - Digital Communication Systems - Spring 2011
% % Prof. Robert Morelos-Zaragoza. San Jose State University
% 
% clear all;
% set(1)=2; set(2)=4; set(3)=8; set(4)=16; set(5)=32;  % Constellation sizes
% for j=1:5
%     i = 1;
%     M = set(j); 
%     for ebnodb=0:1:30
%         ebno = 10^(ebnodb/10);
%         Pe(j,i) = (2*(M-1)/(log2(M)*M))*qfunc(sqrt((6*log2(M)/(M^2-1))*ebno)); 
%         snr(i) = ebnodb;
%         i=i+1;
%     end
% end
% semilogy(snr,Pe(1,:),'-^',snr,Pe(2,:),'-o',snr,Pe(3,:),'-s', ...
%                                     snr,Pe(4,:),'-*',snr,Pe(5,:),'-+')
% axis( [0 30 1e-6 1e-1]), grid on
% xlabel('E_b/N_0 (dB)'); ylabel('Bit error probability')
% legend('2-PAM','4-PAM','8-PAM','16-PAM','32-PAM',1);
% title('Average probability of a bit error for M-PAM')
% hold off

% Comparison of the bit error rates of two bits-to-symbol mappings for a
% 4-PAM constellation
% EE 161 - Digital Communication Systems - Spring 2012
% Prof. Robert Morelos-Zaragoza. San Jose State University

close all; clear all; clc;

id = input('Enter your student ID (tower card) number: ');
rand('state',id), randn('state',id)

% 4-PAM signal constellation with Natural mapping:
% Natural Mapping:
% B1B2  Amp
% 00   -3/sqrt(5)
% 01   -1/sqrt(5)
% 10    1/sqrt(5)
% 11    3/sqrt(5)
S_nat =  [-3/sqrt(5) -1/sqrt(5) 1/sqrt(5) 3/sqrt(5)];

% 4-PAM signal constellation with Gray mapping:
% Gray Mapping:
% B1B2  A
% 00   -3/sqrt(5)
% 01   -1/sqrt(5)
% 11    1/sqrt(5)
% 10    3/sqrt(5)
S_gra = [-3/sqrt(5) -1/sqrt(5) 3/sqrt(5) 1/sqrt(5)];

Nsim = 1000000;                     % Number of simulated symbols
i = 1;                              % Index for arrays with results

for EbNo=0:1:13
    
    error_n = 0; error_g = 0;       % Error counters
    
    for n = 1:Nsim
        
        No = 10^(-EbNo/10)/2;
        sigma= sqrt(No/2);          % Standard deviation of AWGN samples
        
        % Random bits
        B1n=round(rand); B2n=round(rand);       % Natural mapping  
        B1g=round(rand); B2g=round(rand);       % Gray mapping 

        % Amplitudes (symbols)
        A_nat = S_nat(2*B1n+B2n+1);             % Natural mapping              
        A_gra = S_gra(2*B1g+B2g+1);             % Gray mapping              
        
        % Add AWGN:
        Ynat = A_nat + sigma*randn;             % Natural mapping 
        Ygra = A_gra + sigma*randn;             % Gray mapping
        
        % Decision device: Natural
        if Ynat < -2/sqrt(5),                        B1nhat=0; B2nhat=0;
        elseif ((Ynat>-2/sqrt(5)) && (Ynat<=0)),     B1nhat=0; B2nhat=1;
        elseif ((Ynat>0) && (Ynat<=2/sqrt(5))),      B1nhat=1; B2nhat=0;
        else B1nhat=1; B2nhat=1;
        end

        % Decision device: Gray
        if Ygra < -2/sqrt(5),                        B1ghat=0; B2ghat=0;
        elseif ((Ygra>-2/sqrt(5)) && (Ygra<=0)),     B1ghat=0; B2ghat=1;
        elseif ((Ygra>0) && (Ygra<=2/sqrt(5))),      B1ghat=1; B2ghat=1;
        else B1ghat=1; B2ghat=0;
        end

        % Count errors
        if B1nhat ~= B1n, error_n = error_n + 1; end
        if B2nhat ~= B2n, error_n = error_n + 1; end
        if B1ghat ~= B1g, error_g = error_g + 1; end
        if B2ghat ~= B2g, error_g = error_g + 1; end

    end % Nsim
    
    snr(i) = EbNo;
    ber_n(i) = error_n/(2*Nsim);
    ber_g(i) = error_g/(2*Nsim);
    fprintf('%3.1f\t%10.7e\t%10.7e\n', snr(i), ber_n(i), ber_g(i));
    i = i + 1;
        
end % EbNo


[Pb_t,Ps_t] = berawgn(snr,'pam',4);

% Plotting commands
figure(1)
semilogy(snr,ber_n,'-^k'), hold on
semilogy(snr,ber_g,'-ok'), hold on
Pb = (3/4) * qfunc(sqrt((4/5)*10.^(snr/10)));
semilogy(snr,Pb,'-*k')
semilogy([0:9],qfunc(sqrt(2*10.^([0:9]/10))),'-sk'); hold on;
semilogy(snr,Pb_t,'ob');
legend('Simulated BER, Natural','Simulated BER, Gray','Approximation','2-PAM');
title('Digital communication using 4-PAM over an AWGN channel');
xlabel('E_b/N_0 (dB)'); ylabel('Bit error rate');
grid on, axis([3 13 1e-5 1e-1]), hold off

% [1] http://www.engr.sjsu.edu/rmorelos/ee161s12/index.html