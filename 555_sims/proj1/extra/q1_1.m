close all; clear; clc;

%% #1.1

% initial values
    K = 10^5; % 10^2, 10^3, 10^4, 10^5 % # of instances (memory limitation)
    rnd = randn(1,K); % randomly generated realizations of N.Gaussian
    dt = 1/100; temp =5;
    x = -temp:dt:(temp-dt);

% generate Q_est(x) using q_est_gen()
    Q_est = zeros(1,length(x));
    for a = 1:length(x)
        Q_est(a) = Q_est_gen(x(a), rnd, K); % my qfunc estimator
    end

% plot Q(x) & Q_est(x)
    figure(1);
    semilogy(x,qfunc(x),x,Q_est); grid on;

% %% #1.2
% 
%     K = 10^6; % # of instances (memory limitation)
% 
% % run BER trials
%     % select transmission and channel parameters
%     EbN0dBv=0:10;
%     EbN0v=10.^(EbN0dBv/10);
%     N0=1;
%     Ebv=EbN0v*N0;
%     Esv=Ebv;
%     
%     % generate BER vector
%     BER = zeros(1,length(Esv));
% %     for a = 1:length(Esv)
% %         BER(a) = BER_gen(Esv(a), N0, K);
% %     end
%     a = 1;
%     BER(a) = BER_gen(Esv(a), N0, K);
%     while (BER(a) > 100/K)
%         a = a+1;
%         BER(a) = BER_gen(Esv(a), N0, K);
%     end
%     
% % generate actual curve
%     p = qfunc(sqrt(2*EbN0v));
%     
% % plot results
%     figure(2);
%     semilogy(EbN0dBv,BER,'-o',EbN0dBv,p); %grid on;
    
    
    