close all; clear; clc;

%% #1.1

% initial values
    K = 10^4; % 10^2, 10^3, 10^4, 10^5 % # of instances (memory limitation)
    rnd = randn(1,K); % randomly generated realizations of N.Gaussian
    dt = 1/100; temp = 3;
    x = -temp:dt:(temp-dt);

% generate Q_est(x) using q_est_gen()
    q_est = zeros(1,length(x));
    for a = 1:length(x)
        q_est(a) = q_est_gen(x(a), rnd, K); % my qfunc estimator
    end

% plot Q(x) & Q_est(x)
    figure(1);
    semilogy(x,qfunc(x),x,q_est);

%% #1.2

% generate random bit stream b
    b = zeros(1,K); % random bit vector    
    for a = 1:K
        % i'm declaring Es = 1 -> b = 1 <------------------------
        % and Es = -1 -> b = 0      <--------------------
        if (rnd(a) > 0)
            b(a) = 1;
        else
            b(a) = 0;
        end
    end


% run BER trials
    % select transmission and channel parameters
%     Es = 1; % arbitrarily selected value <---------------
    Es = 0:0.5:10; % arbitrarily selected values <---------------
    Eb = Es;
    N0 = 1; % arbitrarily selected value <---------------
    
    % generate BER vector
    BER = zeros(1,length(Es));
    for a = 1:length(Es)
        BER(a) = BER_gen(Es(a), N0, b, K);
    end
    
% generate actual curve
    p = qfunc(sqrt(2*Eb/N0));
    
% plot results
    figure(2);
    semilogy(Eb/N0,BER,Eb/N0,p);
    
    
    