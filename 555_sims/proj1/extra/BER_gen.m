function [ BER ] = BER_gen( Es, N0, K )
% Generates Bit Error Rate
% Es = symbol energy
% N0 = noise energy
% K = # of instances

% % generate random bit and signal vectors
    temp = sign(randn(1,K));
    b = (temp+1)/2;
    s = sqrt(Es)*temp;
% 
% generate Gaussian vector of variance N0/2
    n = 0 + sqrt(N0/2)*randn(1,K); % N.Gaussian vector w/ variance N0/2 and 0 mean

% generate observation equation
    z = s + n; % observation equation

% MAP detection to determine b_est
    b_est = (sign(z)+1)/2;

% count differences b/w b & b_est to determine BER
    BER = sum(abs(b - b_est))/(K*1.0);

end