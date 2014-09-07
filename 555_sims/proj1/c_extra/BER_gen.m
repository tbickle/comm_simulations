function [ BER ] = BER_gen( Es, N0, b, K )
% Generates Bit Error Rate
% Es = symbol energy
% N0 = noise energy
% b = transmitted bit stream
% K = # of instances

    s = zeros(1,length(b));
    for a = 1:K
        % i'm declaring Es = 1 -> b = 1 <------------------------
        % and Es = -1 -> b = 0      <--------------------
        if (b(a) == 1)
            s(a) = sqrt(Es);
        else
            s(a) = -sqrt(Es);
        end
    end

% generate Gaussian vector of variance N0/2
    n = 0 + sqrt(N0/2)*randn(1,K); % N.Gaussian vector w/ variance N0/2 and 0 mean

% generate observation equation
    z = s + n; % observation equation

% MAP detection to determine b_est
    b_est = zeros(1,length(z));
    for a = 1:length(b_est)
        if (z(a) > 0)
            b_est(a) = 1;
        end
    end

% count differences b/w b & b_est to determine BER
    temp = 0;
    for a = 1:length(b)
        if (b(a) ~= b_est(a))
            temp = temp + 1;
        end
    end
    BER = temp/K;


end

