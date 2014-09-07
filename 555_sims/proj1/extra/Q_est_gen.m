function [ out ] = Q_est_gen( x, rnd, K )
% An single estimation generator of the Q(.) function.
% theory: Q(x) ~ (# of realizations exceeding x)/K
% x = single input value
% rnd = is the randomly generated vector of N.Guassian realizations
% K = # of instances

    count = 0;

    for n = 1:K
        if (rnd(n) > x)
            count = count + 1;
        end
    end
    
    out = count/K;

end

