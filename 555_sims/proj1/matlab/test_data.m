function [b] = test_data(M)
    % test bit streams (gray coded)
    switch M
        case 4
            % QPSK test: [lsb,msb]
            b = [0,0, 1,0, 1,1, 0,1];
        case 8
            % 8PSK test: [lsb,x,msb]
            b = [0,0,0, 1,0,0, 1,0,1, 0,0,1, 0,1,1, 1,1,1, 1,1,0, 0,1,0];
%             b = [1,1,0, 1,1,0, 1,1,0, 1,1,0, 1,1,0, 1,1,0, 1,1,0, 1,1,0]; % test
        case 16
            % 16PSK test: [lsb,x,x,msb]
            b = [0,0,0,0, 0,0,1,0, 0,0,1,1, 0,0,0,1, ...
                 1,0,0,1, 1,0,1,1, 1,0,1,0, 1,0,0,0, ...
                 1,1,0,0, 1,1,1,0, 1,1,1,1, 1,1,0,1, ...
                 0,1,0,1, 0,1,1,1, 0,1,1,0, 0,1,0,0];
%             b = [0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, ...
%                  0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0];
        otherwise
            b = [];
    end
end

