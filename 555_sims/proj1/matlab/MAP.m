function [s_est] = MAP(r,M)
    % MAP detection estimation
    
    temp = zeros(log2(M),1);

    switch M
        case 2
            temp = r-1;
        case 4
            % [lsb,msb]
            switch r
                case 1,  temp(:,1) = [0;0];
                case 2,  temp(:,1) = [1;0];
                case 3,  temp(:,1) = [1;1];
                case 4,  temp(:,1) = [0;1];
                otherwise, temp(:,1) = [0;0];
            end
        case 8
            % [lsb;x;msb]
            switch r
                case 1,  temp(:,1) = [0;0;0];
                case 2,  temp(:,1) = [1;0;0];
                case 3,  temp(:,1) = [1;0;1];
                case 4,  temp(:,1) = [0;0;1];
                case 5,  temp(:,1) = [0;1;1];
                case 6,  temp(:,1) = [1;1;1];
                case 7,  temp(:,1) = [1;1;0];
                case 8,  temp(:,1) = [0;1;0];
                otherwise, temp(:,1) = [0;0;0];
            end
        case 16
            % [lsb,x,x,msb]
            switch r
                case 1,  temp(:,1) = [0;0;0;0];
                case 2,  temp(:,1) = [0;0;1;0];
                case 3,  temp(:,1) = [0;0;1;1];
                case 4,  temp(:,1) = [0;0;0;1];
                case 5,  temp(:,1) = [1;0;0;1];
                case 6,  temp(:,1) = [1;0;1;1];
                case 7,  temp(:,1) = [1;0;1;0];
                case 8,  temp(:,1) = [1;0;0;0];
                case 9,  temp(:,1) = [1;1;0;0];
                case 10, temp(:,1) = [1;1;1;0];
                case 11, temp(:,1) = [1;1;1;1];
                case 12, temp(:,1) = [1;1;0;1];
                case 13, temp(:,1) = [0;1;0;1];
                case 14, temp(:,1) = [0;1;1;1];
                case 15, temp(:,1) = [0;1;1;0];
                case 16, temp(:,1) = [0;1;0;0];
                otherwise, temp(:,1) = [0;0;0;0];
            end
    end

    s_est = temp;

end

