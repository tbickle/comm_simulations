function [y] = gray_code(x,M)
    % cross reference for gray code

    temp = zeros(1,length(x));
    switch M
        case 2
            temp = x;
        case 4
            for i = 1:length(x)
                switch x(i)
                    case 4,  temp(i) = 3;
                    case 3,  temp(i) = 4;
                    otherwise, temp(i) = x(i);
                end
            end
        case 8
            for i = 1:length(x)
                switch x(i)
                    case 6,  temp(i) = 3;
                    case 5,  temp(i) = 4;
                    case 7,  temp(i) = 5;
                    case 8,  temp(i) = 6;
                    case 4,  temp(i) = 7;
                    case 3,  temp(i) = 8;
                    otherwise, temp(i) = x(i);
                end
            end
        case 16
            for i = 1:length(x)
                switch x(i)
                    case 1,  temp(i) = 1;
                    case 5,  temp(i) = 2;
                    case 13, temp(i) = 3;
                    case 9,  temp(i) = 4;
                    case 10, temp(i) = 5;
                    case 14, temp(i) = 6;
                    case 6,  temp(i) = 7;
                    case 2,  temp(i) = 8;
                    case 4,  temp(i) = 9;
                    case 8,  temp(i) = 10;
                    case 16, temp(i) = 11;
                    case 12, temp(i) = 12;
                    case 11, temp(i) = 13;
                    case 15, temp(i) = 14;
                    case 7,  temp(i) = 15;
                    case 3,  temp(i) = 16;
                    otherwise, temp(i) = x(i);
                end
            end
    end
    
    y = temp;

end
