function [sym_tbl] = sym_table(M,sch)

    if(sch == 0)
        fprintf('%1.1d-ary PAM\n', M);      % M-ary PAM
        switch M
            case 2,  v = [-1, 1]; % 2PAM symbols
            case 4,  v = [-3, -1, 1, 3]; % 4PAM symbols
            case 8,  v = [-7, -5, -3, -1, 1, 3, 5, 7]; % 8PAM symbols
            case 16, v = [-15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15]; % 16PAM symbols
        end % switch
    elseif (sch == 1)
        fprintf('%1.1d-ary PSK\n', M);      % M-ary PSK
        switch M
            case 2,  v = [-1, 1]; % BPSK symbols
            case 4,  v = [1, 1j, -1, -1j]; % QPSK symbols
            case 8,  v = [1, (1+1j)/sqrt(2), 1j, (-1+1j)/sqrt(2), -1, (-1-1j)/sqrt(2), -1j, (1-1j)/sqrt(2)]; % 8PSK symbols
            case 16, v = [cos(0*pi/8)+1j*sin(0*pi/8), cos(1*pi/8)+1j*sin(1*pi/8), cos(2*pi/8)+1j*sin(2*pi/8), cos(3*pi/8)+1j*sin(3*pi/8), ...
                           cos(4*pi/8)+1j*sin(4*pi/8), cos(5*pi/8)+1j*sin(5*pi/8), cos(6*pi/8)+1j*sin(6*pi/8), cos(7*pi/8)+1j*sin(7*pi/8), ...
                           cos(8*pi/8)+1j*sin(8*pi/8), cos(9*pi/8)+1j*sin(9*pi/8), cos(10*pi/8)+1j*sin(10*pi/8), cos(11*pi/8)+1j*sin(11*pi/8), ...
                           cos(12*pi/8)+1j*sin(12*pi/8), cos(13*pi/8)+1j*sin(13*pi/8), cos(14*pi/8)+1j*sin(14*pi/8), cos(15*pi/8)+1j*sin(15*pi/8)]; % 16PSK symbols
        end % switch
    elseif (sch == 2)                       % M-ary QAM
        fprintf('%1.1d-ary QAM\n', M);
        switch M
            case 2,  v = [-1, 1]; % 2QAM symbols
            case 4,  v = [1, 1j, -1, -1j]; % 4QAM symbols
            case 8,  v = [3+1j, 1+1j, -1+1j, -3+1j, -3-1j, -1-1j, 1-1j, 3-1j]; % 8QAM symbols
            case 16, v = [-3+3j, -1+3j,  1+3j,  3+3j, 3+1j,  1+1j, -1+1j, -3+1j, -3-1j, -1-1j,  1-1j,  3-1j, 3-3j,  1-3j, -1-3j, -3-3j]; % 16QAM symbols
        end % switch
    else                                    % M-ary Orthogonal    
        fprintf('%1.1d-ary Orthogonal\n', M);
        v = 1/sqrt(2^log2(M))*hadamard(M);
    end % if
    
    if (sch==3), sym_tbl = v;
    else
        rms = sqrt(mean(abs(v).^2));
        sym_tbl = v/rms;
        clear v;
    end % if
                
end % fcn
