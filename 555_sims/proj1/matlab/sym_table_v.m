function [sym_tbl] = sym_table_v(M,sch)

    if(sch == 0)
        fprintf('%3.1d-ary PAM\n', M);      % M-ary PAM
        switch M
            case 2,  v = [-1, 1];
            case 4,  v = [-3, -1, 1, 3];
            case 8,  v = [-7, -5, -3, -1, 1, 3, 5, 7];
            case 16, v = [-15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15];
        end % switch
    elseif (sch == 1)
        fprintf('%3.1d-ary PSK\n', M);      % M-ary PSK
        a = 1/sqrt(2);
        switch M
            case 2,  v = [-1, 1];
            case 4,  v = [1, 0, -1,  0; ... 
                          0, 1,  0, -1];
            case 8,  v = [1, a, 0, -a, -1, -a,  0,  a; ... % NOT ACCURATE
                          0, a, 1,  a,  0, -a, -1, -a];
%             case 16, v = [cos(0*pi/8)+1j*sin(0*pi/8), cos(1*pi/8)+1j*sin(1*pi/8), cos(2*pi/8)+1j*sin(2*pi/8), cos(3*pi/8)+1j*sin(3*pi/8), ...
%                            cos(4*pi/8)+1j*sin(4*pi/8), cos(5*pi/8)+1j*sin(5*pi/8), cos(6*pi/8)+1j*sin(6*pi/8), cos(7*pi/8)+1j*sin(7*pi/8), ...
%                            cos(8*pi/8)+1j*sin(8*pi/8), cos(9*pi/8)+1j*sin(9*pi/8), cos(10*pi/8)+1j*sin(10*pi/8), cos(11*pi/8)+1j*sin(11*pi/8), ...
%                            cos(12*pi/8)+1j*sin(12*pi/8), cos(13*pi/8)+1j*sin(13*pi/8), cos(14*pi/8)+1j*sin(14*pi/8), cos(15*pi/8)+1j*sin(15*pi/8)];
        end % switch
    elseif (sch == 2)
        fprintf('%3.1d-ary QAM\n', M);      % M-ary QAM
        switch M
            case 2,  v = [-1, 1];
            case 4,  v = [1, 0, -1,  0; ... 
                          0, 1,  0, -1];
            case 8,  v = [3, 1, -1, -3, -3, -1,  1,  3; ...
                          1, 1,  1,  1, -1, -1, -1, -1];
            case 16, v = [-3, -1, 1, 3, 3, 1, -1, -3, -3, -1,  1,  3,  3,  1, -1, -3;
                           3,  3, 3, 3, 1, 1,  1,  1, -1, -1, -1, -1, -3, -3, -3, -3];
        end % switch
    else                                    % M-ary Orthogonal
        fprintf('%3.1d-ary Orthogonal\n', M);
        switch M
            case 2,  v = [-1, 1];
            case 4,  v = 1/sqrt(2^log2(M))*hadamard(M);
            case 8,  v = [3+1j, 1+1j, -1+1j, -3+1j, -3-1j, -1-1j, 1-1j, 3-1j];
            case 16, v = [-3+3j, -1+3j,  1+3j,  3+3j, 3+1j,  1+1j, -1+1j, -3+1j, -3-1j, -1-1j,  1-1j,  3-1j, 3-3j,  1-3j, -1-3j, -3-3j];
        end % switch
    end % if
    
%     vs = abs(v).^2;
%     if (M~=2&&sch~=0), vs = sum(vs); end
    if (M==2||sch==0), vs = abs(v).^2;
    else vs = sum(v.^2); end
    rms = sqrt(mean(vs));
    sym_tbl = v/rms;
                
end % fcn
