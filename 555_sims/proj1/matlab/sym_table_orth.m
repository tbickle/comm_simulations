function [sym_tbl] = sym_table_orth(M)
    % M-ary Orthogonal
    fprintf('%1.1d-ary Orthogonal\n', M);
    sym_tbl = 1/sqrt(2^log2(M))*hadamard(M);
end % fcn