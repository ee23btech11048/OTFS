function H_parity=select_LDPC_parity_matrix(LDPC_rate,LDPC_codeword_length)
%% generating H parity matrix;
if (LDPC_rate==1/2) && (LDPC_codeword_length==672)
    H_parity=load('H_336_672.mat');
elseif (LDPC_rate==1/2) && (LDPC_codeword_length==3840)
    H_parity=load('H_1920_3840.mat');
else
    msg = 'Choose available half rate codeword length-(672 / 3840)';
    error(msg)
end
H_parity=H_parity.H_parity;
end