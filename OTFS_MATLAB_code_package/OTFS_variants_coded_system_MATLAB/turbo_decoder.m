function [data,par_check]= turbo_decoder(x_data,M_mod,N_bits_perfram,N_syms_perfram,LDPC_trans_blocks,LDPC_codeword_length,hDec_coded_hard,RND_Perm,rev_RND_Perm)
    M_bits=log2(M_mod);
    data_info_est_coded=qamdemod(x_data,M_mod,'gray','OutputType','llr');
    data_info_est_coded_perm = data_info_est_coded(rev_RND_Perm).';
    est_coded_bits=zeros(N_bits_perfram,1);
    par_check=0;
    for i_block=1:1:LDPC_trans_blocks
        t = (i_block-1)*LDPC_codeword_length + 1;
        [est_coded_bits(t:t+LDPC_codeword_length-1,1),parity]= step(hDec_coded_hard, data_info_est_coded_perm(t:t+LDPC_codeword_length-1,1));
        par_check=par_check+sum(parity);
    end
    est_coded_bits=est_coded_bits(RND_Perm);
    data=qammod(reshape(est_coded_bits.',M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');  
end