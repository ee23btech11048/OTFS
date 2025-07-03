%References
%  [R1]. Y. Hong, T. Thaj, E. Viterbo, ``Delay-Doppler Communications: Principles and Applications'', Academic Press, 2022, ISBN:9780323850285



function [est_bits,ite,x_data] = MRC_turbo_decoder_OTFSvariants(N,M,M_mod,no,data_grid,Y_tilda,H_tf,n_ite_MRC,omega,r,Fn,L_set,nu_ml_tilda,init_estimate,LDPC_rate,LDPC_codeword_length,RND_Perm,rev_RND_Perm,hDec_coded_hard,hDec,variant)
%% Initial assignments
%Number of symbols per frame
N_syms_perfram=sum(sum((data_grid>0)));
%Arranging the delay-Doppler grid symbols into an array
data_array=reshape(data_grid,1,N*M);
%finding position of data symbols in the array
[~,data_index]=find(data_array>0);
M_bits=log2(M_mod);
N_bits_perfram = N_syms_perfram*M_bits;
L_set=L_set+1;
error=zeros(n_ite_MRC);

%% LDPC code parameters
LDPC_info_length = LDPC_codeword_length*LDPC_rate;
LDPC_trans_blocks = floor(N_bits_perfram/LDPC_codeword_length);

%% OTFS variant parameters  (in this code we combine the MRC detector for all variants (RCP/RZP/ZP/CP) into a single code using these parameters)
% Check Section 4.5 in Chapter 4, [R1] to see the differences between
% variants
ZP_flag=0;              
cshift=0;
if(strcmp(variant,'ZP'))
    ZP_flag=1;
elseif(strcmp(variant,'CP'))
    cshift=0;
else
    cshift=1;
end
 

%% initial estimate using single tap TF equalizer
if(init_estimate==1)
    Y_tf=fft(Y_tilda).'; % delay-time to frequency-time domain                                                                      % Section 6.5.5 in Chapter 6, [R1]  
    X_tf=conj(H_tf).*Y_tf./(H_tf.*conj(H_tf)+no); % single tap equalizer                                                            
    X_est = ifft(X_tf.')*Fn; % SFFT                                                                                                 
    X_est=qammod(qamdemod(X_est,M_mod,'gray'),M_mod,'gray');
    X_est=X_est.*data_grid;
    X_tilda_est=X_est*Fn';
else
    X_est=zeros(M,N);
    X_tilda_est=X_est*Fn';
end
x_m=X_est.';
x_m_tilda=X_tilda_est.';


%% MRC detector     Algotithm 3 in Chapter 6, [R2])
%% initial computation
d_m_tilda=zeros(N,M);
y_m_tilda=reshape(r,M,N).';
delta_y_m_tilda=y_m_tilda;
for m=1:M
    for l=L_set
        if(m+l-1<=M)
            d_m_tilda(:,m)=d_m_tilda(:,m)+abs(nu_ml_tilda(:,m+(l-1),l).^2);        % equation (6.49) in Chapter 6, [R1]
        elseif(not(ZP_flag))
            d_m_tilda(:,m)=d_m_tilda(:,m)+circshift(abs(nu_ml_tilda(:,mod(m+(l-1)-1,M)+1,l).^2),-cshift);  % See Section 4.5, Chapter 4 and 6.5.6, Chapter 6 in [R1]
        end
    end
end
for m=1:M
    for l=L_set
        if(m>=l)
            delta_y_m_tilda(:,m)=delta_y_m_tilda(:,m)-nu_ml_tilda(:,m,l).*x_m_tilda(:,m-(l-1));  % Line 6 of Algorithm 3 in Chapter 6, [R1]
        elseif(not(ZP_flag))
            delta_y_m_tilda(:,m)=delta_y_m_tilda(:,m)-nu_ml_tilda(:,m,l).*circshift(x_m_tilda(:,mod(m-(l-1)-1,M)+1),cshift); % See Section 4.5, Chapter 4 and 6.5.6, Chapter 6 in [R1]
        end
    end
end
x_m_tilda_old=x_m_tilda;
c_m_tilda=x_m_tilda;

%% iterative computation
for ite=1:n_ite_MRC
    %% MRC detector (one iteration)
    
    delta_g_m_tilda=zeros(N,M);
    for m=1:M
         for l=L_set
            if(m+(l-1)<=M)
                delta_g_m_tilda(:,m)=delta_g_m_tilda(:,m)+conj(nu_ml_tilda(:,mod(m+(l-1)-1,M)+1,l)).*delta_y_m_tilda(:,mod(m+(l-1)-1,M)+1);  % Line 9 of Algorithm 3 in Chapter 6, in [R1]
            elseif(not(ZP_flag))
                delta_g_m_tilda(:,m)=delta_g_m_tilda(:,m)+circshift(conj(nu_ml_tilda(:,mod(m+(l-1)-1,M)+1,l)).*delta_y_m_tilda(:,mod(m+(l-1)-1,M)+1),-cshift); % See Section 4.5, Chapter 4 and 6.5.6, Chapter 6 in [R1]
            end
        end
        c_m_tilda(:,m)=x_m_tilda_old(:,m)+delta_g_m_tilda(:,m)./d_m_tilda(:,m);                                             % Line 10 of Algorithm 3 in Chapter 6, in [R1]
        x_m_tilda(:,m)=c_m_tilda(:,m);
        for l=L_set                                                                                                         % Line 12 of Algorithm 3 in Chapter 12, in [R1]
            if(m+(l-1)<=M)
                delta_y_m_tilda(:,m+(l-1))=delta_y_m_tilda(:,m+(l-1))-nu_ml_tilda(:,m+(l-1),l).*(x_m_tilda(:,m)-x_m_tilda_old(:,m));     % Line 13 of Algorithm 3 in Chapter 12, in [R1]
            elseif(not(ZP_flag))
                delta_y_m_tilda(:,mod(m+(l-1)-1,M)+1)=delta_y_m_tilda(:,mod(m+(l-1)-1,M)+1)-nu_ml_tilda(:,mod(m+(l-1)-1,M)+1,l).*circshift((x_m_tilda(:,m)-x_m_tilda_old(:,m)),cshift); % See Section 4.5, Chapter 4 and 6.5.6, Chapter 6 in [R1]
            end                                                                     
        end                                                                                                                 % Line 14 of Algorithm 3 in Chapter 12, in [R1]
        x_m_tilda_old(:,m)=x_m_tilda(:,m);
    end
    
    %% turbo decoder   (see Section 6.6)
    
    X_est=(Fn*x_m_tilda).';
    x_est=reshape(X_est,1,N*M);
    x_data_old=x_est(data_index);    
    [x_data_new,par_check]= turbo_decoder(x_data_old,M_mod,N_bits_perfram,N_syms_perfram,LDPC_trans_blocks,LDPC_codeword_length,hDec_coded_hard,RND_Perm,rev_RND_Perm);    
    X = Generate_2D_data_grid(N,M,x_data_new,data_grid);
    x_m_tilda=(1-omega)*x_m_tilda_old+omega*Fn'*X.';
    if(par_check==0) %% par_check=0 implies all the LDPC codewords were decoded correctly.
        break;
    end
    
    %% updating the residual error vectors
    
    for m=1:M
        for l=L_set
            if(m+(l-1)<=M)
                delta_y_m_tilda(:,m+(l-1))=delta_y_m_tilda(:,m+(l-1))-nu_ml_tilda(:,m+(l-1),l).*(x_m_tilda(:,m)-x_m_tilda_old(:,m));
            elseif(not(ZP_flag))
                delta_y_m_tilda(:,mod(m+(l-1)-1,M)+1)=delta_y_m_tilda(:,mod(m+(l-1)-1,M)+1)-nu_ml_tilda(:,mod(m+(l-1)-1,M)+1,l).*circshift((x_m_tilda(:,m)-x_m_tilda_old(:,m)),1);
            end
        end
        
    end
    x_m_tilda_old=x_m_tilda;
    
    %% convergence criteria
    error(ite)=norm(delta_y_m_tilda);
    if(ite>1)
        if(error(ite)>=error(ite-1))
            break;
        end
    end
end
if(n_ite_MRC==0)
    ite=0;
end
%% detector output bits
X_est=(Fn*x_m_tilda).';
x_est=reshape(X_est,1,N*M);
x_data=x_est(data_index);
data_info_est_coded=qamdemod(x_data,M_mod,'gray','OutputType','llr');
data_info_est_coded_perm = data_info_est_coded(rev_RND_Perm).';
for i_block=1:1:LDPC_trans_blocks
    t1 = (i_block-1)*LDPC_info_length + 1;
    t2 = (i_block-1)*LDPC_codeword_length + 1;
    est_bits(t1:t1+LDPC_info_length-1,1) = step(hDec, data_info_est_coded_perm(t2:t2+LDPC_codeword_length-1,1));  % LDPC decoder to get the estimated information bits
end
end
