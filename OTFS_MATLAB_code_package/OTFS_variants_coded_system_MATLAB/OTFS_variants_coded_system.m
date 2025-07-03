%References
%  [R1]. Y. Hong, T. Thaj, E. Viterbo, ``Delay-Doppler Communications: Principles and Applications'', Academic Press, 2022, ISBN:9780323850285
%  [R2]. T. Thaj and E. Viterbo, "Low Complexity Iterative Rake Decision Feedback Equalizer for Zero-Padded OTFS Systems," in IEEE Transactions on Vehicular Technology, vol. 69, no. 12, pp. 15606-15622, Dec. 2020, doi: 10.1109/TVT.2020.3044276.





close all
clear all
rng('shuffle')

%% OTFS variant: (RZP / RCP / CP / ZP)
variant='RCP';

%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 64;
% M: number of subcarriers in frequency
M = 64;
% M_mod: size of QAM constellation
M_mod = 64;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));

%% delay-Doppler grid symbol placement
if(strcmp(variant,'ZP'))         
    length_ZP = M/16; % ZP length (required only for ZP-OTFS)
    length_CP = 0;
elseif(strcmp(variant,'CP'))
    length_ZP = 0;
    length_CP = M/16; % CP length (required only for CP-OTFS) 
elseif(strcmp(variant,'RCP')||strcmp(variant,'RZP'))
    length_ZP=0;
    length_CP=0;
else
    msg = 'Choose a valid OTFS variant: (RZP / RCP / CP / ZP)';
    error(msg)
end
M_data=M-length_ZP;
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
data_grid=zeros(M,N);
data_grid(1:M_data,1:N)=1;
% number of symbols per frame
N_syms_perfram = sum(sum(data_grid));
% number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;


% Time and frequency resources
car_fre = 4*10^9;% Carrier frequency
delta_f = 15*10^3; % subcarrier spacing: 15 KHz
T = 1/delta_f; %one time symbol duration in OTFS frame

% SNR and variance of the noise
% SNR = P/\sigma^2; P: avg. power of albhabet transmitted
SNR_dB = 15:2.5:25;
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;

%% Error-correcting code parameters
% LDPC code parameters
LDPC_rate = 1/2;% can be changed to 3/4
LDPC_codeword_length = 672; %(672 / 3840)
LDPC_info_length = LDPC_codeword_length*LDPC_rate;
% LDPC_trans_blocks: number of LDPC blocks transmitted in each frame
% we transmit zero's for the remaining bits in N_bits_perfram
LDPC_trans_blocks = floor(N_bits_perfram/LDPC_codeword_length);
trans_bits_length = LDPC_trans_blocks*LDPC_codeword_length;
trans_info_bits_length = LDPC_trans_blocks*LDPC_info_length;
trans_symbols_tot = trans_info_bits_length/M_bits;
[hEnc,hDec,hDec_coded_soft,hDec_coded_hard]=LDPC_system_objects(LDPC_rate,LDPC_codeword_length);


%% Initializing simulation error count variables

N_fram = 1000;
est_info_bits_MRC=zeros(N_bits_perfram,1);
err_ber_MRC = zeros(1,length(SNR_dB));
avg_ber_MRC = zeros(1,length(SNR_dB));
err_fer_MRC = zeros(1,length(SNR_dB));
avg_fer_MRC = zeros(1,length(SNR_dB));
no_of_detetor_iterations_MRC= zeros(length(SNR_dB),1);
avg_no_of_iterations_MRC=zeros(1,length(SNR_dB));


%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix
current_frame_number=zeros(1,length(SNR_dB));


for iesn0 = 1:length(SNR_dB)
    for ifram = 1:N_fram
        current_frame_number(iesn0)=ifram;
        %% random input bits generation%%%%%                      
        trans_info_bits = randi([0,1],trans_info_bits_length,1);
        %% LDPC encoder %%%%%%
        coded_trans_info_bit = zeros(N_bits_perfram,1);% bits other trans_bits_length are considered as zeros
        for i_block=1:1:LDPC_trans_blocks
            t1 = (i_block-1)*LDPC_info_length + 1;
            t2 = (i_block-1)*LDPC_codeword_length + 1;
            coded_trans_info_bit(t2:t2+(LDPC_codeword_length)-1,1) = step(hEnc, trans_info_bits(t1:t1+LDPC_info_length-1));
        end
       
       %%%% Random Permutation of bits across the frame
        state=randi(10000);
        RND_Perm=randintrlv(1:N_bits_perfram,state);
        rev_RND_Perm=randdeintrlv(1:N_bits_perfram,state);
        coded_trans_info_bit_perm = coded_trans_info_bit(RND_Perm); 
        data_info_est_coded=coded_trans_info_bit_perm.';
        
        %%2D QAM symbols generation %%%%%%%%
        data=qammod(reshape(coded_trans_info_bit_perm,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');        
        X = Generate_2D_data_grid(N,M,data,data_grid);
        
        
        %% OTFS modulation%%%%
        X_tilda=X*Fn';               
        s = reshape(X_tilda,N*M,1); 
               
        
   %% OTFS channel generation%%%% The user can use either the synthetic channel model or the 3GPP channel by uncommenting the corresonding piece of code.        
        %% synthetic channel model with equal power paths with delays [0,l_max] and Dopplers [-k_max,k_max]
%                 taps=4;
%                 l_max=delay_spread;
%                 k_max=4;
%                 chan_coef=1/sqrt(2)*(randn(1,taps)+1i.*randn(1,taps));
%                 delay_taps=randi(taps,[1,l_max+1]);  
%                 delay_taps=sort(delay_taps-min(delay_taps));  %% random delay shifts in the range [0,l_max] 
%                 Doppler_taps=k_max-2*k_max*rand(1,taps);   %% uniform Doppler profile [-k_max,k_max]
        
        % channel model following 3GPP standard
        max_speed=500;  % km/hr
        [chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed);
        L_set=unique(delay_taps);
        l_max=max(L_set);
        
        %% channel output%%%%%                      
        gs=Gen_time_domain_channel_OTFSvariants(N,M,delay_taps,Doppler_taps,chan_coef,length_CP,variant);
        r=Gen_channel_output_OTFSvariants(N,M,L_set,gs,s,sigma_2(iesn0),length_CP,variant);              
        
        %% OTFS demodulation%%%%
        Y_tilda=reshape(r,M,N);     
        Y = Y_tilda*Fn;             
             
        %% delay-time channel vectors
        [nu_ml_tilda]=Gen_delay_time_channel_vectors_OTFSvariants(N,M,l_max,gs,length_CP,variant);
        
        %% Generate block-wise time-frequency domain channel
        [H_t_f]=Generate_time_frequency_channel_OTFSvariants(N,M,gs,L_set,length_CP,variant);
        
        %% MRC turbo detection in [R1,R2,R3]
        
        n_ite_MRC=50; % maximum number of MRC detector iterations        
        omega=1; %damping parameter - optimizing omega improves error performance
        init_estimate=0;    %1-use the TF single tap estimate as the initial estimate for MRC detection, 0-initialize the symbol estimates to 0 at the start of MRC iteration 
        %(Note: it is recommended to set init_estimate to 0 for higher order modulation schemes like 64-QAM or 256-QAM as the single tap equalizer estimate may be less accurate)         
        
        [est_info_bits_MRC,det_iters_MRC,data_MRC] = MRC_turbo_decoder_OTFSvariants(N,M,M_mod,sigma_2(iesn0),data_grid,Y_tilda,H_t_f,n_ite_MRC,omega,r,Fn,L_set,nu_ml_tilda,init_estimate,LDPC_rate,LDPC_codeword_length,RND_Perm,rev_RND_Perm,hDec_coded_hard,hDec,variant);

        %% errors count%%%%%     
        errors_MRC = sum(xor(est_info_bits_MRC,trans_info_bits));
        if(errors_MRC>0)
            err_ber_MRC(1,iesn0) = err_ber_MRC(1,iesn0) + errors_MRC;           % Bit error rate
            err_fer_MRC(1,iesn0) = err_fer_MRC(1,iesn0) + 1;                    % Frame error rate: frame error happends if even one bit in a frame is wronf
        end
        avg_ber_MRC(1,iesn0)=err_ber_MRC(1,iesn0).'/length(trans_info_bits)/ifram;
        avg_fer_MRC(1,iesn0)=err_fer_MRC(1,iesn0).'/ifram;

        
        %%  iterations count
        
        no_of_detetor_iterations_MRC(iesn0)=no_of_detetor_iterations_MRC(iesn0)+det_iters_MRC;
        avg_no_of_iterations_MRC(iesn0)=no_of_detetor_iterations_MRC(iesn0)/ifram;
        
        
        %%         DISP error performance details
        clc
        fprintf('%s%s',variant,'-OTFS')
        fprintf('(N,M,QAM size,LDPC codeword length)');disp([N,M,M_mod,LDPC_codeword_length]);
        display(current_frame_number,'Number of frames');
        display(avg_ber_MRC,'Average coded BER - Delay-time domain Maximal Ratio Combining (MRC)');
        display(avg_fer_MRC,'Average coded FER - Delay-time domain Maximal Ratio Combining (MRC)');
        display(avg_no_of_iterations_MRC,'Average number of iterations for the MRC turbo decoder');
        
        
        
    end
    
end

figure(1)
semilogy(SNR_dB,avg_fer_MRC,'-x','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,avg_ber_MRC,'-x','LineWidth',2,'MarkerSize',8)
legend('coded FER','coded BER')
grid on
xlabel('SNR(dB)')
ylabel('BER/FER')