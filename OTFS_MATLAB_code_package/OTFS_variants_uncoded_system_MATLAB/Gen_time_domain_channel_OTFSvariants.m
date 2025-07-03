%References
%  [R1]. Y. Hong, T. Thaj, E. Viterbo, ``Delay-Doppler Communications: Principles and Applications'', Academic Press, 2022, ISBN:9780323850285


function [gs]=Gen_time_domain_channel_OTFSvariants(N,M,delay_taps,Doppler_taps,chan_coef,L_cp,variant)
z=exp(1i*2*pi/N/M);
P=length(delay_taps);
l_max=max(delay_taps);
if(strcmp(variant,'CP'))  
    frame_size=N*(M+L_cp);
else
    frame_size=N*M;
end
gs=zeros(l_max+1,frame_size);
for q=0:(frame_size-1)
    for i=1:P
        g_i=chan_coef(i);
        l_i=delay_taps(i);
        k_i=Doppler_taps(i);  
        gs(l_i+1,q+1)=gs(l_i+1,q+1)+g_i*z^(k_i*(q-l_i));  %(see Appendix C.4 in [R1])
    end    
end
end