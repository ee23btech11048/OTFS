%References
%  [R1]. Y. Hong, T. Thaj, E. Viterbo, ``Delay-Doppler Communications: Principles and Applications'', Academic Press, 2022, ISBN:9780323850285

function [nu_ml_tilda]=Gen_delay_time_channel_vectors_OTFSvariants(N,M,l_max,gs,L_cp,variant)
nu_ml_tilda=zeros(N,M,l_max+1);
for n=1:N
    for m=1:M
        for ell=1:(l_max+1)
            if(strcmp(variant,'RZP'))   %Section 4.5.1, Chapter 4, [R1]
                if(not(m<ell && n==1))
                    nu_ml_tilda(n,m,ell)=gs(ell,m+(n-1)*M);
                end
                
            elseif(strcmp(variant,'RCP'))  %Section 4.5.2, Chapter 4, [R1]
                nu_ml_tilda(n,m,ell)=gs(ell,m+(n-1)*M);
                
            elseif(strcmp(variant,'CP'))  %Section 4.5.3, Chapter 4, [R1]  (type(b) in Fig 4.13)
                nu_ml_tilda(n,m,ell)=gs(ell,m+(n-1)*(M+L_cp));
                
            elseif(strcmp(variant,'ZP'))  %Section 4.5.4, Chapter 4, [R1]  (type(c) in Fig 4.13)
                if(m>=ell)
                    nu_ml_tilda(n,m,ell)=gs(ell,m+(n-1)*M);
                end
                
            end
        end
    end
end
end